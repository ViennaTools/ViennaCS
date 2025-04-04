#pragma once

#include "csDenseCellSet.hpp"
#include "csTracePath.hpp"
#include "csTracingParticle.hpp"

#include <rayBoundary.hpp>
#include <rayGeometry.hpp>
#include <raySource.hpp>

namespace viennacs {

using namespace viennacore;

template <typename T, int D> class TracingKernel {
public:
  TracingKernel(RTCDevice &pDevice, viennaray::Geometry<T, D> &pRTCGeometry,
                viennaray::Boundary<T, D> &pRTCBoundary,
                viennaray::Source<T> &pSource,
                std::unique_ptr<AbstractParticle<T>> &pParticle,
                const size_t pNumOfRayPerPoint, const size_t pNumOfRayFixed,
                const bool pUseRandomSeed, const size_t pRunNumber,
                SmartPointer<DenseCellSet<T, D>> passedCellSet,
                int passedExclude)
      : mDevice(pDevice), mGeometry(pRTCGeometry), mBoundary(pRTCBoundary),
        mSource(pSource), mParticle(pParticle->clone()),
        mNumRays(pNumOfRayFixed == 0
                     ? pSource.getNumPoints() * pNumOfRayPerPoint
                     : pNumOfRayFixed),
        mUseRandomSeeds(pUseRandomSeed), mRunNumber(pRunNumber),
        cellSet(passedCellSet), excludeMaterial(passedExclude),
        mGridDelta(cellSet->getGridDelta()) {
    assert(rtcGetDeviceProperty(mDevice, RTC_DEVICE_PROPERTY_VERSION) >=
               30601 &&
           "Error: The minimum version of Embree is 3.6.1");
  }

  void apply() {
    auto rtcScene = rtcNewScene(mDevice);
    rtcSetSceneFlags(rtcScene, RTC_SCENE_FLAG_NONE);
    auto bbquality = RTC_BUILD_QUALITY_HIGH;
    rtcSetSceneBuildQuality(rtcScene, bbquality);
    auto rtcGeometry = mGeometry.getRTCGeometry();
    auto rtcBoundary = mBoundary.getRTCGeometry();

    auto boundaryID = rtcAttachGeometry(rtcScene, rtcBoundary);
    auto geometryID = rtcAttachGeometry(rtcScene, rtcGeometry);
    assert(rtcGetDeviceError(mDevice) == RTC_ERROR_NONE &&
           "Embree device error");

    const auto meanFreePath = mParticle->getMeanFreePath();

    auto myCellSet = cellSet;

#pragma omp parallel shared(myCellSet)
    {
      rtcJoinCommitScene(rtcScene);

      alignas(128) auto rayHit =
          RTCRayHit{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      unsigned int seed = mRunNumber;
      if (mUseRandomSeeds) {
        std::random_device rd;
        seed = rd();
      }

      // thread-local particle object
      auto particle = mParticle->clone();

      // thread local path
      TracePath<T> path;
      // if (!traceOnPath)
      path.useGridData(myCellSet->getNumberOfCells());

#if VIENNARAY_EMBREE_VERSION < 4
      auto rtcContext = RTCIntersectContext{};
      rtcInitIntersectContext(&rtcContext);
#endif

#pragma omp for schedule(dynamic)
      for (long long idx = 0; idx < mNumRays; ++idx) {
        // particle specific RNG seed
        auto particleSeed = tea<3>(idx, seed);
        RNG RngState(particleSeed);

        particle->initNew(RngState);

        auto originAndDirection = mSource.getOriginAndDirection(idx, RngState);
        rayInternal::fillRayPosition(rayHit.ray, originAndDirection[0]);
        rayInternal::fillRayDirection(rayHit.ray, originAndDirection[1]);

#ifdef VIENNARAY_USE_RAY_MASKING
        rayHit.ray.mask = -1;
#endif

        bool reflect = false;
        bool hitFromBack = false;
        do {
          rayHit.ray.tfar =
              std::numeric_limits<rayInternal::rtcNumericType>::max();
          rayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
          rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

#if VIENNARAY_EMBREE_VERSION < 4
          // Run the intersection
          rtcIntersect1(rtcScene, &rtcContext, &rayHit);
#else
          rtcIntersect1(rtcScene, &rayHit);
#endif

          /* -------- No hit -------- */
          if (rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            reflect = false;
            break;
          }

          /* -------- Boundary hit -------- */
          if (rayHit.hit.geomID == boundaryID) {
            mBoundary.processHit(rayHit, reflect);
            continue;
          }

          // Calculate point of impact
          const auto &ray = rayHit.ray;
          const rayInternal::rtcNumericType xx =
              ray.org_x + ray.dir_x * ray.tfar;
          const rayInternal::rtcNumericType yy =
              ray.org_y + ray.dir_y * ray.tfar;
          const rayInternal::rtcNumericType zz =
              ray.org_z + ray.dir_z * ray.tfar;
          auto hitPoint = Vec3D<rayInternal::rtcNumericType>{xx, yy, zz};

          /* -------- Hit from back -------- */
          const auto rayDir = Vec3D<T>{ray.dir_x, ray.dir_y, ray.dir_z};
          const auto geomNormal = mGeometry.getPrimNormal(rayHit.hit.primID);
          if (DotProduct(rayDir, geomNormal) > 0) {
            // If the dot product of the ray direction and the surface normal is
            // greater than zero, then we hit the back face of the disk.
            if (hitFromBack) {
              // if hitFromback == true, then the ray hits the back of a disk
              // the second time. In this case we ignore the ray.
              break;
            }
            hitFromBack = true;
            // Let ray through, i.e., continue.
            reflect = true;
#ifdef ARCH_X86
            reinterpret_cast<__m128 &>(rayHit.ray) =
                _mm_set_ps(1e-4f, zz, yy, xx);
#else
            rayHit.ray.org_x = xx;
            rayHit.ray.org_y = yy;
            rayHit.ray.org_z = zz;
            rayHit.ray.tnear = 1e-4f;
#endif
            // keep ray direction as it is
            continue;
          }

          /* -------- Surface hit -------- */
          assert(rayHit.hit.geomID == geometryID && "Geometry hit ID invalid");

          // get fill and reflection
          const auto fillnDirection =
              particle->surfaceHit(rayDir, geomNormal, reflect, RngState);

          if (mGeometry.getMaterialId(rayHit.hit.primID) != excludeMaterial) {
            // trace in cell set
            std::vector<VolumeParticle<T>> particleStack;
            std::normal_distribution<T> normalDist{meanFreePath[0],
                                                   meanFreePath[1]};
            Vec3D<T> hitPoint{xx, yy, zz};
            particleStack.emplace_back(VolumeParticle<T>{
                hitPoint, rayDir, fillnDirection.first, 0., -1, 0});

            while (!particleStack.empty()) {
              auto volumeParticle = std::move(particleStack.back());
              particleStack.pop_back();

              // trace particle
              while (volumeParticle.energy >= 0) {
                volumeParticle.distance = -1;
                while (volumeParticle.distance < 0)
                  volumeParticle.distance = normalDist(RngState);
                auto traveledDist =
                    volumeParticle.direction * volumeParticle.distance;
                volumeParticle.position =
                    traveledDist + volumeParticle.position;

                if (!checkBoundsPeriodic(volumeParticle.position))
                  break;

                auto newIdx = myCellSet->getIndex(volumeParticle.position);
                if (newIdx < 0)
                  break;

                if (newIdx != volumeParticle.cellId) {
                  volumeParticle.cellId = newIdx;
                  auto fill = particle->collision(volumeParticle, RngState,
                                                  particleStack);
                  path.addGridData(newIdx, fill);
                }
              }
            }
          }

          if (!reflect) {
            break;
          }

          // Update ray direction and origin
          rayInternal::fillRayPosition(rayHit.ray, hitPoint);
          rayInternal::fillRayDirection(rayHit.ray, fillnDirection.second);

        } while (reflect);
      } // end ray tracing for loop

#pragma omp critical
      myCellSet->mergePath(path, mNumRays);
    } // end parallel section

    rtcReleaseGeometry(rtcGeometry);
    rtcReleaseGeometry(rtcBoundary);
  }

private:
  bool checkBounds(const Vec3D<T> &hitPoint) const {
    const auto &min = cellSet->getCellGrid()->minimumExtent;
    const auto &max = cellSet->getCellGrid()->maximumExtent;

    return hitPoint[0] >= min[0] && hitPoint[0] <= max[0] &&
           hitPoint[1] >= min[1] && hitPoint[1] <= max[1] &&
           hitPoint[2] >= min[2] && hitPoint[2] <= max[2];
  }

  bool checkBoundsPeriodic(Vec3D<T> &hitPoint) const {
    const auto &min = cellSet->getCellGrid()->minimumExtent;
    const auto &max = cellSet->getCellGrid()->maximumExtent;

    if constexpr (D == 3) {
      if (hitPoint[2] < min[2] || hitPoint[2] > max[2])
        return false;

      if (hitPoint[0] < min[0]) {
        hitPoint[0] = max[0] - mGridDelta / 2.;
      } else if (hitPoint[0] > max[0]) {
        hitPoint[0] = min[0] + mGridDelta / 2.;
      }

      if (hitPoint[1] < min[1]) {
        hitPoint[1] = max[1] - mGridDelta / 2.;
      } else if (hitPoint[1] > max[1]) {
        hitPoint[1] = min[1] + mGridDelta / 2.;
      }
    } else {
      if (hitPoint[1] < min[1] || hitPoint[1] > max[1])
        return false;

      if (hitPoint[0] < min[0]) {
        hitPoint[0] = max[0] - mGridDelta / 2.;
      } else if (hitPoint[0] > max[0]) {
        hitPoint[0] = min[0] + mGridDelta / 2.;
      }
    }

    return true;
  }

private:
  RTCDevice &mDevice;
  viennaray::Geometry<T, D> &mGeometry;
  viennaray::Boundary<T, D> &mBoundary;
  viennaray::Source<T> &mSource;
  std::unique_ptr<AbstractParticle<T>> const mParticle = nullptr;
  const long long mNumRays;
  const bool mUseRandomSeeds;
  const size_t mRunNumber;
  SmartPointer<DenseCellSet<T, D>> cellSet = nullptr;
  const T mGridDelta = 0.;
  const int excludeMaterial = -1;
};

} // namespace viennacs
