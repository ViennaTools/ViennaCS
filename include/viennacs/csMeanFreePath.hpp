#pragma once

#include "csDenseCellSet.hpp"

#include <lsDomain.hpp>
#include <lsToDiskMesh.hpp>

// #include <rayGeometry.hpp>
#include <rayReflection.hpp>

#include <vcKDTree.hpp>
#include <vcUtil.hpp>

namespace viennacs {

using namespace viennacore;

template <class NumericType, int D> class MeanFreePath {
public:
  MeanFreePath(const SmartPointer<DenseCellSet<NumericType, D>> passedCellSet)
      : cellSet(passedCellSet), materialIds(cellSet->getScalarData("Material")),
        numCells(cellSet->getNumberOfCells()) {}

  void setBulkLambda(const NumericType passedBulkLambda) {
    bulkLambda = passedBulkLambda;
  }

  void setMaterial(const int passedMaterial) { material = passedMaterial; }

  void setNumRaysPerCell(const NumericType passedNumRaysPerCell) {
    numRaysPerCell = passedNumRaysPerCell;
  }

  void setReflectionLimit(const int passedReflectionLimit) {
    reflectionLimit = passedReflectionLimit;
  }

  void setRngSeed(const unsigned int passedSeed) { seed = passedSeed; }

  void disableSmoothing() { smoothing = false; }

  void enableSmoothing() { smoothing = true; }

  void apply() {
    Logger::getInstance().addInfo("Calculating mean free path ...").print();
    initGeometry();
    runKernel();
  }

private:
  void runKernel() {
    // thread local data storage
    unsigned numThreads = 1;
#ifdef _OPENMP
    numThreads = omp_get_max_threads();
#endif
    std::vector<std::vector<NumericType>> threadLocalData(numThreads);
    std::vector<std::vector<unsigned>> threadLocalHitCount(numThreads);

#pragma omp parallel
    {
      int threadNum = 0;
#ifdef _OPENMP
      threadNum = omp_get_thread_num();
#endif
      auto &data = threadLocalData[threadNum];
      data.resize(numCells, 0.);
      auto &hitCount = threadLocalHitCount[threadNum];
      hitCount.resize(numCells, 0);
      std::uniform_int_distribution<unsigned> pointDist(0, numPoints - 1);

#pragma omp for schedule(dynamic)
      for (long long idx = 0; idx < numRays; ++idx) {

        if (threadNum == 0 && Logger::getLogLevel() >= 4) {
          util::ProgressBar(idx, numRays);
#ifdef VIENNAPS_PYTHON_BUILD
          if (PyErr_CheckSignals() != 0)
            throw pybind11::error_already_set();
#endif
        }

        // particle specific RNG seed
        auto particleSeed = tea<3>(idx, seed);
        RNG RngState(particleSeed);

        auto pointIdx = pointDist(RngState);
        auto direction = viennaray::ReflectionDiffuse<NumericType, D>(
            surfaceNormals[pointIdx], RngState);
        auto cellIdx = getStartingCell(surfacePoints[pointIdx]);
        auto origin = cellSet->getCellCenter(cellIdx);

        unsigned numReflections = 0;
        while (true) {

          /* -------- Cell Marching -------- */
          std::vector<int> hitCells(1, cellIdx);
          NumericType distance = 0;
          int prevIdx = -1;
          bool hitState = false; // -1 bulk hit, 1 material hit
          // invert direction for faster computation in intersectLineBox
          for (int i = 0; i < D; ++i) {
            direction[i] = 1. / direction[i];
          }

          while (true) {

            hitState = false;
            int currentCell = hitCells.back();
            int nextCell = -1;

            const auto &neighbors = cellSet->getNeighbors(currentCell);
            for (const auto &n : neighbors) {
              if (n < 0 || n == prevIdx) {
                continue;
              }

              if (static_cast<int>(materialIds->at(n)) != material) {
                hitState = true; // could be a hit
                continue;
              }

              auto &cellMin = cellSet->getNode(cellSet->getElement(n)[0]);
              auto &cellMax =
                  cellSet->getNode(cellSet->getElement(n)[D == 2 ? 2 : 6]);

              if (intersectLineBox(origin, direction, cellMin, cellMax,
                                   distance)) {
                nextCell = n;
                break;
              }
            }

            if (nextCell < 0 && hitState) {
              // hit a different material
              cellIdx = currentCell;
              break;
            }

            if (nextCell < 0) {
              // no hit
              distance = bulkLambda;
              break;
            }

            if (distance > bulkLambda) {
              // gas phase hit
              cellIdx = currentCell;
              break;
            }

            prevIdx = currentCell;
            hitCells.push_back(nextCell);
          }

          /* -------- Add to cells -------- */
          for (const auto &c : hitCells) {
            data[c] += distance + gridDelta;
            hitCount[c]++;
          }

          /* -------- Reflect -------- */
          if (!hitState)
            break;

          if (++numReflections >= reflectionLimit)
            break;

          // update origin
          origin = cellSet->getCellCenter(cellIdx);

          // update direction
          if (distance > bulkLambda) {
            // gas phase scatter
            randomDirection(direction, RngState);
          } else {
            // material reflection
            auto closestSurfacePoint = kdTree.findNearest(origin);
            assert(closestSurfacePoint->second < gridDelta);
            direction = viennaray::ReflectionDiffuse<NumericType, D>(
                surfaceNormals[closestSurfacePoint->first], RngState);
          }
        }
      }
    }

    // reduce data
    std::vector<NumericType> result(numCells, 0);
    for (const auto &data : threadLocalData) {
#pragma omp parallel for
      for (int i = 0; i < numCells; ++i) {
        result[i] += data[i];
      }
    }

    // reduce hit counts
    std::vector<NumericType> hitCounts(numCells, 0);
    for (const auto &data : threadLocalHitCount) {
#pragma omp parallel for
      for (int i = 0; i < numCells; ++i) {
        hitCounts[i] += data[i];
      }
    }

    // normalize data
#pragma omp parallel for
    for (int i = 0; i < numCells; ++i) {
      if (hitCounts[i] > 0)
        result[i] = result[i] / hitCounts[i];
      else
        result[i] = 0.;
    }

    // smooth result
    auto finalResult = cellSet->addScalarData("MeanFreePath");
    materialIds = cellSet->getScalarData("Material");
#pragma omp parallel for
    for (int i = 0; i < numCells; i++) {
      if (static_cast<int>(materialIds->at(i)) != material)
        continue;

      if (smoothing) {
        const auto &neighbors = cellSet->getNeighbors(i);
        NumericType sum = 0;
        unsigned count = 0;
        for (const auto &n : neighbors) {
          if (n < 0 || static_cast<int>(materialIds->at(i)) != material)
            continue;
          sum += result[n];
          count++;
        }
        if (count > 0)
          finalResult->at(i) = sum / count;
      } else {
        finalResult->at(i) = result[i];
      }
    }

    if (Logger::getLogLevel() >= 4) {
      std::cout << std::endl;
    }
  }

  void initGeometry() {
    auto mesh = SmartPointer<viennals::Mesh<NumericType>>::New();
    viennals::ToDiskMesh<NumericType, D>(cellSet->getSurface(), mesh).apply();
    surfacePoints = mesh->getNodes();
    surfaceNormals = *mesh->getCellData().getVectorData("Normals");
    numPoints = surfacePoints.size();

    kdTree.setPoints(surfacePoints);
    kdTree.build();

    gridDelta = cellSet->getGridDelta();
    numRays = static_cast<long long>(numCells * numRaysPerCell);
  }

  int getStartingCell(const Vec3D<NumericType> &origin) const {
    int cellIdx = cellSet->getIndex(origin);
    if (cellIdx < 0) {
      Logger::getInstance()
          .addError("No starting cell found for ray " +
                    std::to_string(origin[0]) + " " +
                    std::to_string(origin[1]) + " " + std::to_string(origin[2]))
          .print();
    }
    if (static_cast<int>(materialIds->at(cellIdx)) != material) {
      const auto &neighbors = cellSet->getNeighbors(cellIdx);
      for (const auto &n : neighbors) {
        if (n >= 0 && static_cast<int>(materialIds->at(n)) == material) {
          cellIdx = n;
          break;
        }
      }
    }
    return cellIdx;
  }

  // https://gamedev.stackexchange.com/a/18459
  static bool intersectLineBox(const Vec3D<NumericType> &origin,
                               const Vec3D<NumericType> &direction,
                               const Vec3D<NumericType> &min,
                               const Vec3D<NumericType> &max,
                               NumericType &distance) {
    Vec3D<NumericType> t1, t2;
    for (int i = 0; i < D; ++i) {
      // direction is inverted
      t1[i] = (min[i] - origin[i]) * direction[i];
      t2[i] = (max[i] - origin[i]) * direction[i];
    }
    NumericType tmin, tmax;
    if constexpr (D == 2) {
      tmin = std::max(std::min(t1[0], t2[0]), std::min(t1[1], t2[1]));
      tmax = std::min(std::max(t1[0], t2[0]), std::max(t1[1], t2[1]));
    } else {
      tmin = std::max(std::max(std::min(t1[0], t2[0]), std::min(t1[1], t2[1])),
                      std::min(t1[2], t2[2]));
      tmax = std::min(std::min(std::max(t1[0], t2[0]), std::max(t1[1], t2[1])),
                      std::max(t1[2], t2[2]));
    }

    if (tmax > 0 && tmin < tmax) {
      // ray intersects box
      distance = tmin;
      return true;
    }

    return false;
  }

  static void randomDirection(Vec3D<NumericType> &direction, RNG &rngState) {
    std::uniform_real_distribution<NumericType> dist(-1, 1);
    for (int i = 0; i < D; ++i) {
      direction[i] = dist(rngState);
    }
    rayInternal::Normalize(direction);
  }

private:
  SmartPointer<DenseCellSet<NumericType, D>> const cellSet = nullptr;
  KDTree<NumericType, Vec3D<NumericType>> kdTree;

  std::vector<NumericType> const *materialIds;
  std::vector<Vec3D<NumericType>> surfaceNormals;
  std::vector<Vec3D<NumericType>> surfacePoints;

  NumericType bulkLambda = 0;
  NumericType gridDelta = 0;
  unsigned int const numCells = 0;
  unsigned int numPoints = 0;
  unsigned int seed = 15235135;
  unsigned int reflectionLimit = 100;
  long long numRays = 0;
  NumericType numRaysPerCell = 1000;
  bool smoothing = true;
  int material = 0;
};

/* Old version

template <class NumericType, int D> class MeanFreePath {

private:
  using levelSetsType =
      SmartPointer<std::vector<SmartPointer<viennals::Domain<NumericType, D>>>>;
  using cellSetType = SmartPointer<DenseCellSet<NumericType, D>>;

public:
  MeanFreePath() : traceDevice(rtcNewDevice("hugepages=1")) {
    static_assert(D == 2 &&
                  "Mean free path calculation only implemented for 2D");
  }

  ~MeanFreePath() { rtcReleaseDevice(traceDevice); }

  void setLevelSets(levelSetsType passedLevelSets) {
    levelSets = passedLevelSets;
  }

  void setCellSet(cellSetType passedCellSet) { cellSet = passedCellSet; }

  void setBulkLambda(const NumericType passedBulkLambda) {
    bulkLambda = passedBulkLambda;
  }

  template <class Material> void setMaterial(const Material passedMaterial) {
    material = static_cast<int>(passedMaterial);
  }

  void setNumRaysPerCell(const int passedNumRaysPerCell) {
    numRaysPerCell = passedNumRaysPerCell;
  }

  NumericType getMaxLambda() const { return maxLambda; }

  void apply() {
    Logger::getInstance().addInfo("Calculating mean free path ...").print();
    cellSet->addScalarData("MeanFreePath", 0.);
    runKernel();
  }

private:
  void runKernel() {
#ifdef ARCH_X86
    // for best performance set FTZ and DAZ flags in MXCSR control and status
    // register
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
    unsigned numCells = cellSet->getElements().size();
    auto data = cellSet->getScalarData("MeanFreePath");
    auto materials = cellSet->getScalarData("Material");

    auto traceGeometry = rayGeometry<NumericType, D>();
    {
      auto mesh = SmartPointer<lsMesh<NumericType>>::New();
      lsToDiskMesh<NumericType, D>(levelSets->back(), mesh).apply();
      auto &points = mesh->getNodes();
      auto normals = mesh->getCellData().getVectorData("Normals");
      gridDelta = levelSets->back()->getGrid().getGridDelta();
      traceGeometry.initGeometry(traceDevice, points, *normals,
                                 gridDelta * rayInternal::DiskFactor<D>);
    }

    auto rtcScene = rtcNewScene(traceDevice);
    rtcSetSceneFlags(rtcScene, RTC_SCENE_FLAG_NONE);
    rtcSetSceneBuildQuality(rtcScene, RTC_BUILD_QUALITY_HIGH);
    auto rtcGeometry = traceGeometry.getRTCGeometry();
    auto geometryID = rtcAttachGeometry(rtcScene, rtcGeometry);
    assert(rtcGetDeviceError(traceDevice) == RTC_ERROR_NONE &&
           "Embree device error");

    maxLambda = 0.;

#pragma omp parallel reduction(max : maxLambda)
    {
      rtcJoinCommitScene(rtcScene);

      alignas(128) auto rayHit =
          RTCRayHit{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

#if VIENNARAY_EMBREE_VERSION < 4
      auto rtcContext = RTCIntersectContext{};
      rtcInitIntersectContext(&rtcContext);
#endif

#pragma omp for
      for (int idx = 0; idx < numCells; ++idx) {
        if (static_cast<int>(materials->at(idx)) != material)
          continue;

        auto cellCenter = cellSet->getCellCenter(idx);
        auto &ray = rayHit.ray;
#ifdef ARCH_X86
        reinterpret_cast<__m128 &>(ray) =
            _mm_set_ps(1e-4f, (float)cellCenter[2], (float)cellCenter[1],
                       (float)cellCenter[0]);
#else
        ray.org_x = (float)cellCenter[0];
        ray.org_y = (float)cellCenter[1];
        ray.org_z = (float)cellCenter[2];
        ray.tnear = 1e-4f;
#endif

#ifdef VIENNARAY_USE_RAY_MASKING
        ray.mask = -1;
#endif
        for (unsigned cIdx = 0; cIdx < numRaysPerCell; ++cIdx) {

          auto direction = getDirection(cIdx);

#ifdef ARCH_X86
          reinterpret_cast<__m128 &>(ray.dir_x) =
              _mm_set_ps(0.0f, (float)direction[2], (float)direction[1],
                         (float)direction[0]);
#else
          ray.dir_x = (float)direction[0];
          ray.dir_y = (float)direction[1];
          ray.dir_z = (float)direction[2];
          ray.time = 0.0f;
#endif

          ray.tfar = std::numeric_limits<rayInternal::rtcNumericType>::max();
          rayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
          rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

#if VIENNARAY_EMBREE_VERSION < 4
          // Run the intersection
          rtcIntersect1(rtcScene, &rtcContext, &rayHit);
#else
          rtcIntersect1(rtcScene, &rayHit);
#endif

          if (rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            data->at(idx) += bulkLambda;
            continue;
          }

          const auto rayDir =
              Vec3D<NumericType>{ray.dir_x, ray.dir_y, ray.dir_z};
          auto geomNormal = traceGeometry.getPrimNormal(rayHit.hit.primID);

          if (DotProduct(rayDir, geomNormal) > 0) {
            continue;
          }

          data->at(idx) += ray.tfar;
        }

        data->at(idx) /= numRaysPerCell;
        maxLambda = std::max(maxLambda, data->at(idx));
      }
    } // end of parallel section

    traceGeometry.releaseGeometry();
    rtcReleaseScene(rtcScene);
  }

  Vec3D<NumericType> getDirection(const unsigned int idx) {
    Vec3D<NumericType> direction;
    NumericType theta = idx * 2. * M_PI / numRaysPerCell;
    direction[0] = std::cos(theta);
    direction[1] = std::sin(theta);
    direction[2] = 0.;
    return direction;
  }

private:
  levelSetsType levelSets = nullptr;
  cellSetType cellSet = nullptr;
  RTCDevice traceDevice;

  NumericType gridDelta = 0;
  NumericType bulkLambda = 0;
  NumericType maxLambda = 0.;
  long numRaysPerCell = 100;
  int material = -1;
};
*/

} // namespace viennacs
