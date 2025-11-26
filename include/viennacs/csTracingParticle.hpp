#pragma once

#include <rayUtil.hpp>

#include <vcRNG.hpp>
#include <vcVectorType.hpp>

namespace viennacs {

using namespace viennacore;

template <typename T> struct VolumeParticle {
  Vec3D<T> position;
  Vec3D<T> direction;
  T energy;
  T distance;
  int cellId;
  int scattered;
};

template <typename T> class AbstractParticle {
public:
  virtual ~AbstractParticle() = default;
  virtual std::unique_ptr<AbstractParticle> clone() const = 0;

  virtual void initNew(RNG &rng) = 0;

  virtual std::pair<T, Vec3D<T>> surfaceHit(const Vec3D<T> &rayDir,
                                            const Vec3D<T> &geomNormal,
                                            bool &reflect, RNG &rng) = 0;
  virtual T getSourceDistributionPower() const = 0;
  virtual std::array<T, 2> getMeanFreePath() const = 0;
  virtual T collision(VolumeParticle<T> &particle, RNG &rng,
                      std::vector<VolumeParticle<T>> &particleStack) = 0;
};

template <typename Derived, typename T>
class Particle : public AbstractParticle<T> {
public:
  std::unique_ptr<AbstractParticle<T>> clone() const final {
    return std::make_unique<Derived>(static_cast<Derived const &>(*this));
  }
  void initNew(RNG &rng) override {}
  std::pair<T, Vec3D<T>> surfaceHit(const Vec3D<T> &rayDir,
                                    const Vec3D<T> &geomNormal, bool &reflect,
                                    RNG &rng) override {
    reflect = false;
    return std::pair<T, Vec3D<T>>{1., Vec3D<T>{0., 0., 0.}};
  }
  T getSourceDistributionPower() const override { return 1.; }
  std::array<T, 2> getMeanFreePath() const override { return {1., 1.}; }
  T collision(VolumeParticle<T> &particle, RNG &rng,
              std::vector<VolumeParticle<T>> &particleStack) override {
    return 0.;
  }

protected:
  // We make clear Particle class needs to be inherited
  Particle() = default;
  Particle(const Particle &) = default;
  Particle(Particle &&) = default;
};

} // namespace viennacs
