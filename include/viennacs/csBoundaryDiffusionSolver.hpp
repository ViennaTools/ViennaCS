#pragma once

#include "csDenseCellSet.hpp"
#include "csDiffusionSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace viennacs {

enum class BoundaryConditionType { Neumann, Dirichlet, Robin };

template <class NumericType> struct BoundaryCondition {
  BoundaryConditionType type = BoundaryConditionType::Neumann;
  NumericType value = 0.;
  NumericType transferCoefficient = 0.;

  static BoundaryCondition dirichlet(NumericType boundaryValue) {
    BoundaryCondition condition;
    condition.type = BoundaryConditionType::Dirichlet;
    condition.value = boundaryValue;
    return condition;
  }

  static BoundaryCondition neumann(NumericType outwardFlux = 0.) {
    BoundaryCondition condition;
    condition.type = BoundaryConditionType::Neumann;
    condition.value = outwardFlux;
    return condition;
  }

  static BoundaryCondition robin(NumericType transfer,
                                 NumericType exteriorValue) {
    BoundaryCondition condition;
    condition.type = BoundaryConditionType::Robin;
    condition.transferCoefficient = transfer;
    condition.value = exteriorValue;
    return condition;
  }
};

/// Cell-centred finite-difference diffusion with embedded level-set boundary
/// points on Cartesian cell edges.
///
/// For regular same-material neighbours this uses the standard cell-centred
/// distance dx. If an axis direction leaves the diffusive material through an
/// embedded boundary point, the stencil uses the sub-grid centre-to-boundary
/// distance and applies the configured boundary condition at that point. This
/// mirrors the nonuniform boundary-side treatment used by the ViennaLS
/// oxidation diffusion solve, adapted from grid nodes to ViennaCS cells.
template <class NumericType, int D> class BoundaryDiffusionSolver {
  using CellSetType = DenseCellSet<NumericType, D>;
  using CellSetPointer = SmartPointer<CellSetType>;
  using BoundaryPoint = typename CellSetType::EmbeddedBoundaryPoint;

  struct StencilSide {
    NumericType distance = 1.;
    NumericType nodeCoefficient = 1.;
    NumericType constant = 0.;
  };

public:
  void setCellSet(CellSetPointer cellSet) { cellSet_ = cellSet; }

  void setMode(DiffusionSolverMode mode) { mode_ = mode; }

  void setMaxIterations(int n) { maxIterations_ = std::max(n, 1); }

  void setTolerance(NumericType tol) {
    tolerance_ = std::max(tol, NumericType(1e-12));
  }

  void setRelaxation(NumericType omega) {
    relaxation_ = std::clamp(omega, NumericType(0.1), NumericType(1.95));
  }

  void setClampNonNegative(bool value) { clampNonNegative_ = value; }

  /// Per-cell volumetric source term S(x) added to the RHS: ∂c/∂t = D∇²c + S.
  /// The pointer must remain valid for the lifetime of the step call.
  void setSourceField(const std::vector<NumericType> *source) {
    source_ = source;
  }
  void clearSourceField() { source_ = nullptr; }

  /// Returns the effective minimum stencil distance the solver would use at
  /// boundary cells: max(rawMin, minBoundaryDistanceFactor * dx).
  /// Use this as the effective dx in a CFL stability calculation for explicit
  /// time stepping — it can be much smaller than dx near sharp boundaries.
  NumericType getEffectiveMinBoundaryDistance(NumericType dx) const {
    if (!cellSet_)
      return dx;
    return std::max(cellSet_->getMinFaceBoundaryDistance(),
                    minBoundaryDistanceFactor_ * dx);
  }

  void setDiffusiveMaterials(const std::vector<int> &materials) {
    diffusiveMaterials_ = {materials.begin(), materials.end()};
  }

  void setSourceMaterials(const std::vector<int> &materials) {
    sourceMaterials_ = {materials.begin(), materials.end()};
  }

  void setBlockedMaterials(const std::vector<int> &materials) {
    blockedMaterials_ = {materials.begin(), materials.end()};
  }

  void setBoundaryCondition(unsigned levelSetIndex,
                            const BoundaryCondition<NumericType> &condition) {
    boundaryConditions_[levelSetIndex] = condition;
  }

  void clearBoundaryConditions() { boundaryConditions_.clear(); }

  void
  setDefaultBoundaryCondition(const BoundaryCondition<NumericType> &condition) {
    defaultBoundaryCondition_ = condition;
    useDefaultBoundaryCondition_ = true;
  }

  void step(std::vector<NumericType> &field,
            const std::vector<NumericType> &materials, NumericType dx,
            NumericType dt, NumericType diffusivity) {
    auto diffusion = [diffusivity](int) { return diffusivity; };
    stepVariable(field, materials, dx, dt, diffusion);
  }

  template <typename DiffusivityFn>
  void stepVariable(std::vector<NumericType> &field,
                    const std::vector<NumericType> &materials, NumericType dx,
                    NumericType dt, const DiffusivityFn &diffusivityAt) {
    if (!cellSet_)
      return;

    switch (mode_) {
    case DiffusionSolverMode::Explicit:
      stepExplicit(field, materials, dx, dt, diffusivityAt);
      break;
    case DiffusionSolverMode::GaussSeidel:
      stepGaussSeidel(field, materials, dx, dt, diffusivityAt);
      break;
#ifdef VIENNACS_USE_EIGEN
    case DiffusionSolverMode::EigenSparseLU:
      Logger::getInstance()
          .addWarning("BoundaryDiffusionSolver: EigenSparseLU is not "
                      "available for embedded-boundary stencils; using "
                      "GaussSeidel.")
          .print();
      stepGaussSeidel(field, materials, dx, dt, diffusivityAt);
      break;
#endif
    }
  }

private:
  bool isDiffusive(int material) const {
    if (diffusiveMaterials_.empty())
      return true;
    return diffusiveMaterials_.count(material) > 0;
  }

  bool isBlocked(int material) const {
    return blockedMaterials_.count(material) > 0;
  }

  bool isSource(int material) const {
    return sourceMaterials_.count(material) > 0;
  }

  bool contributes(int material) const {
    return !isBlocked(material) &&
           (isDiffusive(material) || isSource(material));
  }

  template <typename DiffusivityFn>
  void stepExplicit(std::vector<NumericType> &field,
                    const std::vector<NumericType> &materials, NumericType dx,
                    NumericType dt, const DiffusivityFn &diffusivityAt) {
    buffer_ = field;

#pragma omp parallel for
    for (int cellId = 0; cellId < static_cast<int>(field.size()); ++cellId) {
      const int material = static_cast<int>(materials[cellId]);
      if (!isDiffusive(material) || isBlocked(material))
        continue;

      const NumericType diffusivity =
          std::max(diffusivityAt(cellId), NumericType(0));
      NumericType rhs = 0.;
      for (unsigned direction = 0; direction < D; ++direction) {
        const auto negative = makeStencilSide(cellId, direction, -1, field,
                                              materials, dx, diffusivity);
        const auto positive = makeStencilSide(cellId, direction, 1, field,
                                              materials, dx, diffusivity);
        addAxisContribution(rhs, negative, positive, diffusivity,
                            field[cellId]);
      }

      const NumericType src = source_ ? (*source_)[cellId] : NumericType(0);
      const NumericType updated = field[cellId] + dt * (rhs + src);
      buffer_[cellId] = (clampNonNegative_ && updated < NumericType(0))
                            ? NumericType(0)
                            : updated;
    }

    field.swap(buffer_);
  }

  template <typename DiffusivityFn>
  void stepGaussSeidel(std::vector<NumericType> &field,
                       const std::vector<NumericType> &materials,
                       NumericType dx, NumericType dt,
                       const DiffusivityFn &diffusivityAt) {
    buffer_ = field;

    for (int iteration = 0; iteration < maxIterations_; ++iteration) {
      NumericType maxDelta = 0.;
      NumericType maxValue = 0.;

      for (int cellId = 0; cellId < static_cast<int>(field.size()); ++cellId) {
        const int material = static_cast<int>(materials[cellId]);
        if (!isDiffusive(material) || isBlocked(material))
          continue;

        const NumericType diffusivity =
            std::max(diffusivityAt(cellId), NumericType(0));
        NumericType rightHandSide = field[cellId];
        NumericType diagonal = 1.;

        for (unsigned direction = 0; direction < D; ++direction) {
          const auto negative = makeStencilSide(cellId, direction, -1, buffer_,
                                                materials, dx, diffusivity);
          const auto positive = makeStencilSide(cellId, direction, 1, buffer_,
                                                materials, dx, diffusivity);
          addImplicitAxisContribution(rightHandSide, diagonal, negative,
                                      positive, diffusivity, dt);
        }

        if (source_)
          rightHandSide += dt * (*source_)[cellId];
        const NumericType gsValue = rightHandSide / diagonal;
        const NumericType previous = buffer_[cellId];
        const NumericType relaxed =
            previous + relaxation_ * (gsValue - previous);
        buffer_[cellId] = (clampNonNegative_ && relaxed < NumericType(0))
                              ? NumericType(0)
                              : relaxed;
        maxDelta = std::max(maxDelta, std::abs(buffer_[cellId] - previous));
        maxValue = std::max(maxValue, std::abs(buffer_[cellId]));
      }

      if (maxDelta / std::max(maxValue, NumericType(1e-30)) < tolerance_)
        break;
    }

    field.swap(buffer_);
  }

  StencilSide makeStencilSide(int cellId, unsigned direction, int offset,
                              const std::vector<NumericType> &field,
                              const std::vector<NumericType> &materials,
                              NumericType dx, NumericType diffusivity) const {
    const unsigned faceIdx = direction * 2u + (offset > 0 ? 1u : 0u);
    const int pointId = cellSet_->getFaceBoundaryPointId(cellId, faceIdx);
    if (pointId >= 0) {
      const auto &point = cellSet_->getEmbeddedBoundaryPoints()[pointId];
      if (hasBoundaryCondition(point.levelSetIndex)) {
        const NumericType dist =
            std::max(cellSet_->getFaceBoundaryDistance(cellId, faceIdx),
                     minBoundaryDistanceFactor_ * dx);
        return boundarySide(point.levelSetIndex, dist, diffusivity, dx);
      }
    }

    const auto &neighbors = cellSet_->getNeighbors(cellId);
    const int neighbor = neighbors[2 * direction + (offset > 0 ? 1 : 0)];
    if (neighbor >= 0 && contributes(static_cast<int>(materials[neighbor])))
      return {dx, 0., field[neighbor]};

    return zeroFluxSide(dx);
  }

  StencilSide boundarySide(unsigned levelSetIndex, NumericType distance,
                           NumericType diffusivity, NumericType dx) const {
    const auto condition = boundaryCondition(levelSetIndex);
    switch (condition.type) {
    case BoundaryConditionType::Dirichlet:
      return {distance, 0., condition.value};
    case BoundaryConditionType::Neumann:
      if (std::abs(condition.value) <=
          std::numeric_limits<NumericType>::epsilon())
        return zeroFluxSide(dx);
      return {distance, 1.,
              -condition.value * distance /
                  std::max(diffusivity, NumericType(1e-30))};
    case BoundaryConditionType::Robin: {
      const NumericType conductance = diffusivity / distance;
      const NumericType denominator =
          conductance + condition.transferCoefficient;
      if (denominator <= std::numeric_limits<NumericType>::epsilon())
        return zeroFluxSide(dx);
      return {distance, conductance / denominator,
              condition.transferCoefficient * condition.value / denominator};
    }
    }
    return zeroFluxSide(dx);
  }

  BoundaryCondition<NumericType>
  boundaryCondition(unsigned levelSetIndex) const {
    const auto found = boundaryConditions_.find(levelSetIndex);
    if (found != boundaryConditions_.end())
      return found->second;
    return defaultBoundaryCondition_;
  }

  bool hasBoundaryCondition(unsigned levelSetIndex) const {
    return boundaryConditions_.find(levelSetIndex) !=
               boundaryConditions_.end() ||
           useDefaultBoundaryCondition_;
  }

  StencilSide zeroFluxSide(NumericType dx) const { return {dx, 1., 0.}; }

  void addAxisContribution(NumericType &rhs, const StencilSide &negative,
                           const StencilSide &positive, NumericType diffusivity,
                           NumericType nodeValue) const {
    const NumericType distanceSum = negative.distance + positive.distance;
    if (distanceSum <= std::numeric_limits<NumericType>::epsilon())
      return;
    addSideContribution(rhs, negative, distanceSum, diffusivity, nodeValue);
    addSideContribution(rhs, positive, distanceSum, diffusivity, nodeValue);
  }

  void addSideContribution(NumericType &rhs, const StencilSide &side,
                           NumericType distanceSum, NumericType diffusivity,
                           NumericType nodeValue) const {
    if (side.distance <= std::numeric_limits<NumericType>::epsilon())
      return;

    const NumericType coefficient =
        NumericType(2) * diffusivity / (side.distance * distanceSum);
    rhs += coefficient * (side.constant +
                          (side.nodeCoefficient - NumericType(1)) * nodeValue);
  }

  void addImplicitAxisContribution(NumericType &rightHandSide,
                                   NumericType &diagonal,
                                   const StencilSide &negative,
                                   const StencilSide &positive,
                                   NumericType diffusivity,
                                   NumericType dt) const {
    const NumericType distanceSum = negative.distance + positive.distance;
    if (distanceSum <= std::numeric_limits<NumericType>::epsilon())
      return;
    addImplicitSideContribution(rightHandSide, diagonal, negative, distanceSum,
                                diffusivity, dt);
    addImplicitSideContribution(rightHandSide, diagonal, positive, distanceSum,
                                diffusivity, dt);
  }

  void
  addImplicitSideContribution(NumericType &rightHandSide, NumericType &diagonal,
                              const StencilSide &side, NumericType distanceSum,
                              NumericType diffusivity, NumericType dt) const {
    if (side.distance <= std::numeric_limits<NumericType>::epsilon())
      return;

    const NumericType coefficient =
        NumericType(2) * diffusivity / (side.distance * distanceSum);
    rightHandSide += dt * coefficient * side.constant;
    diagonal += dt * coefficient * (NumericType(1) - side.nodeCoefficient);
  }

  CellSetPointer cellSet_ = nullptr;
  const std::vector<NumericType> *source_ = nullptr;
  DiffusionSolverMode mode_ = DiffusionSolverMode::GaussSeidel;
  bool clampNonNegative_ = true;
  int maxIterations_ = 400;
  NumericType tolerance_ = NumericType(1e-6);
  NumericType relaxation_ = NumericType(1);
  NumericType minBoundaryDistanceFactor_ = NumericType(1e-6);
  std::unordered_set<int> diffusiveMaterials_;
  std::unordered_set<int> sourceMaterials_;
  std::unordered_set<int> blockedMaterials_;
  std::unordered_map<unsigned, BoundaryCondition<NumericType>>
      boundaryConditions_;
  BoundaryCondition<NumericType> defaultBoundaryCondition_{};
  bool useDefaultBoundaryCondition_ = false;
  std::vector<NumericType> buffer_;
};

} // namespace viennacs
