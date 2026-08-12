#pragma once

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "csDenseCellSet.hpp"
#include <vcLogger.hpp>

#ifdef VIENNACS_USE_EIGEN
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#endif

namespace viennacs {

enum class DiffusionSolverMode {
  Explicit,    // forward Euler, stability-limited dt
  GaussSeidel, // backward Euler, matrix-free iterative solve
#ifdef VIENNACS_USE_EIGEN
  EigenSparseLU, // backward Euler, direct sparse LU factorisation
#endif
};

// Finite-difference diffusion solver on a DenseCellSet neighbourhood graph.
//
// Three material roles (set by the corresponding setters):
//   Diffusive  — solved each step (default: all materials when list is empty)
//   Source     — Dirichlet / fixed-value cells: not updated but contribute
//                flux to diffusive neighbours
//   Blocked    — hard walls: not updated, invisible to neighbours
template <class NumericType, int D> class DiffusionSolver {
public:
  void setCellSet(SmartPointer<DenseCellSet<NumericType, D>> cellSet) {
    cellSet_ = cellSet;
    invalidateActiveCache();
    invalidateEigenCache();
  }

  void setMode(DiffusionSolverMode mode) { mode_ = mode; }

  // GaussSeidel convergence options
  void setMaxIterations(int n) { maxIterations_ = std::max(n, 1); }
  void setTolerance(NumericType tol) {
    tolerance_ = std::max(tol, NumericType(1e-12));
  }
  void setRelaxation(NumericType omega) {
    relaxation_ = std::clamp(omega, NumericType(0.1), NumericType(1.95));
  }

  void setClampNonNegative(bool v) { clampNonNegative_ = v; }

  /// Per-cell volumetric source term S(x) added to the RHS: ∂c/∂t = D∇²c + S.
  /// The pointer must remain valid for the lifetime of the step call.
  /// Pass nullptr (or call clearSourceField) to disable.
  void setSourceField(const std::vector<NumericType> *source) {
    source_ = source;
  }
  void clearSourceField() { source_ = nullptr; }

  void setDiffusiveMaterials(const std::vector<int> &m) {
    diffusiveMaterials_ = {m.begin(), m.end()};
    invalidateActiveCache();
    invalidateEigenCache();
  }

  void setSourceMaterials(const std::vector<int> &m) {
    sourceMaterials_ = {m.begin(), m.end()};
    invalidateActiveCache();
    invalidateEigenCache();
  }

  void setBlockedMaterials(const std::vector<int> &m) {
    blockedMaterials_ = {m.begin(), m.end()};
    invalidateActiveCache();
    invalidateEigenCache();
  }

  // Advance `field` in-place by `dt` with a spatially constant diffusivity.
  // The Eigen factorisation is cached and reused when dt and diffusivity are
  // unchanged between calls.
  void step(std::vector<NumericType> &field,
            const std::vector<NumericType> &materials, NumericType dx,
            NumericType dt, NumericType diffusivity) {
    auto Dfn = [d = diffusivity](int) { return d; };
    dispatch(field, materials, dx, dt, Dfn, false, diffusivity, dt);
  }

  // Advance `field` in-place by `dt` with a per-cell diffusivity given by
  // `diffusivityAt(cellIndex)`. The Eigen factorisation is always rebuilt.
  template <typename DiffusivityFn>
  void stepVariable(std::vector<NumericType> &field,
                    const std::vector<NumericType> &materials, NumericType dx,
                    NumericType dt, const DiffusivityFn &diffusivityAt) {
    dispatch(field, materials, dx, dt, diffusivityAt, true, NumericType(-1),
             dt);
  }

private:
  // -----------------------------------------------------------------------
  // Material predicates
  // -----------------------------------------------------------------------
  bool isDiffusive(int mat) const {
    if (diffusiveMaterials_.empty())
      return true;
    return diffusiveMaterials_.count(mat) > 0;
  }
  bool isBlocked(int mat) const { return blockedMaterials_.count(mat) > 0; }
  bool isSource(int mat) const { return sourceMaterials_.count(mat) > 0; }

  bool contributes(int mat) const {
    return !isBlocked(mat) && (isDiffusive(mat) || isSource(mat));
  }

  void invalidateActiveCache() {
    activeCacheDirty_ = true;
    activeMaterialsPtr_ = nullptr;
    activeMaterialsSize_ = 0;
  }

  void ensureActiveCache(const std::vector<NumericType> &materials) {
    if (!activeCacheDirty_ && activeMaterialsPtr_ == materials.data() &&
        activeMaterialsSize_ == materials.size())
      return;

    activeCells_.clear();
    activeNeighbors_.clear();
    if (!cellSet_)
      return;

    activeCells_.reserve(materials.size());
    activeNeighbors_.reserve(materials.size());
    for (int i = 0; i < static_cast<int>(materials.size()); ++i) {
      const auto mat = static_cast<int>(materials[i]);
      if (!isDiffusive(mat) || isBlocked(mat))
        continue;

      std::vector<int> neighbors;
      neighbors.reserve(2 * D);
      for (const auto n : cellSet_->getNeighbors(i)) {
        if (n < 0)
          continue;
        if (!contributes(static_cast<int>(materials[n])))
          continue;
        neighbors.push_back(n);
      }
      activeCells_.push_back(i);
      activeNeighbors_.push_back(std::move(neighbors));
    }

    activeMaterialsPtr_ = materials.data();
    activeMaterialsSize_ = materials.size();
    activeCacheDirty_ = false;
  }

  // -----------------------------------------------------------------------
  // Dispatch to the selected solver
  // -----------------------------------------------------------------------
  template <typename DiffusivityFn>
  void dispatch(std::vector<NumericType> &field,
                const std::vector<NumericType> &materials, NumericType dx,
                NumericType dt, const DiffusivityFn &Dfn, bool variableD,
                NumericType constantD, NumericType requestedDt) {
    switch (mode_) {
    case DiffusionSolverMode::Explicit:
      stepExplicit(field, materials, dx, dt, Dfn);
      break;
    case DiffusionSolverMode::GaussSeidel:
      stepGaussSeidel(field, materials, dx, dt, Dfn);
      break;
#ifdef VIENNACS_USE_EIGEN
    case DiffusionSolverMode::EigenSparseLU:
      stepEigen(field, materials, dx, dt, Dfn, variableD, constantD,
                requestedDt);
      break;
#endif
    }
  }

  // -----------------------------------------------------------------------
  // Explicit forward Euler
  // -----------------------------------------------------------------------
  template <typename DiffusivityFn>
  void stepExplicit(std::vector<NumericType> &field,
                    const std::vector<NumericType> &materials, NumericType dx,
                    NumericType dt, const DiffusivityFn &Dfn) {
    const auto invDx2 = NumericType(1) / (dx * dx);
    ensureActiveCache(materials);
    buffer_ = field;

#pragma omp parallel for
    for (int k = 0; k < static_cast<int>(activeCells_.size()); ++k) {
      const auto i = activeCells_[k];
      const auto c = field[i];

      NumericType lap = NumericType(0);
      const auto &neighbors = activeNeighbors_[k];
      for (const auto n : neighbors)
        lap += field[n] - c;
      if (neighbors.empty())
        continue;

      const NumericType src = source_ ? (*source_)[i] : NumericType(0);
      auto v = c + dt * (Dfn(i) * invDx2 * lap + src);
      buffer_[i] =
          (clampNonNegative_ && v < NumericType(0)) ? NumericType(0) : v;
    }

    field.swap(buffer_);
  }

  // -----------------------------------------------------------------------
  // Implicit Gauss-Seidel (backward Euler, matrix-free)
  // -----------------------------------------------------------------------
  template <typename DiffusivityFn>
  void stepGaussSeidel(std::vector<NumericType> &field,
                       const std::vector<NumericType> &materials,
                       NumericType dx, NumericType dt,
                       const DiffusivityFn &Dfn) {
    const auto invDx2 = NumericType(1) / (dx * dx);
    ensureActiveCache(materials);
    buffer_ = field;

    for (int iter = 0; iter < maxIterations_; ++iter) {
      NumericType maxDelta = NumericType(0);
      NumericType maxVal = NumericType(0);

      for (int k = 0; k < static_cast<int>(activeCells_.size()); ++k) {
        const auto i = activeCells_[k];
        const auto &neighbors = activeNeighbors_[k];
        if (neighbors.empty()) {
          maxVal = std::max(maxVal, std::abs(buffer_[i]));
          continue;
        }

        NumericType nsum = NumericType(0);
        for (const auto n : neighbors)
          nsum += buffer_[n];

        const auto alpha = dt * Dfn(i) * invDx2;
        const NumericType src = source_ ? (*source_)[i] : NumericType(0);
        const auto gs = (field[i] + alpha * nsum + dt * src) /
                        (NumericType(1) + alpha * neighbors.size());
        const auto prev = buffer_[i];
        const auto upd = prev + relaxation_ * (gs - prev);
        buffer_[i] =
            (clampNonNegative_ && upd < NumericType(0)) ? NumericType(0) : upd;
        maxDelta = std::max(maxDelta, std::abs(buffer_[i] - prev));
        maxVal = std::max(maxVal, std::abs(buffer_[i]));
      }

      if (maxDelta / std::max(maxVal, NumericType(1e-30)) < tolerance_)
        break;
    }

    field.swap(buffer_);
  }

#ifdef VIENNACS_USE_EIGEN
  // -----------------------------------------------------------------------
  // Implicit Eigen SparseLU (backward Euler, direct factorisation)
  // -----------------------------------------------------------------------
  void invalidateEigenCache() {
    eigenMapDirty_ = true;
    eigenFactorized_ = false;
  }

  void buildEigenMap(const std::vector<NumericType> &materials) {
    eigenCellMap_.clear();
    eigenNumCells_ = 0;
    for (int i = 0; i < static_cast<int>(materials.size()); ++i) {
      const auto mat = static_cast<int>(materials[i]);
      if (isDiffusive(mat) || isSource(mat))
        eigenCellMap_[i] = eigenNumCells_++;
    }
    eigenMapDirty_ = false;
    eigenFactorized_ = false;
  }

  template <typename DiffusivityFn>
  void assembleAndFactorise(const std::vector<NumericType> &materials,
                            NumericType dx, NumericType dt,
                            const DiffusivityFn &Dfn) {
    const auto invDx2 = NumericType(1) / (dx * dx);
    std::vector<Eigen::Triplet<NumericType>> triplets;
    triplets.reserve(static_cast<std::size_t>(2 * D * eigenNumCells_));

    for (const auto &[ci, si] : eigenCellMap_) {
      const auto mat = static_cast<int>(materials[ci]);
      if (!isDiffusive(mat) || isBlocked(mat)) {
        triplets.push_back({si, si, NumericType(1)}); // Dirichlet identity row
        continue;
      }

      NumericType alphaSum = NumericType(0);
      for (const auto n : cellSet_->getNeighbors(ci)) {
        if (n < 0 || !contributes(static_cast<int>(materials[n])))
          continue;
        const auto it = eigenCellMap_.find(n);
        if (it == eigenCellMap_.end())
          continue;
        const auto alpha = dt * std::max(Dfn(ci), NumericType(0)) * invDx2;
        triplets.emplace_back(si, it->second, -alpha);
        alphaSum += alpha;
      }
      triplets.emplace_back(si, si, NumericType(1) + alphaSum);
    }

    Eigen::SparseMatrix<NumericType> A(eigenNumCells_, eigenNumCells_);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    eigenSolver_.compute(A);
    if (eigenSolver_.info() != Eigen::Success)
      Logger::getInstance()
          .addWarning("DiffusionSolver: Eigen factorisation failed.")
          .print();
    eigenFactorized_ = true;
  }

  template <typename DiffusivityFn>
  void stepEigen(std::vector<NumericType> &field,
                 const std::vector<NumericType> &materials, NumericType dx,
                 NumericType dt, const DiffusivityFn &Dfn, bool variableD,
                 NumericType constantD, NumericType requestedDt) {
    if (!cellSet_)
      return;

    if (eigenMapDirty_)
      buildEigenMap(materials);

    const bool dtChanged = (requestedDt != eigenCachedDt_);
    const bool dChanged = !variableD && (constantD != eigenCachedD_);
    if (!eigenFactorized_ || dtChanged || dChanged || variableD) {
      assembleAndFactorise(materials, dx, dt, Dfn);
      eigenCachedDt_ = requestedDt;
      eigenCachedD_ = variableD ? NumericType(-1) : constantD;
    }

    eigenRhs_.resize(eigenNumCells_);
    for (const auto &[ci, si] : eigenCellMap_) {
      const NumericType src = source_ ? (*source_)[ci] : NumericType(0);
      eigenRhs_[si] = field[ci] + dt * src;
    }

    eigenRhs_ = eigenSolver_.solve(eigenRhs_);
    if (eigenSolver_.info() != Eigen::Success) {
      Logger::getInstance()
          .addWarning("DiffusionSolver: Eigen solve failed.")
          .print();
      return;
    }

    for (const auto &[ci, si] : eigenCellMap_) {
      const auto v = eigenRhs_[si];
      field[ci] =
          (clampNonNegative_ && v < NumericType(0)) ? NumericType(0) : v;
    }
  }

  bool eigenMapDirty_ = true;
  bool eigenFactorized_ = false;
  NumericType eigenCachedDt_ = NumericType(-1);
  NumericType eigenCachedD_ = NumericType(-1);
  int eigenNumCells_ = 0;
  std::unordered_map<int, int> eigenCellMap_;
  Eigen::Matrix<NumericType, Eigen::Dynamic, 1> eigenRhs_;
  Eigen::SparseLU<Eigen::SparseMatrix<NumericType>> eigenSolver_;

#else
  void invalidateEigenCache() {}
#endif

  // -----------------------------------------------------------------------
  // Members
  // -----------------------------------------------------------------------
  SmartPointer<DenseCellSet<NumericType, D>> cellSet_;
  DiffusionSolverMode mode_ = DiffusionSolverMode::GaussSeidel;
  bool clampNonNegative_ = true;
  int maxIterations_ = 400;
  NumericType tolerance_ = NumericType(1e-6);
  NumericType relaxation_ = NumericType(1);
  std::unordered_set<int> diffusiveMaterials_;
  std::unordered_set<int> sourceMaterials_;
  std::unordered_set<int> blockedMaterials_;
  const std::vector<NumericType> *source_ = nullptr;
  std::vector<NumericType> buffer_;
  bool activeCacheDirty_ = true;
  const NumericType *activeMaterialsPtr_ = nullptr;
  std::size_t activeMaterialsSize_ = 0;
  std::vector<int> activeCells_;
  std::vector<std::vector<int>> activeNeighbors_;
};

} // namespace viennacs
