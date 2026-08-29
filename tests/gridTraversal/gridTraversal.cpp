#include <csGridTraversal.hpp>

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

// A ray walking a voxel grid has to visit exactly the cells it passes through,
// in order, with no gaps and no repeats. That is checkable without a simulator
// and without a ray tracer, so it is checked here: every later comparison
// against embree, and every profile the voxel method produces, rests on it.

using T = double;

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
  if (!detail.empty())
    std::cout << "   " << detail;
  std::cout << "\n";
  if (!ok)
    ++failures;
}

template <int D> std::array<T, D> normalize(std::array<T, D> v) {
  T n = 0;
  for (int d = 0; d < D; ++d)
    n += v[d] * v[d];
  n = std::sqrt(n);
  for (int d = 0; d < D; ++d)
    v[d] /= n;
  return v;
}

// The reference: march the ray in very small steps and record each distinct
// cell. Slow, obviously correct, and independent of the algorithm under test.
template <int D>
std::vector<std::array<int, D>>
bruteForce(const std::array<T, D> &origin, const std::array<T, D> &dir,
           const std::array<T, D> &min, const std::array<int, D> &dims,
           T gridDelta) {
  std::vector<std::array<int, D>> cells;
  const T ds = gridDelta * 1e-4;
  const T tEnd =
      gridDelta * 4 * D * (*std::max_element(dims.begin(), dims.end()));
  for (T t = 0; t < tEnd; t += ds) {
    std::array<int, D> idx{};
    bool inside = true;
    for (int d = 0; d < D; ++d) {
      const T x = origin[d] + dir[d] * t;
      idx[d] = static_cast<int>(std::floor((x - min[d]) / gridDelta));
      if (idx[d] < 0 || idx[d] >= dims[d])
        inside = false;
    }
    if (!inside)
      continue;
    if (cells.empty() || cells.back() != idx)
      cells.push_back(idx);
  }
  return cells;
}

template <int D>
void compareWithBruteForce(const std::string &label, int numRays, unsigned seed,
                           T gridDelta) {
  std::array<int, D> dims{};
  std::array<T, D> min{};
  for (int d = 0; d < D; ++d) {
    dims[d] = 12;
    min[d] = -6.0 * gridDelta;
  }

  viennacs::GridTraversal<T, D> traversal(min, dims, gridDelta);

  std::mt19937 rng(seed);
  std::uniform_real_distribution<T> uni(-1.0, 1.0);

  int mismatches = 0, emptyWalks = 0;
  T worstGap = 0;
  for (int r = 0; r < numRays; ++r) {
    std::array<T, D> dir{}, origin{};
    for (int d = 0; d < D; ++d)
      dir[d] = uni(rng);
    T norm = 0;
    for (int d = 0; d < D; ++d)
      norm += dir[d] * dir[d];
    if (norm < 1e-6) {
      --r;
      continue;
    }
    dir = normalize<D>(dir);
    // Start well outside, aimed through the box.
    for (int d = 0; d < D; ++d)
      origin[d] = -dir[d] * 20.0 * gridDelta + uni(rng) * 2.0 * gridDelta;

    std::vector<std::array<int, D>> walked;
    T prevExit = -1;
    bool contiguous = true, midpointInside = true, boundedSteps = true;
    const T maxChord = gridDelta * std::sqrt(T(D)) * (1 + 1e-9);
    traversal.traverse(origin, dir, [&](viennacs::GridStep<T, D> s) {
      walked.push_back(s.index);
      if (prevExit >= 0 && std::abs(s.tEntry - prevExit) > 1e-9)
        contiguous = false;
      prevExit = s.tExit;
      const T len = s.tExit - s.tEntry;
      if (len <= 0 || len > maxChord)
        boundedSteps = false;
      // The midpoint of the segment inside a cell must lie in that cell. This
      // is the invariant that catches an off-by-one in the step or the initial
      // index, which a length check alone would miss.
      const T tm = 0.5 * (s.tEntry + s.tExit);
      for (int d = 0; d < D; ++d) {
        const T x = origin[d] + dir[d] * tm;
        const int i = static_cast<int>(std::floor((x - min[d]) / gridDelta));
        if (i != s.index[d])
          midpointInside = false;
      }
      return true;
    });

    if (walked.empty()) {
      ++emptyWalks;
      continue;
    }
    if (!contiguous || !midpointInside || !boundedSteps) {
      ++mismatches;
      continue;
    }
    auto reference = bruteForce<D>(origin, dir, min, dims, gridDelta);
    if (walked != reference)
      ++mismatches;
    worstGap = std::max(worstGap, T(0));
  }

  check(label + ": cell sequence matches a fine-stepped reference",
        mismatches == 0,
        std::to_string(numRays - emptyWalks) + " rays through the grid, " +
            std::to_string(mismatches) + " mismatched");
}

template <int D> void checkPathLength(const std::string &label) {
  const T gridDelta = 0.5;
  std::array<int, D> dims{};
  std::array<T, D> min{};
  for (int d = 0; d < D; ++d) {
    dims[d] = 20;
    min[d] = 0.0;
  }
  viennacs::GridTraversal<T, D> traversal(min, dims, gridDelta);

  std::mt19937 rng(7);
  std::uniform_real_distribution<T> uni(-1.0, 1.0);
  T worst = 0;
  for (int r = 0; r < 500; ++r) {
    std::array<T, D> dir{}, origin{};
    for (int d = 0; d < D; ++d)
      dir[d] = uni(rng);
    T norm = 0;
    for (int d = 0; d < D; ++d)
      norm += dir[d] * dir[d];
    if (norm < 1e-6) {
      --r;
      continue;
    }
    dir = normalize<D>(dir);
    for (int d = 0; d < D; ++d)
      origin[d] = 5.0 - dir[d] * 40.0;

    T tNear, tFar;
    if (!traversal.clip(origin, dir, tNear, tFar))
      continue;
    T summed = 0;
    traversal.traverse(origin, dir, [&](viennacs::GridStep<T, D> s) {
      summed += s.tExit - s.tEntry;
      return true;
    });
    worst = std::max(worst, std::abs(summed - (tFar - tNear)));
  }
  // NOTE: this telescopes by construction -- each entry is the previous exit
  // -- so it verifies the CLIPPING against the grid bounds, not the walk. The
  // per-step bound in the reference comparison is what constrains the walk.
  check(label + ": the walk spans exactly the clipped chord", worst < 1e-9,
        "worst deviation " + std::to_string(worst));
}

template <int D> void checkAxisAligned(const std::string &label) {
  // A direction with an exactly zero component is where a traversal divides by
  // zero if it is careless, and it is also the commonest ray in a process
  // simulation: straight down from the source.
  const T gridDelta = 1.0;
  std::array<int, D> dims{};
  std::array<T, D> min{};
  for (int d = 0; d < D; ++d) {
    dims[d] = 8;
    min[d] = 0.0;
  }
  viennacs::GridTraversal<T, D> traversal(min, dims, gridDelta);

  std::array<T, D> dir{}, origin{};
  for (int d = 0; d < D; ++d) {
    dir[d] = 0;
    origin[d] = 3.5;
  }
  dir[D - 1] = -1.0;
  origin[D - 1] = 20.0;

  std::vector<int> depths;
  traversal.traverse(origin, dir, [&](viennacs::GridStep<T, D> s) {
    depths.push_back(s.index[D - 1]);
    return true;
  });

  bool descending = depths.size() == 8;
  for (size_t i = 0; i < depths.size(); ++i)
    if (depths[i] != static_cast<int>(7 - i))
      descending = false;
  check(label + ": a straight-down ray visits every cell of its column once",
        descending, std::to_string(depths.size()) + " cells visited");
}

template <int D> void checkMiss(const std::string &label) {
  const T gridDelta = 1.0;
  std::array<int, D> dims{};
  std::array<T, D> min{};
  for (int d = 0; d < D; ++d) {
    dims[d] = 4;
    min[d] = 0.0;
  }
  viennacs::GridTraversal<T, D> traversal(min, dims, gridDelta);

  std::array<T, D> dir{}, origin{};
  for (int d = 0; d < D; ++d) {
    dir[d] = 0;
    origin[d] = -10.0; // well off to the side
  }
  dir[D - 1] = -1.0;
  origin[D - 1] = 20.0;

  size_t visited = traversal.traverse(
      origin, dir, [](viennacs::GridStep<T, D>) { return true; });
  check(label + ": a ray that misses the grid visits nothing", visited == 0,
        std::to_string(visited) + " cells");
}

template <int D> void checkEarlyStop(const std::string &label) {
  // An interaction stops the walk. The visitor returning false must end it
  // immediately -- a voxel method relies on this to find the FIRST cell that
  // interacts, not all of them.
  const T gridDelta = 1.0;
  std::array<int, D> dims{};
  std::array<T, D> min{};
  for (int d = 0; d < D; ++d) {
    dims[d] = 10;
    min[d] = 0.0;
  }
  viennacs::GridTraversal<T, D> traversal(min, dims, gridDelta);

  std::array<T, D> dir{}, origin{};
  for (int d = 0; d < D; ++d) {
    dir[d] = 0;
    origin[d] = 4.5;
  }
  dir[D - 1] = -1.0;
  origin[D - 1] = 20.0;

  int seen = 0;
  size_t visited =
      traversal.traverse(origin, dir, [&](viennacs::GridStep<T, D>) {
        ++seen;
        return seen < 3;
      });
  check(label + ": the walk stops when the visitor says so",
        visited == 3 && seen == 3, std::to_string(visited) + " cells visited");
}

int main() {
  std::cout << "Grid traversal: a ray walking the cells of a voxel grid\n\n";

  std::cout << "1) against a fine-stepped reference\n";
  compareWithBruteForce<2>("2D, delta 1.0 ", 400, 1, 1.0);
  compareWithBruteForce<3>("3D, delta 1.0 ", 200, 2, 1.0);
  // A spacing other than one, or a tDelta that has lost its factor of
  // gridDelta walks the grid at the wrong rate undetected.
  compareWithBruteForce<2>("2D, delta 0.25", 400, 3, 0.25);
  compareWithBruteForce<3>("3D, delta 2.50", 200, 4, 2.5);

  std::cout << "\n2) geometry the walk must conserve\n";
  checkPathLength<2>("2D");
  checkPathLength<3>("3D");

  std::cout << "\n3) the awkward cases\n";
  checkAxisAligned<2>("2D");
  checkAxisAligned<3>("3D");
  checkMiss<2>("2D");
  checkMiss<3>("3D");
  checkEarlyStop<2>("2D");
  checkEarlyStop<3>("3D");

  std::cout << "\n";
  if (failures) {
    std::cout << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
