#include <csGridTraversal.hpp>

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>

#include <cmath>
#include <iostream>
#include <set>

// A cell set stores its cells compactly, in the order the level-set iterator
// found them, so a lattice coordinate says nothing about where a cell lives in
// the array. LatticeMap is the table that bridges the two, and a ray walking
// the grid trusts it once per cell it crosses. If it is wrong -- by a rounding
// convention, by a bounding box that spans nodes where cells were assumed, by
// a node that is not the corner it was taken to be -- every flux the voxel
// method computes lands in the wrong cell, quietly.
//
// So it is checked against the cell set's own BVH lookup, which is slow, well
// used and independent of it.

namespace ls = viennals;
namespace cs = viennacs;

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

/// A substrate, optionally with a trench cut into it.
///
/// NOTE on `depth`: it is an absolute COORDINATE of the cut plane, not a
/// thickness. A depth above the surface generates the gas region too and the
/// lattice comes out full -- which is what a simulation wants, so that a
/// deposit has somewhere to grow. A depth below the surface stops at that
/// plane, and a trench then leaves the lattice with holes.
template <int D>
cs::DenseCellSet<T, D> makeCellSet(T gridDelta, bool withTrench, bool withCover,
                                   T depth) {
  ls::BoundaryConditionEnum boundaryConds[D];
  for (int i = 0; i < D - 1; ++i)
    boundaryConds[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  boundaryConds[D - 1] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;

  T bounds[2 * D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -5.;
    bounds[2 * i + 1] = 5.;
  }

  T origin[D] = {};
  T normal[D] = {};
  normal[D - 1] = 1.;

  auto substrate =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryConds, gridDelta);
  ls::MakeGeometry<T, D>(
      substrate, ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
      .apply();

  if (withTrench) {
    T minCorner[D], maxCorner[D];
    for (int i = 0; i < D - 1; ++i) {
      minCorner[i] = -2.;
      maxCorner[i] = 2.;
    }
    minCorner[D - 1] = -3.;
    maxCorner[D - 1] = 1.;
    auto trench = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryConds,
                                                          gridDelta);
    ls::MakeGeometry<T, D>(
        trench, ls::SmartPointer<ls::Box<T, D>>::New(minCorner, maxCorner))
        .apply();
    ls::BooleanOperation<T, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  auto levelSets = std::vector<ls::SmartPointer<ls::Domain<T, D>>>{substrate};

  cs::DenseCellSet<T, D> cellSet;
  cellSet.setCellSetPosition(withCover);
  cellSet.setCoverMaterial(0);
  cellSet.fromLevelSets(levelSets, nullptr, depth);
  return cellSet;
}

template <int D>
void checkMapping(const std::string &label, bool withTrench, bool withCover,
                  T depth) {
  const T gridDelta = 0.5;
  auto cellSet = makeCellSet<D>(gridDelta, withTrench, withCover, depth);

  cs::LatticeMap<T, D> lattice(cellSet);

  const auto numCells = cellSet.getNumberOfCells();
  const auto &dims = lattice.dims();
  const auto &minCorner = lattice.minCorner();

  size_t sites = 1;
  for (int d = 0; d < D; ++d)
    sites *= static_cast<size_t>(dims[d]);

  // 1. every cell must map to itself, and to what the BVH says
  int roundTripFails = 0, bvhDisagreements = 0;
  for (size_t c = 0; c < numCells; ++c) {
    const auto centre = cellSet.getCellCenter(c);
    std::array<int, D> idx{};
    for (int d = 0; d < D; ++d)
      idx[d] = static_cast<int>(
          std::floor((centre[d] - minCorner[d]) / gridDelta));

    if (lattice.cellId(idx) != static_cast<int>(c))
      ++roundTripFails;
    if (lattice.cellId(idx) != cellSet.getIndex(centre))
      ++bvhDisagreements;
  }
  check(label + ": every cell round-trips through its lattice coordinate",
        roundTripFails == 0,
        std::to_string(numCells) + " cells, " + std::to_string(roundTripFails) +
            " wrong");
  check(label + ": the table agrees with the cell set's own BVH lookup",
        bvhDisagreements == 0,
        std::to_string(bvhDisagreements) + " disagreements");

  // 2. the mapping is injective: no two lattice sites claim the same cell
  std::set<int> claimed;
  size_t occupied = 0, duplicates = 0;
  std::array<int, D> idx{};
  for (size_t flat = 0; flat < sites; ++flat) {
    size_t rem = flat;
    for (int d = 0; d < D; ++d) {
      idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
      rem /= static_cast<size_t>(dims[d]);
    }
    const int id = lattice.cellId(idx);
    if (id < 0)
      continue;
    ++occupied;
    if (!claimed.insert(id).second)
      ++duplicates;
  }
  check(label + ": no cell is claimed by two lattice sites", duplicates == 0,
        std::to_string(duplicates) + " duplicates");
  check(label + ": every cell of the set appears in the table",
        occupied == numCells,
        std::to_string(occupied) + " mapped of " + std::to_string(numCells));

  // 3. a sparse set must actually have holes, or this test proves nothing
  //    about the -1 path
  const bool hasHoles = occupied < sites;
  if (depth > 0)
    check(label + ": a set spanning the gas region fills its lattice", !hasHoles,
          std::to_string(sites - occupied) + " empty sites");
  else
    check(label + ": a trench below the cut plane leaves holes (the -1 path)",
          hasHoles, std::to_string(sites - occupied) + " empty sites");

  // 4. out of range is -1, not a wrapped or clamped answer
  std::array<int, D> outside{};
  for (int d = 0; d < D; ++d)
    outside[d] = dims[d] + 3;
  bool outOfRange = lattice.cellId(outside) == -1;
  for (int d = 0; d < D; ++d)
    outside[d] = -2;
  outOfRange = outOfRange && lattice.cellId(outside) == -1;
  check(label + ": a coordinate off the lattice returns -1", outOfRange);
}

/// The two pieces together, which is how they will be used: a ray walks the
/// grid, and each cell it crosses is resolved to a cell of the set.
template <int D> void checkTraversalThroughLattice(const std::string &label) {
  const T gridDelta = 0.5;
  auto cellSet = makeCellSet<D>(gridDelta, /*withTrench*/ true,
                                /*withCover*/ false, /*depth*/ -4.);
  cs::LatticeMap<T, D> lattice(cellSet);
  cs::GridTraversal<T, D> traversal(lattice);

  // Straight down, outside the trench, where the substrate is solid.
  std::array<T, D> origin{}, dir{};
  for (int d = 0; d < D; ++d) {
    origin[d] = -3.75; // a column clear of the trench
    dir[d] = 0;
  }
  origin[D - 1] = 20.;
  dir[D - 1] = -1.;

  std::vector<int> ids;
  std::vector<T> depths;
  traversal.traverse(origin, dir, [&](cs::GridStep<T, D> s) {
    const int id = lattice.cellId(s.index);
    if (id >= 0) {
      ids.push_back(id);
      depths.push_back(cellSet.getCellCenter(id)[D - 1]);
    }
    return true;
  });

  check(label + ": a ray down a solid column resolves to real cells",
        !ids.empty(), std::to_string(ids.size()) + " cells hit");

  bool descending = true;
  for (size_t i = 1; i < depths.size(); ++i)
    if (depths[i] >= depths[i - 1])
      descending = false;
  check(label + ": they come back top to bottom, in order", descending);

  // Each resolved id must be the cell the BVH finds at the midpoint of the
  // segment the ray spent inside it. This is the check that would catch the
  // lattice and the traversal disagreeing about which cell is which.
  int mismatches = 0;
  traversal.traverse(origin, dir, [&](cs::GridStep<T, D> s) {
    const int id = lattice.cellId(s.index);
    if (id < 0)
      return true;
    viennacore::Vec3D<T> mid{0., 0., 0.};
    const T tm = 0.5 * (s.tEntry + s.tExit);
    for (int d = 0; d < D; ++d)
      mid[d] = origin[d] + dir[d] * tm;
    if (cellSet.getIndex(mid) != id)
      ++mismatches;
    return true;
  });
  check(label + ": the cell the ray is in is the cell the BVH finds there",
        mismatches == 0, std::to_string(mismatches) + " mismatches");

  // And a column through the trench must find the gap: no cells above the
  // trench floor, cells below it.
  std::vector<int> throughTrench;
  for (int d = 0; d < D - 1; ++d)
    origin[d] = 0.; // the trench centre
  traversal.traverse(origin, dir, [&](cs::GridStep<T, D> s) {
    if (lattice.cellId(s.index) >= 0)
      throughTrench.push_back(lattice.cellId(s.index));
    return true;
  });
  check(label + ": a column through the trench sees fewer cells than a solid one",
        throughTrench.size() < ids.size(),
        std::to_string(throughTrench.size()) + " vs " +
            std::to_string(ids.size()));
}

int main() {
  std::cout << "Lattice map: addressing a compact cell array by grid "
               "coordinate\n\n";

  std::cout << "1) a set spanning the gas region, whose lattice is full\n";
  checkMapping<2>("2D", /*trench*/ false, /*cover*/ true, /*depth*/ 4.);
  checkMapping<3>("3D", /*trench*/ false, /*cover*/ true, /*depth*/ 4.);

  std::cout << "\n2) a trench cut below the plane, whose lattice has holes\n";
  checkMapping<2>("2D", true, false, -4.);
  checkMapping<3>("3D", true, false, -4.);

  std::cout << "\n3) the traversal and the map together\n";
  checkTraversalThroughLattice<2>("2D");
  checkTraversalThroughLattice<3>("3D");

  std::cout << "\n";
  if (failures) {
    std::cout << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
