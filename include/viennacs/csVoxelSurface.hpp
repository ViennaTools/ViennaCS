#pragma once

#include "csGridTraversal.hpp"

#include <map>

namespace viennacs {

using namespace viennacore;

/// The exposed faces of a voxel geometry: every cell face with solid on one
/// side and not on the other.
///
/// This is the voxel method's surface made explicit -- a staircase, with
/// normals quantised to the 2D face directions. It exists for two reasons.
/// First, it is what a ray tracer over triangles needs, so the same geometry
/// can be put through embree and through a grid walk and the two compared.
/// Second, the face normal is one of the two candidate answers to "what normal
/// does an ion see here", the other being a gradient of the filling fractions.
/// Which of them a voxel method uses is the question the comparison is for, so
/// both have to be available.
///
/// In 2D a face is a segment, in 3D a quad split into two triangles. Faces
/// carry the cell they belong to, so a ray that hits one knows which cell it
/// struck.
template <class T, int D> struct VoxelSurface {
  std::vector<Vec3D<T>> nodes;
  std::vector<std::array<unsigned, D>> faces; ///< segments in 2D, triangles in 3D
  std::vector<int> cellId;                    ///< owning cell, per face
  std::vector<Vec3D<T>> normal;               ///< outward, per face

  size_t numFaces() const { return faces.size(); }
};

/// Extracts the exposed faces. `isSolid(cellId)` decides what counts as solid;
/// a cell id of -1 means the lattice has no cell there, which is never solid.
///
/// The default treats every cell of the set as solid, which is what a cell set
/// built without a gas region gives. A set that spans the gas region must pass
/// a predicate on the material or the filling fraction instead.
template <class T, int D, class SolidPredicate>
VoxelSurface<T, D> extractVoxelSurface(const LatticeMap<T, D> &lattice,
                                       SolidPredicate &&isSolid) {
  VoxelSurface<T, D> surface;
  const auto &dims = lattice.dims();
  const auto &min = lattice.minCorner();
  const T delta = lattice.gridDelta();

  std::map<std::array<int, D>, unsigned> nodeIndex;
  auto nodeAt = [&](std::array<int, D> corner) {
    auto found = nodeIndex.find(corner);
    if (found != nodeIndex.end())
      return found->second;
    Vec3D<T> p{0, 0, 0};
    for (int d = 0; d < D; ++d)
      p[d] = min[d] + delta * static_cast<T>(corner[d]);
    const auto idx = static_cast<unsigned>(surface.nodes.size());
    surface.nodes.push_back(p);
    nodeIndex.emplace(corner, idx);
    return idx;
  };

  auto solidAt = [&](std::array<int, D> idx) {
    const int id = lattice.cellId(idx);
    return id >= 0 && isSolid(id);
  };

  size_t sites = 1;
  for (int d = 0; d < D; ++d)
    sites *= static_cast<size_t>(dims[d]);

  std::array<int, D> idx{};
  for (size_t flat = 0; flat < sites; ++flat) {
    size_t rem = flat;
    for (int d = 0; d < D; ++d) {
      idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
      rem /= static_cast<size_t>(dims[d]);
    }
    if (!solidAt(idx))
      continue;
    const int owner = lattice.cellId(idx);

    for (int axis = 0; axis < D; ++axis) {
      for (int sign = -1; sign <= 1; sign += 2) {
        auto neighbour = idx;
        neighbour[axis] += sign;
        if (solidAt(neighbour))
          continue; // buried: not a surface

        Vec3D<T> outward{0, 0, 0};
        outward[axis] = static_cast<T>(sign);

        // The face sits on the low corner of the cell for sign -1, and one
        // cell further along the axis for sign +1.
        auto base = idx;
        if (sign > 0)
          base[axis] += 1;

        if constexpr (D == 2) {
          // A segment spanning the other axis. Wind it so that ViennaRay's
          // normal, (-dy, dx) for the direction p1 - p0, comes out outward.
          const int other = 1 - axis;
          auto a = base, b = base;
          b[other] += 1;
          const bool lowFirst = (axis == 0) ? (sign < 0) : (sign > 0);
          surface.faces.push_back(lowFirst
                                      ? std::array<unsigned, 2>{nodeAt(a), nodeAt(b)}
                                      : std::array<unsigned, 2>{nodeAt(b), nodeAt(a)});
          surface.cellId.push_back(owner);
          surface.normal.push_back(outward);
        } else {
          // A quad, split into two triangles. With u, v the two axes after
          // `axis` cyclically, (e_u, e_v, e_axis) is right-handed, so walking
          // the corners in that order gives a normal along +axis.
          const int u = (axis + 1) % 3;
          const int v = (axis + 2) % 3;
          std::array<int, 3> c0 = base, c1 = base, c2 = base, c3 = base;
          c1[u] += 1;
          c2[u] += 1;
          c2[v] += 1;
          c3[v] += 1;
          unsigned n0 = nodeAt(c0), n1 = nodeAt(c1), n2 = nodeAt(c2),
                   n3 = nodeAt(c3);
          if (sign < 0)
            std::swap(n1, n3); // reverse the winding for the far side
          surface.faces.push_back({n0, n1, n2});
          surface.cellId.push_back(owner);
          surface.normal.push_back(outward);
          surface.faces.push_back({n0, n2, n3});
          surface.cellId.push_back(owner);
          surface.normal.push_back(outward);
        }
      }
    }
  }
  return surface;
}

template <class T, int D>
VoxelSurface<T, D> extractVoxelSurface(const LatticeMap<T, D> &lattice) {
  return extractVoxelSurface<T, D>(lattice, [](int) { return true; });
}

} // namespace viennacs
