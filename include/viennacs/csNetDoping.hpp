#pragma once

// csNetDoping — net doping field and metallurgical junction utilities.
//
// Computes the net doping concentration (donor sum − acceptor sum) from
// existing scalar fields in a DenseCellSet and writes the result to a new
// field.  Also provides junction-depth extraction and carrier-concentration
// estimation from the charge balance equation.
//
// Typical use after implanting P and B and running Anneal with solid activation:
//
//   C++
//   ---
//   auto nd = NetDoping<double, 2>{};
//   nd.setCellSet(cellSet);
//   nd.addDonorLabel("P_active");
//   nd.addAcceptorLabel("B_active");
//   nd.apply();
//   double xj = nd.junctionDepth();         // nm, positive into substrate
//   double xj2 = nd.junctionDepths()[0];    // same, via vector API
//
//   Python (vps or vcs module)
//   --------------------------
//   nd = vps.NetDoping()
//   nd.setCellSet(domain.getCellSet())
//   nd.addDonorLabel("P_active")
//   nd.addAcceptorLabel("B_active")
//   nd.apply()
//   xj = nd.junctionDepth()                 # nm, positive into substrate

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <csDenseCellSet.hpp>

namespace viennacs {

template <typename NumericType, int D>
class NetDoping {
    using CellSet = DenseCellSet<NumericType, D>;

    SmartPointer<CellSet> cellSet_;
    std::vector<std::string> donorLabels_;
    std::vector<std::string> acceptorLabels_;
    std::string outputLabel_ = "net_doping";
    int depthAxis_ = D - 1;
    NumericType surfacePosition_ = NumericType(0);

public:
    NetDoping() = default;
    explicit NetDoping(SmartPointer<CellSet> cs) : cellSet_(std::move(cs)) {}

    // ── configuration ─────────────────────────────────────────────────────────

    void setCellSet(SmartPointer<CellSet> cs) { cellSet_ = std::move(cs); }

    // Append a single donor (n-type) concentration field.
    void addDonorLabel(const std::string &label) {
        donorLabels_.push_back(label);
    }

    // Append a single acceptor (p-type) concentration field.
    void addAcceptorLabel(const std::string &label) {
        acceptorLabels_.push_back(label);
    }

    // Replace the entire donor list at once.
    void setDonorLabels(const std::vector<std::string> &labels) {
        donorLabels_ = labels;
    }

    // Replace the entire acceptor list at once.
    void setAcceptorLabels(const std::vector<std::string> &labels) {
        acceptorLabels_ = labels;
    }

    // Cell-set field to write the result into (default: "net_doping").
    void setOutputLabel(const std::string &label) { outputLabel_ = label; }

    // Cell-centre axis index for depth (default: D−1 — y for 2-D, z for 3-D).
    void setDepthAxis(int axis) { depthAxis_ = axis; }

    // Coordinate of the wafer surface along depthAxis_. Depth is computed as
    // surfacePosition - coordinate, matching SIMS' positive-into-substrate
    // convention for ViennaPS domains with the substrate at y < 0.
    void setSurfacePosition(NumericType surfacePosition) {
        surfacePosition_ = surfacePosition;
    }

    // ── computation ───────────────────────────────────────────────────────────

    // Compute net_doping = Σ donors − Σ acceptors for every cell and write
    // the result to the field named outputLabel.  Missing labels are silently
    // skipped.
    void apply() {
        if (!cellSet_) return;

        const int n = static_cast<int>(cellSet_->getNumberOfCells());

        // Create or reset the output field first — addScalarData may reallocate
        // the internal scalar-data container, which would invalidate any pointers
        // collected before the call.
        cellSet_->addScalarData(outputLabel_, NumericType(0));

        // Collect donor/acceptor pointers after any potential reallocation.
        std::vector<const std::vector<NumericType> *> donors, acceptors;
        for (const auto &lbl : donorLabels_) {
            if (auto *d = cellSet_->getScalarData(lbl))
                donors.push_back(d);
        }
        for (const auto &lbl : acceptorLabels_) {
            if (auto *d = cellSet_->getScalarData(lbl))
                acceptors.push_back(d);
        }

        auto *out = cellSet_->getScalarData(outputLabel_);

#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            NumericType val = NumericType(0);
            for (const auto *d : donors)    val += (*d)[i];
            for (const auto *a : acceptors) val -= (*a)[i];
            (*out)[i] = val;
        }
    }

    // Return the shallowest metallurgical junction depth [domain length units,
    // positive into the substrate].
    // The junction is the depth where the laterally-averaged net_doping changes
    // sign (donor-to-acceptor or vice versa).
    // Returns infinity if no junction is found or apply() has not been called.
    NumericType junctionDepth() const {
        const auto jd = junctionDepths();
        return jd.empty() ? std::numeric_limits<NumericType>::infinity()
                          : jd.front();
    }

    // Return all junction depths [domain length units, positive into the
    // substrate], sorted ascending.
    // Useful for retrograde doping profiles that have multiple sign changes.
    std::vector<NumericType> junctionDepths() const {
        auto [depths, vals] = _depthProfile();
        std::vector<NumericType> result;
        for (std::size_t i = 0; i + 1 < depths.size(); ++i) {
            if (vals[i] * vals[i + 1] < NumericType(0)) {
                // Linear interpolation to the zero crossing.
                const NumericType t = vals[i] / (vals[i] - vals[i + 1]);
                result.push_back(depths[i] + t * (depths[i + 1] - depths[i]));
            }
        }
        return result;
    }

    // Number of sign changes (junctions) in the net_doping depth profile.
    int junctionCount() const {
        return static_cast<int>(junctionDepths().size());
    }

    // ── Lateral junction (vertical PN junction) ────────────────────────────────
    // For structures where the junction runs vertically (P on the left, B on
    // the right), the sign change is along the lateral axis, not the depth axis.
    // These methods scan along the axis perpendicular to depthAxis at the depth
    // slice nearest to atDepth.

    // Return the shallowest lateral position where net_doping changes sign
    // at the given positive substrate depth.  Returns infinity if no lateral
    // junction exists.
    NumericType lateralJunctionPosition(NumericType atDepth) const {
        const auto v = lateralJunctionPositions(atDepth);
        return v.empty() ? std::numeric_limits<NumericType>::infinity()
                         : v.front();
    }

    // All lateral positions where net_doping changes sign at the given positive
    // substrate depth.
    std::vector<NumericType> lateralJunctionPositions(NumericType atDepth) const {
        auto [xs, vals] = _lateralProfile(atDepth);
        std::vector<NumericType> result;
        for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
            if (vals[i] * vals[i + 1] < NumericType(0)) {
                const NumericType t = vals[i] / (vals[i] - vals[i + 1]);
                result.push_back(xs[i] + t * (xs[i + 1] - xs[i]));
            }
        }
        return result;
    }

private:
    // Build a depth-sorted, laterally-averaged net_doping profile.
    // Returns {depths, averagedNetDopings}.
    std::pair<std::vector<NumericType>, std::vector<NumericType>>
    _depthProfile() const {
        if (!cellSet_) return {};
        auto *nd = cellSet_->getScalarData(outputLabel_);
        if (!nd) return {};

        const NumericType delta =
            static_cast<NumericType>(cellSet_->getGridDelta());
        const int n = static_cast<int>(cellSet_->getNumberOfCells());

        // Bin by rounded positive depth, average laterally. Cells above the
        // wafer surface are excluded from junction-depth extraction.
        std::map<NumericType, std::pair<NumericType, int>> bins;
        for (int i = 0; i < n; ++i) {
            const auto center = cellSet_->getCellCenter(i);
            const NumericType z =
                surfacePosition_ - static_cast<NumericType>(center[depthAxis_]);
            if (z < NumericType(0))
                continue;
            const NumericType key =
                std::round(static_cast<double>(z / delta)) * delta;
            auto &slot = bins[key];
            slot.first  += (*nd)[i];
            slot.second += 1;
        }

        std::vector<NumericType> depths, vals;
        depths.reserve(bins.size());
        vals.reserve(bins.size());
        for (const auto &[d, sv] : bins) {
            depths.push_back(d);
            vals.push_back(sv.first / static_cast<NumericType>(sv.second));
        }
        return {depths, vals};
    }

    // Build a lateral profile at the depth slice nearest to atDepth.
    // Scans along the axis perpendicular to depthAxis, averaging the remaining
    // axes (for D > 2).  Returns {lateral_coords, averagedNetDopings}.
    std::pair<std::vector<NumericType>, std::vector<NumericType>>
    _lateralProfile(NumericType atDepth) const {
        if (!cellSet_) return {};
        auto *nd = cellSet_->getScalarData(outputLabel_);
        if (!nd) return {};

        const NumericType delta =
            static_cast<NumericType>(cellSet_->getGridDelta());
        const int n = static_cast<int>(cellSet_->getNumberOfCells());

        // Primary lateral axis: the first axis that is not depthAxis_.
        const int lateralAxis = (depthAxis_ == 0) ? 1 : 0;

        // Accept cells within half a grid cell of the requested positive depth.
        const NumericType depthKey =
            std::round(static_cast<double>(atDepth / delta)) * delta;

        std::map<NumericType, std::pair<NumericType, int>> bins;
        for (int i = 0; i < n; ++i) {
            const auto center = cellSet_->getCellCenter(i);
            const NumericType z =
                surfacePosition_ - static_cast<NumericType>(center[depthAxis_]);
            if (z < NumericType(0))
                continue;
            const NumericType zKey =
                std::round(static_cast<double>(z / delta)) * delta;
            if (std::abs(zKey - depthKey) > delta * NumericType(0.5))
                continue;
            const NumericType x = static_cast<NumericType>(center[lateralAxis]);
            const NumericType xKey =
                std::round(static_cast<double>(x / delta)) * delta;
            auto &slot = bins[xKey];
            slot.first  += (*nd)[i];
            slot.second += 1;
        }

        std::vector<NumericType> xs, vals;
        xs.reserve(bins.size());
        vals.reserve(bins.size());
        for (const auto &[x, sv] : bins) {
            xs.push_back(x);
            vals.push_back(sv.first / static_cast<NumericType>(sv.second));
        }
        return {xs, vals};
    }
};

} // namespace viennacs
