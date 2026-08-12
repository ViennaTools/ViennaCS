#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "csBoundaryDiffusionSolver.hpp"
#include "csDenseCellSet.hpp"
#include "csDiffusionSolver.hpp"
#include <vcLogger.hpp>

namespace viennacs {

using AnnealMode = DiffusionSolverMode;

template <class NumericType, int D> class Anneal {
public:
  struct TemperatureStep {
    NumericType duration = NumericType(0);
    NumericType startTemperatureK = NumericType(300);
    NumericType endTemperatureK = NumericType(300);
  };

  struct DefectDiagnosticsRow {
    int step = 0;
    NumericType time_s = NumericType(0);
    NumericType temperature_K = NumericType(300);
    NumericType I_mean = NumericType(0);
    NumericType V_mean = NumericType(0);
    NumericType I_min = NumericType(0);
    NumericType I_max = NumericType(0);
    NumericType V_min = NumericType(0);
    NumericType V_max = NumericType(0);
    NumericType I_over_Ieq_mean = NumericType(0);
    NumericType V_over_Veq_mean = NumericType(0);
    NumericType IV_over_IeqVeq_mean = NumericType(0);
    NumericType Ieq = NumericType(0);
    NumericType Veq = NumericType(0);
  };

  void setCellSet(SmartPointer<DenseCellSet<NumericType, D>> passedCellSet) {
    cellSet_ = passedCellSet;
    diffSolver_.setCellSet(passedCellSet);
    embeddedDiffSolver_.setCellSet(passedCellSet);
  }

  /// Boundary condition applied at embedded level-set surfaces when the cell
  /// set has embedded boundaries enabled. Default: zero-flux Neumann (dopant
  /// does not evaporate through the free surface).
  void
  setSurfaceBoundaryCondition(const BoundaryCondition<NumericType> &condition) {
    surfaceBoundaryCondition_ = condition;
  }

  void setSpeciesLabel(const std::string &label) { speciesLabel_ = label; }

  void setMaterialLabel(const std::string &label) { materialLabel_ = label; }

  void setDuration(const NumericType duration) {
    duration_ = std::max(duration, NumericType(0));
  }

  void setTimeStep(const NumericType dt) { timeStep_ = dt; }

  void setStabilityFactor(const NumericType factor) {
    stabilityFactor_ = std::clamp(factor, NumericType(1e-6), NumericType(1));
  }

  void setClampNonNegative(const bool enable = true) {
    clampNonNegative_ = enable;
    diffSolver_.setClampNonNegative(enable);
    embeddedDiffSolver_.setClampNonNegative(enable);
  }

  /// Per-cell volumetric source term for the dopant concentration field.
  /// S is added to ∂c/∂t = D∇²c + S each time step.
  /// The pointer must remain valid across all calls to apply().
  void setSourceField(const std::vector<NumericType> *source) {
    diffSolver_.setSourceField(source);
    embeddedDiffSolver_.setSourceField(source);
  }
  void clearSourceField() {
    diffSolver_.clearSourceField();
    embeddedDiffSolver_.clearSourceField();
    ownedSourceField_.clear();
  }

  /// Convenience overload that copies the data into internal storage so the
  /// caller does not need to manage the pointer lifetime (used by Python).
  void setSourceField(std::vector<NumericType> data) {
    ownedSourceField_ = std::move(data);
    diffSolver_.setSourceField(&ownedSourceField_);
    embeddedDiffSolver_.setSourceField(&ownedSourceField_);
  }

  void enableSolidActivation(const bool enable = true) {
    activationEnabled_ = enable;
  }

  void setSolidSolubilityArrhenius(const NumericType C0,
                                   const NumericType Ea_eV) {
    solidSolubilityC0_ = std::max(C0, NumericType(0));
    solidSolubilityEa_eV_ = std::max(Ea_eV, NumericType(0));
  }

  // ── Damage-gated (BIC/SPER) activation ─────────────────────────────────────
  // Electrically-active dopant fraction from the local implant-damage field:
  //   f_active = f_floor + (1 - f_floor)·(1 - exp(-(damage/D_amorph)^beta))
  // Amorphizing regions (damage >> D_amorph) regrow by solid-phase epitaxy →
  // full substitutional incorporation (f_active→1); sub-amorphizing
  // (crystalline) regions stay boron-interstitial-cluster limited at the floor
  // f_floor. The dose dependence is thus spatially resolved from the damage
  // field rather than supplied per condition. Composes multiplicatively with
  // the solid-solubility cap when both are enabled. Requires the damage field
  // (setDamageLabels) to be present on the cell set.
  void enableDamageActivation(const bool enable = true) {
    damageActivationEnabled_ = enable;
  }

  void setActivationFloor(const NumericType floor) {
    activationFloor_ = std::clamp(floor, NumericType(0), NumericType(1));
  }

  void setAmorphizationThreshold(const NumericType damageThreshold,
                                 const NumericType beta = NumericType(2)) {
    amorphThreshold_ = std::max(damageThreshold, NumericType(1e-30));
    amorphBeta_ = std::max(beta, NumericType(1e-6));
  }

  // Contiguous buried-amorphous-layer fill (solid-phase epitaxial regrowth).
  // An amorphizing implant forms a connected amorphous layer from the surface
  // down to the amorphous/crystalline interface. When enabled, every diffusive
  // cell shallower than the deepest amorphized cell in its lateral column is
  // treated as amorphized (full SPER activation), not just cells whose own
  // damage exceeds the threshold. Assumes the surface is toward +vertical (the
  // last coordinate axis), the ViennaPS convention.
  void enableAmorphousLayerFill(const bool enable = true) {
    amorphLayerFillEnabled_ = enable;
  }

  // ── Interface segregation of the dopant
  // ───────────────────────────────────── Dopant flux J = v_seg·C into blocking
  // materials (oxide) at material interfaces: boron segregates into SiO2
  // (segregation coefficient < 1), so the Si depletes next to the interface.
  // Applied volumetrically to the interface-adjacent cell: dC/dt =
  // −v_seg·n_faces·C/dx, integrated as an exponential decay (unconditionally
  // stable).  v_seg in [length]/s using the same length unit as the grid.
  void setInterfaceSegregation(const NumericType velocity,
                               const NumericType width = NumericType(0)) {
    segregationVelocity_ = std::max(velocity, NumericType(0));
    segWidth_ = width;
  }

  // Restrict trapping / segregation to interfaces whose blocking neighbour is
  // in the given material set (empty = all blocking materials).  Lets a thin
  // screen oxide (trap → pile-up) and a thick buried oxide (segregation → loss)
  // be treated differently even though both are SiO2 by giving the BOX its own
  // material label.
  void setTrapMaterials(const std::vector<int> &materials) {
    trapMaterials_ =
        std::unordered_set<int>(materials.begin(), materials.end());
    ifDistTrap_.clear();
  }
  void setSegregationMaterials(const std::vector<int> &materials) {
    segMaterials_ = std::unordered_set<int>(materials.begin(), materials.end());
    ifDistSeg_.clear();
  }

  // ── Interface trapping of the dopant
  // ──────────────────────────────────────── Dopant capture J = v_trap·C at
  // interface-adjacent cells into an immobile companion field (label below):
  // boron piles up at Si/SiO2 interfaces during RTA ("interfacial pile-up" /
  // dose loss in the USJ literature) — trapped, electrically inactive, but
  // still visible to SIMS.  Unlike segregation the dose is conserved:
  // SIMS-comparable profile = mobile + trapped. width ≤ 0: capture only in
  // interface-adjacent cells (surface-flux form, dC/dt = −v·n_faces·C/dx).
  // width > 0: capture distributed over a near-interface band, dC/dt =
  // −(v/w)·exp(−d/w)·C with d the distance to the nearest blocking interface —
  // same integrated capture strength, but the pile-up develops the ~w-wide
  // shoulder seen in SIMS instead of a single-cell delta.
  void setInterfaceTrap(const NumericType velocity,
                        const NumericType width = NumericType(0)) {
    trapVelocity_ = std::max(velocity, NumericType(0));
    trapWidth_ = width;
  }

  void setTrappedLabel(const std::string &label) { trappedLabel_ = label; }

  void setActiveLabel(const std::string &label) { activeLabel_ = label; }

  void setMode(const AnnealMode mode) {
    mode_ = mode;
    diffSolver_.setMode(mode);
    embeddedDiffSolver_.setMode(mode);
  }

  void setImplicitSolverOptions(const int maxIterations,
                                const NumericType relativeTolerance,
                                const NumericType relaxation = NumericType(1)) {
    diffSolver_.setMaxIterations(maxIterations);
    diffSolver_.setTolerance(relativeTolerance);
    diffSolver_.setRelaxation(relaxation);
    embeddedDiffSolver_.setMaxIterations(maxIterations);
    embeddedDiffSolver_.setTolerance(relativeTolerance);
    embeddedDiffSolver_.setRelaxation(relaxation);
  }

  void setDiffusionCoefficient(const NumericType diffusionCoefficient) {
    diffusionCoefficient_ = std::max(diffusionCoefficient, NumericType(0));
    useArrhenius_ = false;
  }

  void setArrheniusParameters(const NumericType D0, const NumericType Ea_eV) {
    D0_ = std::max(D0, NumericType(0));
    Ea_eV_ = std::max(Ea_eV, NumericType(0));
    useArrhenius_ = true;
  }

  // Fermi-level (concentration-dependent) diffusivity enhancement.
  //
  // Models B-in-Si Fermi-level effect: dopants diffusing via charged point
  // defects see an effective diffusivity that scales with local carrier
  // concentration.  For a p-type dopant at local concentration C:
  //
  //   p(z) = [C + sqrt(C^2 + 4*ni^2)] / 2       (mass action, single acceptor)
  //   h(z) = f*(p/ni) + (1-f)*(ni/p)             (Fermi enhancement factor)
  //   D_eff(z) = D_arr(T) * h(z)
  //
  // where ni = intrinsic carrier density, f = chargedFraction (fraction of
  // diffusion via the negatively-charged interstitial mechanism; ~0.9 for B).
  //
  // ni(T) is computed from the Varshni bandgap model for Si.  The parameters
  // ni_C0 and ni_Ea_eV are the pre-factor and activation (same form as
  // evalArrhenius) for a simpler exponential fit when needed; pass 0 to use
  // the built-in Si Varshni model.
  void enableFermiEnhancement(const bool enable = true) {
    fermiEnabled_ = enable;
  }

  // chargedFraction: D_{BI^-} / D_{B,intrinsic}; typical value 0.9 for B in Si.
  void setFermiChargedFraction(const NumericType chargedFraction) {
    fermiChargedFraction_ =
        std::clamp(chargedFraction, NumericType(0), NumericType(1));
  }

  void setTemperature(const NumericType temperatureK) {
    temperatureK_ = std::max(temperatureK, NumericType(1));
  }

  void clearTemperatureSchedule() { temperatureSchedule_.clear(); }

  void addIsothermalStep(const NumericType duration,
                         const NumericType temperatureK) {
    TemperatureStep step;
    step.duration = std::max(duration, NumericType(0));
    step.startTemperatureK = std::max(temperatureK, NumericType(1));
    step.endTemperatureK = step.startTemperatureK;
    temperatureSchedule_.push_back(step);
  }

  void addRampStep(const NumericType duration,
                   const NumericType startTemperatureK,
                   const NumericType endTemperatureK) {
    TemperatureStep step;
    step.duration = std::max(duration, NumericType(0));
    step.startTemperatureK = std::max(startTemperatureK, NumericType(1));
    step.endTemperatureK = std::max(endTemperatureK, NumericType(1));
    temperatureSchedule_.push_back(step);
  }

  void setTemperatureSchedule(const std::vector<NumericType> &durations,
                              const std::vector<NumericType> &temperatures) {
    if (durations.empty() || temperatures.empty())
      return;

    clearTemperatureSchedule();
    const auto n = durations.size();
    if (temperatures.size() == n + 1) {
      for (std::size_t i = 0; i < n; ++i)
        addRampStep(durations[i], temperatures[i], temperatures[i + 1]);
    } else if (temperatures.size() == n) {
      for (std::size_t i = 0; i < n; ++i)
        addIsothermalStep(durations[i], temperatures[i]);
    } else {
      Logger::getInstance()
          .addWarning("setTemperatureSchedule: temperatures count (" +
                      std::to_string(temperatures.size()) +
                      ") must equal durations count (" + std::to_string(n) +
                      ") or count+1. Schedule not changed.")
          .print();
    }
  }

  void setDiffusionMaterials(const std::vector<int> &materials) {
    diffusionMaterials_ = materials;
    diffSolver_.setDiffusiveMaterials(materials);
    embeddedDiffSolver_.setDiffusiveMaterials(materials);
  }

  void setBlockingMaterials(const std::vector<int> &materials) {
    blockingMaterials_ = materials;
    diffSolver_.setBlockedMaterials(materials);
    embeddedDiffSolver_.setBlockedMaterials(materials);
  }

  void enableDefectCoupling(const bool enable = true) {
    defectCouplingEnabled_ = enable;
  }

  void resetDefectInitialization() { defectFieldsInitialized_ = false; }

  void setDamageLabels(const std::string &damageLabel,
                       const std::string &lastDamageLabel) {
    damageLabel_ = damageLabel;
    lastDamageLabel_ = lastDamageLabel;
  }

  void setDefectLabels(const std::string &interstitialLabel,
                       const std::string &vacancyLabel) {
    interstitialLabel_ = interstitialLabel;
    vacancyLabel_ = vacancyLabel;
  }

  void setDefectSourceWeights(const NumericType historyWeight,
                              const NumericType lastDamageWeight) {
    damageHistoryWeight_ = std::max(historyWeight, NumericType(0));
    defectFromLastDamageWeight_ = std::max(lastDamageWeight, NumericType(0));
  }

  void setDefectPartition(const NumericType interstitialFraction,
                          const NumericType vacancyFraction) {
    defectToInterstitial_ = std::max(interstitialFraction, NumericType(0));
    defectToVacancy_ = std::max(vacancyFraction, NumericType(0));
  }

  void setDefectPartitionFactors(const NumericType interstitialFactor,
                                 const NumericType vacancyFactor) {
    const auto i = std::max(interstitialFactor, NumericType(0));
    const auto v = std::max(vacancyFactor, NumericType(0));
    const auto sum = i + v;
    if (sum <= NumericType(0)) {
      defectToInterstitial_ = NumericType(0.5);
      defectToVacancy_ = NumericType(0.5);
      return;
    }
    defectToInterstitial_ = i / sum;
    defectToVacancy_ = v / sum;
  }

  void setDefectDiffusivities(const NumericType interstitialDiffusivity,
                              const NumericType vacancyDiffusivity) {
    interstitialDiffusivity_ =
        std::max(interstitialDiffusivity, NumericType(0));
    vacancyDiffusivity_ = std::max(vacancyDiffusivity, NumericType(0));
  }

  void setDefectReactionRates(const NumericType recombinationRate,
                              const NumericType interstitialSinkRate,
                              const NumericType vacancySinkRate) {
    recombinationRate_ = std::max(recombinationRate, NumericType(0));
    interstitialSinkRate_ = std::max(interstitialSinkRate, NumericType(0));
    vacancySinkRate_ = std::max(vacancySinkRate, NumericType(0));
  }

  void enableDefectEquilibrium(const bool enable = true) {
    defectEquilibriumEnabled_ = enable;
  }

  void setDefectEquilibrium(const NumericType interstitialEquilibrium,
                            const NumericType vacancyEquilibrium) {
    interstitialEquilibrium_ =
        std::max(interstitialEquilibrium, NumericType(0));
    vacancyEquilibrium_ = std::max(vacancyEquilibrium, NumericType(0));
    equilibriumArrhEnabled_ = false;
  }

  void setDefectEquilibriumArrhenius(const NumericType interstitialC0,
                                     const NumericType interstitialEa_eV,
                                     const NumericType vacancyC0,
                                     const NumericType vacancyEa_eV) {
    interstitialEquilibriumC0_ = std::max(interstitialC0, NumericType(0));
    interstitialEquilibriumEa_eV_ = std::max(interstitialEa_eV, NumericType(0));
    vacancyEquilibriumC0_ = std::max(vacancyC0, NumericType(0));
    vacancyEquilibriumEa_eV_ = std::max(vacancyEa_eV, NumericType(0));
    equilibriumArrhEnabled_ = true;
  }

  void clearEquilibriumArrhenius() { equilibriumArrhEnabled_ = false; }

  void setDefectEnhancedDiffusion(const NumericType tedCoefficient,
                                  const NumericType normalization) {
    tedCoefficient_ = std::max(tedCoefficient, NumericType(0));
    tedNormalization_ = std::max(normalization, NumericType(1e-30));
  }

  // Instantaneous effective TED: D_eff = D*(1 + factor*damage/norm), driven by
  // the static implant-damage field (no dynamic defect solver). Bounded/stable.
  void setStaticDamageTED(const NumericType factor,
                          const NumericType normalization) {
    tedStaticFactor_ = std::max(factor, NumericType(0));
    tedStaticNorm_ = std::max(normalization, NumericType(1e-30));
  }

  void
  setTEDFromDamageFactor(const NumericType damageFactor,
                         const NumericType coefficientScale = NumericType(0.5),
                         const NumericType normalization = NumericType(1e20)) {
    const auto d = std::max(damageFactor, NumericType(0));
    tedCoefficient_ = std::max(coefficientScale, NumericType(0)) * d;
    tedNormalization_ = std::max(normalization, NumericType(1e-30));
  }

  void enableDefectClustering(const bool enable = true) {
    defectClusteringEnabled_ = enable;
  }

  void setDefectClusterLabel(const std::string &clusterLabel) {
    clusterLabel_ = clusterLabel;
  }

  void setDefectClusterKinetics(const NumericType ikfi, const NumericType ikfc,
                                const NumericType ikr) {
    ikfi_ = std::max(ikfi, NumericType(0));
    ikfc_ = std::max(ikfc, NumericType(0));
    ikr_ = std::max(ikr, NumericType(0));
  }

  void setDefectClusterInitFraction(const NumericType fraction) {
    clusterInitFraction_ = std::clamp(fraction, NumericType(0), NumericType(1));
  }

  void enableDiagnostics(const bool enable = true) {
    diagnosticsEnabled_ = enable;
  }

  void setDiagnosticsMaterialFilter(const int materialId = -1) {
    diagnosticsMaterialId_ = materialId;
  }

  void clearDefectDiagnostics() {
    defectDiagnostics_.clear();
    diagnosticsElapsedTime_ = NumericType(0);
    diagnosticsStepCounter_ = 0;
  }

  const std::vector<DefectDiagnosticsRow> &getDefectDiagnostics() const {
    return defectDiagnostics_;
  }

  NumericType diffusivity() const { return diffusivityAt(temperatureK_); }

  NumericType diffusivityAt(const NumericType temperatureK) const {
    if (!useArrhenius_)
      return diffusionCoefficient_;
    constexpr NumericType kB_eV_per_K = NumericType(8.617333262145e-5);
    const auto clampedT = std::max(temperatureK, NumericType(1));
    return D0_ * std::exp(-Ea_eV_ / (kB_eV_per_K * clampedT));
  }

  void apply() {
    if (!cellSet_) {
      Logger::getInstance().addWarning("No cellSet passed to Anneal.").print();
      return;
    }

    if (duration_ <= NumericType(0) && temperatureSchedule_.empty())
      return;

    auto concentration = cellSet_->getScalarData(speciesLabel_);
    if (!concentration) {
      Logger::getInstance()
          .addWarning("Anneal species data '" + speciesLabel_ + "' not found.")
          .print();
      return;
    }

    auto materials = cellSet_->getScalarData(materialLabel_);
    if (!materials) {
      Logger::getInstance()
          .addWarning("Anneal material data '" + materialLabel_ +
                      "' not found.")
          .print();
      return;
    }

    const auto diffusionMaterialsSet = std::unordered_set<int>(
        diffusionMaterials_.begin(), diffusionMaterials_.end());
    const auto blockingMaterialsSet = std::unordered_set<int>(
        blockingMaterials_.begin(), blockingMaterials_.end());

    auto isDiffusiveMaterial = [&](int mat) {
      if (diffusionMaterialsSet.empty())
        return true;
      return diffusionMaterialsSet.count(mat) > 0;
    };

    auto isBlockedMaterial = [&](int mat) {
      return blockingMaterialsSet.count(mat) > 0;
    };

    cellSet_->buildNeighborhood();

    const auto dx = std::max(cellSet_->getGridDelta(), NumericType(1e-12));

    const bool useEmbedded = cellSet_->hasEmbeddedBoundaries();
    if (useEmbedded)
      embeddedDiffSolver_.setDefaultBoundaryCondition(
          surfaceBoundaryCondition_);

    bool useDefectCoupling = defectCouplingEnabled_;
    bool useDefectClustering = defectClusteringEnabled_;
    std::vector<NumericType> *interstitial = nullptr;
    std::vector<NumericType> *vacancy = nullptr;
    std::vector<NumericType> *cluster = nullptr;
    if (useDefectCoupling) {
      auto *cellData = &cellSet_->getCellGrid()->getCellData();
      const auto hasInterstitial =
          (cellData->getScalarData(interstitialLabel_, true) != nullptr);
      const auto hasVacancy =
          (cellData->getScalarData(vacancyLabel_, true) != nullptr);
      const auto hasCluster =
          (!useDefectClustering ||
           cellData->getScalarData(clusterLabel_, true) != nullptr);

      if (!hasInterstitial) {
        cellSet_->addScalarData(interstitialLabel_, NumericType(0));
      }
      if (!hasVacancy) {
        cellSet_->addScalarData(vacancyLabel_, NumericType(0));
      }
      if (useDefectClustering && !hasCluster) {
        cellSet_->addScalarData(clusterLabel_, NumericType(0));
      }

      // Retrieve pointers only after all potential addScalarData() calls.
      interstitial = cellSet_->getScalarData(interstitialLabel_);
      vacancy = cellSet_->getScalarData(vacancyLabel_);
      cluster = useDefectClustering ? cellSet_->getScalarData(clusterLabel_)
                                    : nullptr;

      concentration = cellSet_->getScalarData(speciesLabel_);
      materials = cellSet_->getScalarData(materialLabel_);
      if (!interstitial || !vacancy || !concentration || !materials ||
          (useDefectClustering && !cluster)) {
        Logger::getInstance()
            .addWarning("Defect-coupled anneal field setup failed. Proceeding "
                        "without defect coupling.")
            .print();
        useDefectCoupling = false;
        useDefectClustering = false;
        interstitial = nullptr;
        vacancy = nullptr;
        cluster = nullptr;
      }

      if (useDefectCoupling) {
        auto damage = cellSet_->getScalarData(damageLabel_);
        auto lastDamage = cellSet_->getScalarData(lastDamageLabel_);
        if (!damage || !lastDamage) {
          Logger::getInstance()
              .addWarning(
                  "Defect-coupled anneal requested but damage datasets '" +
                  damageLabel_ + "'/'" + lastDamageLabel_ +
                  "' are missing. Proceeding without defect coupling.")
              .print();
          useDefectCoupling = false;
          useDefectClustering = false;
          interstitial = nullptr;
          vacancy = nullptr;
          cluster = nullptr;
        } else if (!defectFieldsInitialized_) {
#pragma omp parallel for
          for (int i = 0; i < static_cast<int>(interstitial->size()); ++i) {
            const auto source = damageHistoryWeight_ * (*damage)[i] +
                                defectFromLastDamageWeight_ * (*lastDamage)[i];
            (*interstitial)[i] =
                std::max(NumericType(0), defectToInterstitial_ * source);
            (*vacancy)[i] = std::max(NumericType(0), defectToVacancy_ * source);
            if (cluster != nullptr) {
              const auto trapped = clusterInitFraction_ *
                                   std::max((*interstitial)[i], NumericType(0));
              (*cluster)[i] = trapped;
              (*interstitial)[i] =
                  std::max((*interstitial)[i] - trapped, NumericType(0));
            }
          }
          defectFieldsInitialized_ = true;
        }
      }
    }

    auto integrateSegment = [&](const NumericType segmentDuration,
                                const NumericType startTemperatureK,
                                const NumericType endTemperatureK,
                                const int scheduleStep) {
      if (segmentDuration <= NumericType(0))
        return;

      NumericType time = NumericType(0);
      while (time < segmentDuration) {
        const auto remaining = segmentDuration - time;
        NumericType localTemperatureK = startTemperatureK;
        if (std::abs(endTemperatureK - startTemperatureK) > NumericType(0)) {
          const auto progress =
              std::clamp(time / std::max(segmentDuration, NumericType(1e-30)),
                         NumericType(0), NumericType(1));
          localTemperatureK = startTemperatureK +
                              (endTemperatureK - startTemperatureK) * progress;
        }

        const auto diffusionCoefficient = diffusivityAt(localTemperatureK);
        NumericType dt = timeStep_;
        if (mode_ == AnnealMode::Explicit) {
          if (diffusionCoefficient > NumericType(0)) {
            // When embedded boundaries are active the stencil uses sub-grid
            // distances that can be far smaller than dx, so use the tightest
            // effective spacing for the CFL condition.
            const NumericType effectiveDx =
                useEmbedded
                    ? embeddedDiffSolver_.getEffectiveMinBoundaryDistance(dx)
                    : dx;
            const auto maxStableDt =
                stabilityFactor_ * effectiveDx * effectiveDx /
                (NumericType(2) * static_cast<NumericType>(D) *
                 diffusionCoefficient);
            if (dt <= NumericType(0))
              dt = maxStableDt;
            else
              dt = std::min(dt, maxStableDt);
          } else if (dt <= NumericType(0)) {
            dt = remaining;
          }
        } else if (dt <= NumericType(0)) {
          dt = remaining;
        }
        dt = std::max(dt, NumericType(1e-15));
        if (dt > remaining)
          dt = remaining;
        // With the effective TED active D_eff is enhanced well above the
        // intrinsic value; over a whole-segment implicit step the diffusion
        // number D_eff*dt/dx^2 is too large for the iterative solver to
        // converge. Cap dt (the while-loop sub-steps) so each step stays
        // convergent.
        if (tedStaticFactor_ > NumericType(0) &&
            defectMaxTimeStep_ > NumericType(0))
          dt = std::min(dt, defectMaxTimeStep_);

        if (useDefectCoupling && interstitial && vacancy) {
          if (interstitialDiffusivity_ > NumericType(0)) {
            if (useEmbedded)
              embeddedDiffSolver_.step(*interstitial, *materials, dx, dt,
                                       interstitialDiffusivity_);
            else
              diffSolver_.step(*interstitial, *materials, dx, dt,
                               interstitialDiffusivity_);
          }
          if (vacancyDiffusivity_ > NumericType(0)) {
            if (useEmbedded)
              embeddedDiffSolver_.step(*vacancy, *materials, dx, dt,
                                       vacancyDiffusivity_);
            else
              diffSolver_.step(*vacancy, *materials, dx, dt,
                               vacancyDiffusivity_);
          }

#pragma omp parallel for
          for (int i = 0; i < static_cast<int>(interstitial->size()); ++i) {
            const auto I = (*interstitial)[i];
            const auto V = (*vacancy)[i];
            const auto [Ieq, Veq] = defectEquilibriumAt(localTemperatureK);
            // Pair recombination k·I·V is stiff right after implantation
            // (k·I ≫ 1/dt): explicit Euler overshoots negative and the
            // non-negative clamp then wipes out BOTH populations including
            // the net excess Δ = I−V — even though exact pair kinetics
            // conserve Δ (recombination removes pairs).  Δ is what drives
            // TED, so integrate the pair term exactly:
            //   dV/dt = −k·V·(V+Δ),  Δ const →
            //   V(t) = Δ·V₀ / ((V₀+Δ)·e^{kΔt} − V₀)   (logistic; Δ≠0)
            //   V(t) = V₀ / (1 + k·V₀·t)              (Δ=0)
            // Unconditionally stable and positivity-preserving.  The
            // equilibrium generation term k·Ieq·Veq and the sinks are gentle
            // and stay explicit.
            auto newI = I;
            auto newV = V;
            if (recombinationRate_ > NumericType(0)) {
              const auto delta = I - V;
              const auto kdt = recombinationRate_ * dt;
              if (std::abs(delta) * kdt < NumericType(1e-12)) {
                newV = V / (NumericType(1) + kdt * V);
                newI = newV + delta;
              } else {
                const auto expArg =
                    std::min(recombinationRate_ * delta * dt, NumericType(600));
                if (expArg < NumericType(-600)) {
                  // I side fully consumed: I→0, V→−Δ
                  newI = NumericType(0);
                  newV = -delta;
                } else {
                  const auto e = std::exp(expArg);
                  // den has the sign of Δ and |den| ≥ |Δ| (never zero):
                  // Δ>0 → den = I·e − V ≥ Δ;  Δ<0 → den ≤ Δ.  The quotient
                  // Δ·V/den is therefore ≥ 0 in both cases.
                  const auto den = (V + delta) * e - V;
                  newV = delta * V / den;
                  newI = newV + delta;
                }
              }
            }
            newI += dt * (recombinationRate_ * Ieq * Veq -
                          interstitialSinkRate_ * (I - Ieq));
            newV += dt * (recombinationRate_ * Ieq * Veq -
                          vacancySinkRate_ * (V - Veq));
            if (cluster != nullptr) {
              const auto C = (*cluster)[i];
              const auto capture =
                  dt * (ikfi_ * std::max(newI, NumericType(0)) +
                        ikfc_ * std::max(newI, NumericType(0)) *
                            std::max(C, NumericType(0)));
              const auto release = dt * ikr_ * std::max(C, NumericType(0));
              auto newC = C + capture - release;
              newI = newI - capture + release;
              if (clampNonNegative_)
                newC = std::max(newC, NumericType(0));
              (*cluster)[i] = newC;
            }
            if (clampNonNegative_) {
              newI = std::max(newI, NumericType(0));
              newV = std::max(newV, NumericType(0));
            }
            (*interstitial)[i] = newI;
            (*vacancy)[i] = newV;
          }

          if (diffusionCoefficient > NumericType(0)) {
            if (fermiEnabled_) {
              // Fermi-level enhancement: D_eff(z) = D(T)*h(C(z),T)
              // where h = f*(p/ni) + (1-f)*(ni/p),  p from local [dopant].
              // May be combined with TED if defect coupling is also active.
              const auto ni = intrinsicCarrierDensityAt(localTemperatureK);
              const auto f = fermiChargedFraction_;
              auto combinedFn = [&](const int idx) {
                const auto C = std::max((*concentration)[idx], NumericType(0));
                const auto p =
                    NumericType(0.5) *
                    (C + std::sqrt(C * C + NumericType(4) * ni * ni));
                const auto h = f * (p / ni) + (NumericType(1) - f) * (ni / p);
                NumericType Deff = diffusionCoefficient * h;
                if (tedCoefficient_ > NumericType(0)) {
                  const auto defectTerm = std::max(
                      NumericType(0), ((*interstitial)[idx] - (*vacancy)[idx]) /
                                          tedNormalization_);
                  Deff *= (NumericType(1) + tedCoefficient_ * defectTerm);
                }
                return Deff;
              };
              if (useEmbedded)
                embeddedDiffSolver_.stepVariable(*concentration, *materials, dx,
                                                 dt, combinedFn);
              else
                diffSolver_.stepVariable(*concentration, *materials, dx, dt,
                                         combinedFn);
            } else if (tedCoefficient_ > NumericType(0)) {
              auto tedFn = [&](const int idx) {
                const auto defectTerm = std::max(
                    NumericType(0), ((*interstitial)[idx] - (*vacancy)[idx]) /
                                        tedNormalization_);
                return diffusionCoefficient *
                       (NumericType(1) + tedCoefficient_ * defectTerm);
              };
              if (useEmbedded)
                embeddedDiffSolver_.stepVariable(*concentration, *materials, dx,
                                                 dt, tedFn);
              else
                diffSolver_.stepVariable(*concentration, *materials, dx, dt,
                                         tedFn);
            } else {
              if (useEmbedded)
                embeddedDiffSolver_.step(*concentration, *materials, dx, dt,
                                         diffusionCoefficient);
              else
                diffSolver_.step(*concentration, *materials, dx, dt,
                                 diffusionCoefficient);
            }
          }
        } else {
          if (diffusionCoefficient > NumericType(0)) {
            if (tedStaticFactor_ > NumericType(0)) {
              // Instantaneous effective TED (Cowern/Giles): a static, bounded
              // diffusivity enhancement driven by the implant-damage field,
              //   D_eff(z) = D(T) * (1 + factor * damage(z)/norm).
              // Models the NET TED-enhanced Dt (the transient is short vs. the
              // anneal), not the resolved interstitial dynamics — no I/V
              // solver, so it is stable and bounded. factor/norm set the
              // enhancement magnitude (literature: Cowern, Packan, Stolk boron
              // TED).
              auto dmg = cellSet_->getScalarData(damageLabel_);
              auto tedStaticFn = [&](const int idx) {
                const auto e =
                    dmg ? std::max(NumericType(0), (*dmg)[idx] / tedStaticNorm_)
                        : NumericType(0);
                return diffusionCoefficient *
                       (NumericType(1) + tedStaticFactor_ * e);
              };
              if (useEmbedded)
                embeddedDiffSolver_.stepVariable(*concentration, *materials, dx,
                                                 dt, tedStaticFn);
              else
                diffSolver_.stepVariable(*concentration, *materials, dx, dt,
                                         tedStaticFn);
            } else if (fermiEnabled_) {
              const auto ni = intrinsicCarrierDensityAt(localTemperatureK);
              const auto f = fermiChargedFraction_;
              auto fermiFn = [&](const int idx) {
                const auto C = std::max((*concentration)[idx], NumericType(0));
                const auto p =
                    NumericType(0.5) *
                    (C + std::sqrt(C * C + NumericType(4) * ni * ni));
                const auto h = f * (p / ni) + (NumericType(1) - f) * (ni / p);
                return diffusionCoefficient * h;
              };
              if (useEmbedded)
                embeddedDiffSolver_.stepVariable(*concentration, *materials, dx,
                                                 dt, fermiFn);
              else
                diffSolver_.stepVariable(*concentration, *materials, dx, dt,
                                         fermiFn);
            } else {
              if (useEmbedded)
                embeddedDiffSolver_.step(*concentration, *materials, dx, dt,
                                         diffusionCoefficient);
              else
                diffSolver_.step(*concentration, *materials, dx, dt,
                                 diffusionCoefficient);
            }
          }
        }

        // Interface reactions: trapping (capture into an immobile field —
        // dose-conserving pile-up, visible to SIMS, electrically inactive) and
        // segregation (loss into the blocking material — dose removed from the
        // film).  Each acts at a configurable set of blocking-neighbour
        // materials (empty = all) with an independent near-interface band of
        // rate (v/w)·exp(-d/w); w≤0 falls back to single-cell surface flux.
        const bool trapActive = trapVelocity_ > NumericType(0);
        const bool segActive = segregationVelocity_ > NumericType(0);
        if (trapActive || segActive) {
          std::vector<NumericType> *trapped = nullptr;
          if (trapActive) {
            trapped = cellSet_->getScalarData(trappedLabel_);
            if (!trapped) {
              // Create once — addScalarData ZERO-FILLS existing fields, so it
              // must never run on a field that accumulates across steps.
              cellSet_->addScalarData(trappedLabel_, NumericType(0));
              concentration = cellSet_->getScalarData(speciesLabel_);
              materials = cellSet_->getScalarData(materialLabel_);
              trapped = cellSet_->getScalarData(trappedLabel_);
            }
          }
          const int nCells = static_cast<int>(concentration->size());

          // Multi-source BFS distance (in cells) to the nearest interface whose
          // blocking neighbour is in `trig` (empty = any blocking).  Cached.
          auto buildDistance = [&](std::vector<NumericType> &dist,
                                   const std::unordered_set<int> &trig) {
            if (static_cast<int>(dist.size()) == nCells)
              return;
            dist.assign(nCells, NumericType(-1));
            std::vector<int> queue;
            queue.reserve(nCells);
            for (int i = 0; i < nCells; ++i) {
              const auto mat = static_cast<int>((*materials)[i]);
              if (!isDiffusiveMaterial(mat) || isBlockedMaterial(mat))
                continue;
              for (const auto nb : cellSet_->getNeighbors(i)) {
                if (nb < 0)
                  continue;
                const auto nbmat = static_cast<int>((*materials)[nb]);
                if (isBlockedMaterial(nbmat) &&
                    (trig.empty() || trig.count(nbmat) > 0)) {
                  dist[i] = NumericType(0);
                  queue.push_back(i);
                  break;
                }
              }
            }
            for (std::size_t q = 0; q < queue.size(); ++q) {
              const int i = queue[q];
              for (const auto nb : cellSet_->getNeighbors(i)) {
                if (nb < 0 || dist[nb] >= NumericType(0))
                  continue;
                const auto mat = static_cast<int>((*materials)[nb]);
                if (!isDiffusiveMaterial(mat) || isBlockedMaterial(mat))
                  continue;
                dist[nb] = dist[i] + NumericType(1);
                queue.push_back(nb);
              }
            }
          };
          if (trapActive)
            buildDistance(ifDistTrap_, trapMaterials_);
          if (segActive)
            buildDistance(ifDistSeg_, segMaterials_);

          // Per-cell volumetric rate [1/s] for a band reaction.
          auto bandRate =
              [&](int i, NumericType v, NumericType w,
                  const std::vector<NumericType> &dist,
                  const std::unordered_set<int> &trig) -> NumericType {
            if (w > NumericType(0)) {
              const auto d = dist[i];
              if (d < NumericType(0))
                return NumericType(0);
              return (v / w) * std::exp(-d * dx / w);
            }
            int nf = 0;
            for (const auto nb : cellSet_->getNeighbors(i)) {
              if (nb < 0)
                continue;
              const auto nbmat = static_cast<int>((*materials)[nb]);
              if (isBlockedMaterial(nbmat) &&
                  (trig.empty() || trig.count(nbmat) > 0))
                ++nf;
            }
            return nf > 0 ? v * NumericType(nf) / dx : NumericType(0);
          };

#pragma omp parallel for
          for (int i = 0; i < nCells; ++i) {
            const auto mat = static_cast<int>((*materials)[i]);
            if (!isDiffusiveMaterial(mat) || isBlockedMaterial(mat))
              continue;
            const auto rTrap = trapActive
                                   ? bandRate(i, trapVelocity_, trapWidth_,
                                              ifDistTrap_, trapMaterials_)
                                   : NumericType(0);
            const auto rSeg = segActive
                                  ? bandRate(i, segregationVelocity_, segWidth_,
                                             ifDistSeg_, segMaterials_)
                                  : NumericType(0);
            const auto rTot = rTrap + rSeg;
            if (rTot <= NumericType(0))
              continue;
            const auto decay = std::exp(-rTot * dt);
            const auto removed = (*concentration)[i] * (NumericType(1) - decay);
            (*concentration)[i] -= removed;
            if (trapped)
              (*trapped)[i] += removed * (rTrap / rTot); // seg fraction lost
          }
        }

        time += dt;
      }

      if (diagnosticsEnabled_ && useDefectCoupling && interstitial && vacancy) {
        DefectDiagnosticsRow row;
        row.step = scheduleStep;
        row.time_s = diagnosticsElapsedTime_ + segmentDuration;
        row.temperature_K =
            NumericType(0.5) * (startTemperatureK + endTemperatureK);
        const auto [Ieq, Veq] = defectEquilibriumAt(row.temperature_K);
        row.Ieq = Ieq;
        row.Veq = Veq;

        const auto nanValue = std::numeric_limits<NumericType>::quiet_NaN();
        NumericType sumI = NumericType(0);
        NumericType sumV = NumericType(0);
        NumericType sumIOver = NumericType(0);
        NumericType sumVOver = NumericType(0);
        NumericType sumIVOver = NumericType(0);
        NumericType minI = std::numeric_limits<NumericType>::max();
        NumericType maxI = NumericType(0);
        NumericType minV = std::numeric_limits<NumericType>::max();
        NumericType maxV = NumericType(0);
        NumericType count = NumericType(0);

        for (int i = 0; i < static_cast<int>(interstitial->size()); ++i) {
          if (diagnosticsMaterialId_ >= 0 && materials &&
              static_cast<int>(std::round((*materials)[i])) !=
                  diagnosticsMaterialId_) {
            continue;
          }
          const auto iVal = std::max((*interstitial)[i], NumericType(0));
          const auto vVal = std::max((*vacancy)[i], NumericType(0));
          sumI += iVal;
          sumV += vVal;
          minI = std::min(minI, iVal);
          maxI = std::max(maxI, iVal);
          minV = std::min(minV, vVal);
          maxV = std::max(maxV, vVal);
          if (row.Ieq > NumericType(0))
            sumIOver += iVal / row.Ieq;
          if (row.Veq > NumericType(0))
            sumVOver += vVal / row.Veq;
          if (row.Ieq > NumericType(0) && row.Veq > NumericType(0))
            sumIVOver += (iVal * vVal) / (row.Ieq * row.Veq);
          count += NumericType(1);
        }

        if (count > NumericType(0)) {
          row.I_mean = sumI / count;
          row.V_mean = sumV / count;
          row.I_min = minI;
          row.I_max = maxI;
          row.V_min = minV;
          row.V_max = maxV;
          row.I_over_Ieq_mean =
              row.Ieq > NumericType(0) ? (sumIOver / count) : nanValue;
          row.V_over_Veq_mean =
              row.Veq > NumericType(0) ? (sumVOver / count) : nanValue;
          row.IV_over_IeqVeq_mean =
              (row.Ieq > NumericType(0) && row.Veq > NumericType(0))
                  ? (sumIVOver / count)
                  : nanValue;
          defectDiagnostics_.push_back(row);
        }
      }
      diagnosticsElapsedTime_ += segmentDuration;
    };

    if (temperatureSchedule_.empty()) {
      integrateSegment(duration_, temperatureK_, temperatureK_,
                       ++diagnosticsStepCounter_);
    } else {
      for (const auto &step : temperatureSchedule_) {
        integrateSegment(step.duration, step.startTemperatureK,
                         step.endTemperatureK, ++diagnosticsStepCounter_);
      }
    }

    // --- Solid activation model ---
    applyActivationImpl_();
  }

public:
  // Apply only the solid-activation model (C_active = C_SS·C / (C_SS+C))
  // without running the diffusion solver.  Equivalent to Sentaurus
  // "diffuse time=0": updates the active-concentration field so that
  // SheetResistance and NetDoping see a valid activation state immediately
  // after implantation, before any thermal anneal.
  //
  // Prerequisites: setCellSet() must have been called, and
  //   enableSolidActivation(true) + setSolidSolubilityArrhenius(C0, Ea)
  //   must have been configured.
  void applyActivation() {
    if (!cellSet_) {
      Logger::getInstance()
          .addWarning("Anneal::applyActivation: no cell set attached — "
                      "call setCellSet() first.")
          .print();
      return;
    }
    applyActivationImpl_();
  }

private:
  // Shared implementation for apply() and applyActivation().
  // Re-builds the material filter locally so it can be called independently
  // of the diffusion solver context.
  void applyActivationImpl_() {
    const bool solid =
        activationEnabled_ && solidSolubilityC0_ > NumericType(0);
    const bool damageAct = damageActivationEnabled_;
    if (!solid && !damageAct)
      return;
    if (!cellSet_)
      return;

    // Guard: species and material fields must already exist.
    if (!cellSet_->getScalarData(speciesLabel_) ||
        !cellSet_->getScalarData(materialLabel_))
      return;

    // Create the active field first — addScalarData may reallocate the
    // internal scalar-data container, invalidating any pointers obtained before
    // the call.  Fetch all working pointers only after this point.
    cellSet_->addScalarData(activeLabel_, NumericType(0));

    auto *concentration = cellSet_->getScalarData(speciesLabel_);
    auto *materials = cellSet_->getScalarData(materialLabel_);
    auto *active = cellSet_->getScalarData(activeLabel_);

    // Damage field for the BIC/SPER model (optional; only when enabled).
    std::vector<NumericType> *damage = nullptr;
    if (damageAct) {
      damage = cellSet_->getScalarData(damageLabel_);
      if (!damage)
        Logger::getInstance()
            .addWarning(
                "Anneal::applyActivation: damage activation enabled but "
                "damage field '" +
                damageLabel_ +
                "' not found — crystalline cells pinned at the floor.")
            .print();
    }

    // Rebuild the material filter sets (same logic as apply()).
    const auto diffusionMaterialsSet = std::unordered_set<int>(
        diffusionMaterials_.begin(), diffusionMaterials_.end());
    const auto blockingMaterialsSet = std::unordered_set<int>(
        blockingMaterials_.begin(), blockingMaterials_.end());

    auto isDiffusiveMaterial = [&](int mat) {
      return diffusionMaterialsSet.empty() ||
             diffusionMaterialsSet.count(mat) > 0;
    };
    auto isBlockedMaterial = [&](int mat) {
      return blockingMaterialsSet.count(mat) > 0;
    };

    // Use the configured temperature (or peak temperature for a schedule).
    const NumericType T = temperatureSchedule_.empty()
                              ? temperatureK_
                              : temperatureSchedule_.back().endTemperatureK;
    const auto C_SS =
        solid ? evalArrhenius(solidSolubilityC0_, solidSolubilityEa_eV_, T)
              : NumericType(0);

    // Buried-amorphous-layer detection: per lateral column, find the deepest
    // (minimum vertical coordinate) amorphized cell. Every diffusive cell in
    // that column shallower than this a/c interface is regrown by SPER → full
    // activation. Vertical axis is the last coordinate (surface toward +).
    const bool doFill = damageAct && damage && amorphLayerFillEnabled_;
    std::unordered_map<std::int64_t, NumericType> acInterface;
    const NumericType gridDelta = cellSet_->getGridDelta();
    auto lateralKey = [&](const auto &c) -> std::int64_t {
      std::int64_t key = 0;
      for (int ax = 0; ax < D - 1; ++ax)
        key = key * std::int64_t(1000003) +
              static_cast<std::int64_t>(std::llround(c[ax] / gridDelta));
      return key;
    };
    if (doFill) {
      for (int i = 0; i < static_cast<int>(concentration->size()); ++i) {
        const auto mat = static_cast<int>((*materials)[i]);
        if (!isDiffusiveMaterial(mat) || isBlockedMaterial(mat))
          continue;
        if (std::max((*damage)[i], NumericType(0)) < amorphThreshold_)
          continue;
        const auto c = cellSet_->getCellCenter(i);
        const auto key = lateralKey(c);
        const auto v = c[D - 1];
        auto it = acInterface.find(key);
        if (it == acInterface.end() || v < it->second)
          acInterface[key] = v;
      }
    }

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(concentration->size()); ++i) {
      const auto mat = static_cast<int>((*materials)[i]);
      if (!isDiffusiveMaterial(mat) || isBlockedMaterial(mat)) {
        (*active)[i] = NumericType(0);
        continue;
      }
      const auto C = std::max((*concentration)[i], NumericType(0));
      // Solid-solubility cap (precipitation above C_SS), if enabled.
      const NumericType base =
          (solid && C_SS > NumericType(0)) ? (C_SS * C / (C_SS + C)) : C;
      // Damage-gated active fraction (BIC/SPER), if enabled.
      NumericType frac = NumericType(1);
      if (damageAct && damage) {
        bool amorphized = false;
        if (doFill) {
          const auto c = cellSet_->getCellCenter(i);
          auto it = acInterface.find(lateralKey(c));
          // Shallower than (or at) the a/c interface → within the a-layer.
          if (it != acInterface.end() &&
              c[D - 1] >= it->second - NumericType(1e-9))
            amorphized = true;
        }
        if (amorphized) {
          frac = NumericType(1); // full SPER activation
        } else {
          const auto d = std::max((*damage)[i], NumericType(0));
          const auto x = d / amorphThreshold_;
          const auto S = NumericType(1) - std::exp(-std::pow(x, amorphBeta_));
          frac = activationFloor_ + (NumericType(1) - activationFloor_) * S;
        }
      } else if (damageAct) {
        // Damage requested but field missing: pin crystalline cells at floor.
        frac = activationFloor_;
      }
      (*active)[i] = frac * base;
    }
  }

private:
  std::pair<NumericType, NumericType>
  defectEquilibriumAt(const NumericType temperatureK) const {
    if (!defectEquilibriumEnabled_)
      return {NumericType(0), NumericType(0)};
    if (!equilibriumArrhEnabled_)
      return {interstitialEquilibrium_, vacancyEquilibrium_};
    return {evalArrhenius(interstitialEquilibriumC0_,
                          interstitialEquilibriumEa_eV_, temperatureK),
            evalArrhenius(vacancyEquilibriumC0_, vacancyEquilibriumEa_eV_,
                          temperatureK)};
  }

  NumericType evalArrhenius(const NumericType prefactor,
                            const NumericType activation_eV,
                            const NumericType temperatureK) const {
    constexpr NumericType kB_eV_per_K = NumericType(8.617333262145e-5);
    const auto clampedT = std::max(temperatureK, NumericType(1));
    return std::max(prefactor, NumericType(0)) *
           std::exp(-std::max(activation_eV, NumericType(0)) /
                    (kB_eV_per_K * clampedT));
  }

  // Intrinsic carrier concentration of Si [same units as dopant concentration].
  // Uses Varshni bandgap + Sze effective DOS.  Valid 300–1600 K.
  // At 300 K: ~9.7e9 cm^-3.  At 1303 K (1030 °C): ~1.7e19 cm^-3.
  NumericType intrinsicCarrierDensityAt(const NumericType T_K) const {
    constexpr NumericType kB_eV = NumericType(8.617333262145e-5);
    const NumericType t = T_K / NumericType(300);
    const NumericType Nc = NumericType(2.86e19) * std::pow(t, NumericType(1.5));
    const NumericType Nv = NumericType(3.10e19) * std::pow(t, NumericType(1.5));
    const NumericType Eg = NumericType(1.17) - NumericType(4.37e-4) * T_K *
                                                   T_K /
                                                   (T_K + NumericType(636));
    const NumericType kBT = kB_eV * T_K;
    return std::sqrt(Nc * Nv) * std::exp(-Eg / (NumericType(2) * kBT));
  }

  SmartPointer<DenseCellSet<NumericType, D>> cellSet_;
  std::string speciesLabel_ = "concentration";
  std::string materialLabel_ = "Material";
  NumericType duration_ = NumericType(0);
  NumericType timeStep_ = NumericType(-1);
  NumericType stabilityFactor_ = NumericType(0.45);
  bool clampNonNegative_ = true;
  bool activationEnabled_ = false;
  NumericType solidSolubilityC0_ = NumericType(0);
  NumericType solidSolubilityEa_eV_ = NumericType(0);
  std::string activeLabel_ = "active_concentration";
  AnnealMode mode_ = DiffusionSolverMode::GaussSeidel;
  DiffusionSolver<NumericType, D> diffSolver_;
  BoundaryDiffusionSolver<NumericType, D> embeddedDiffSolver_;
  BoundaryCondition<NumericType> surfaceBoundaryCondition_ =
      BoundaryCondition<NumericType>::neumann(NumericType(0));
  std::vector<NumericType> ownedSourceField_;

  bool useArrhenius_ = false;
  NumericType diffusionCoefficient_ = NumericType(0);
  NumericType D0_ = NumericType(0);
  NumericType Ea_eV_ = NumericType(0);
  NumericType temperatureK_ = NumericType(300);
  std::vector<TemperatureStep> temperatureSchedule_;

  std::vector<int> diffusionMaterials_;
  std::vector<int> blockingMaterials_;

  bool defectCouplingEnabled_ = false;
  bool defectFieldsInitialized_ = false;
  std::string damageLabel_ = "Damage";
  std::string lastDamageLabel_ = "Damage_Last";
  std::string interstitialLabel_ = "Interstitial";
  std::string vacancyLabel_ = "Vacancy";
  NumericType damageHistoryWeight_ = NumericType(0);
  NumericType defectFromLastDamageWeight_ = NumericType(1);
  NumericType defectToInterstitial_ = NumericType(0.5);
  NumericType defectToVacancy_ = NumericType(0.5);
  NumericType interstitialDiffusivity_ = NumericType(0);
  NumericType vacancyDiffusivity_ = NumericType(0);
  NumericType recombinationRate_ = NumericType(0);
  NumericType interstitialSinkRate_ = NumericType(0);
  NumericType vacancySinkRate_ = NumericType(0);
  bool defectEquilibriumEnabled_ = false;
  NumericType interstitialEquilibrium_ = NumericType(0);
  NumericType vacancyEquilibrium_ = NumericType(0);
  bool equilibriumArrhEnabled_ = false;
  NumericType interstitialEquilibriumC0_ = NumericType(0);
  NumericType interstitialEquilibriumEa_eV_ = NumericType(0);
  NumericType vacancyEquilibriumC0_ = NumericType(0);
  NumericType vacancyEquilibriumEa_eV_ = NumericType(0);
  NumericType tedCoefficient_ = NumericType(0);
  NumericType tedNormalization_ = NumericType(1e20);
  NumericType defectMaxTimeStep_ =
      NumericType(0.05);                         // s; caps dt when TED active
  NumericType tedStaticFactor_ = NumericType(0); // instantaneous effective TED
  NumericType tedStaticNorm_ = NumericType(1);
  bool fermiEnabled_ = false;
  NumericType fermiChargedFraction_ = NumericType(0.9);
  bool defectClusteringEnabled_ = false;
  std::string clusterLabel_ = "ICluster";
  NumericType ikfi_ = NumericType(0);
  NumericType ikfc_ = NumericType(0);
  NumericType ikr_ = NumericType(0);
  NumericType clusterInitFraction_ = NumericType(0);
  bool damageActivationEnabled_ = false;
  NumericType activationFloor_ = NumericType(0);
  NumericType amorphThreshold_ = NumericType(1);
  NumericType amorphBeta_ = NumericType(2);
  bool amorphLayerFillEnabled_ = true;
  NumericType segregationVelocity_ = NumericType(0);
  NumericType trapVelocity_ = NumericType(0);
  NumericType trapWidth_ = NumericType(0);
  NumericType segWidth_ = NumericType(0);
  std::unordered_set<int> trapMaterials_; // empty = all blocking
  std::unordered_set<int> segMaterials_;  // empty = all blocking
  std::vector<NumericType> ifDistTrap_;
  std::vector<NumericType> ifDistSeg_;
  std::string trappedLabel_ = "trapped_concentration";
  bool diagnosticsEnabled_ = false;
  int diagnosticsMaterialId_ = -1;
  NumericType diagnosticsElapsedTime_ = NumericType(0);
  int diagnosticsStepCounter_ = 0;
  std::vector<DefectDiagnosticsRow> defectDiagnostics_;
};

} // namespace viennacs
