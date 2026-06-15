#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "csDenseCellSet.hpp"
#include "csDiffusionSolver.hpp"
#include "csBoundaryDiffusionSolver.hpp"
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
  void setSurfaceBoundaryCondition(
      const BoundaryCondition<NumericType> &condition) {
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

  void clearEquilibriumArrhenius() {
    equilibriumArrhEnabled_ = false;
  }

  void setDefectEnhancedDiffusion(const NumericType tedCoefficient,
                                  const NumericType normalization) {
    tedCoefficient_ = std::max(tedCoefficient, NumericType(0));
    tedNormalization_ = std::max(normalization, NumericType(1e-30));
  }

  void setTEDFromDamageFactor(
      const NumericType damageFactor,
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

  NumericType diffusivity() const {
    return diffusivityAt(temperatureK_);
  }

  NumericType diffusivityAt(
      const NumericType temperatureK) const {
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
      embeddedDiffSolver_.setDefaultBoundaryCondition(surfaceBoundaryCondition_);

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

        const auto diffusionCoefficient =
            diffusivityAt(localTemperatureK);
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
              diffSolver_.step(*vacancy, *materials, dx, dt, vacancyDiffusivity_);
          }

#pragma omp parallel for
          for (int i = 0; i < static_cast<int>(interstitial->size()); ++i) {
            const auto I = (*interstitial)[i];
            const auto V = (*vacancy)[i];
            const auto [Ieq, Veq] =
                defectEquilibriumAt(localTemperatureK);
            const auto recombination = recombinationRate_ * (I * V - Ieq * Veq);
            const auto sinkI = interstitialSinkRate_ * (I - Ieq);
            const auto sinkV = vacancySinkRate_ * (V - Veq);
            auto newI = I - dt * (recombination + sinkI);
            auto newV = V - dt * (recombination + sinkV);
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
            if (tedCoefficient_ > NumericType(0)) {
              auto tedFn = [&](const int idx) {
                const auto defectTerm =
                    std::max(NumericType(0),
                             ((*interstitial)[idx] - (*vacancy)[idx]) /
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
            if (useEmbedded)
              embeddedDiffSolver_.step(*concentration, *materials, dx, dt,
                                       diffusionCoefficient);
            else
              diffSolver_.step(*concentration, *materials, dx, dt,
                               diffusionCoefficient);
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
        const auto [Ieq, Veq] =
            defectEquilibriumAt(row.temperature_K);
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
    if (!activationEnabled_ || solidSolubilityC0_ <= NumericType(0))
      return;
    if (!cellSet_) return;

    auto *concentration = cellSet_->getScalarData(speciesLabel_);
    auto *materials     = cellSet_->getScalarData(materialLabel_);
    if (!concentration || !materials) return;

    auto *cellData = &cellSet_->getCellGrid()->getCellData();
    if (cellData->getScalarData(activeLabel_, true) == nullptr)
      cellSet_->addScalarData(activeLabel_, NumericType(0));
    auto *active = cellSet_->getScalarData(activeLabel_);

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
    const NumericType T =
        temperatureSchedule_.empty()
            ? temperatureK_
            : temperatureSchedule_.back().endTemperatureK;
    const auto C_SS = evalArrhenius(solidSolubilityC0_, solidSolubilityEa_eV_, T);

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(concentration->size()); ++i) {
      const auto mat = static_cast<int>((*materials)[i]);
      if (!isDiffusiveMaterial(mat) || isBlockedMaterial(mat)) {
        (*active)[i] = NumericType(0);
        continue;
      }
      const auto C = std::max((*concentration)[i], NumericType(0));
      (*active)[i] = (C_SS > NumericType(0)) ? (C_SS * C / (C_SS + C)) : C;
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
  bool defectClusteringEnabled_ = false;
  std::string clusterLabel_ = "ICluster";
  NumericType ikfi_ = NumericType(0);
  NumericType ikfc_ = NumericType(0);
  NumericType ikr_ = NumericType(0);
  NumericType clusterInitFraction_ = NumericType(0);
  bool diagnosticsEnabled_ = false;
  int diagnosticsMaterialId_ = -1;
  NumericType diagnosticsElapsedTime_ = NumericType(0);
  int diagnosticsStepCounter_ = 0;
  std::vector<DefectDiagnosticsRow> defectDiagnostics_;
};

} // namespace viennacs
