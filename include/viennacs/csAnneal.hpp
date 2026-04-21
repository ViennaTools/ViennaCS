#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_set>
#include <vector>

#include "csDenseCellSet.hpp"
#include <vcLogger.hpp>

namespace viennacs {

enum class AnnealMode {
  Explicit,
  Implicit,
};

template <class NumericType, int D> class Anneal {
public:
  struct TemperatureStep {
    NumericType duration = NumericType(0);
    NumericType startTemperatureK = NumericType(300);
    NumericType endTemperatureK = NumericType(300);
  };

  void setCellSet(SmartPointer<DenseCellSet<NumericType, D>> passedCellSet) {
    cellSet_ = passedCellSet;
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
  }

  void setMode(const AnnealMode mode) { mode_ = mode; }

  void setImplicitSolverOptions(const int maxIterations,
                                const NumericType relativeTolerance) {
    implicitMaxIterations_ = std::max(maxIterations, 1);
    implicitRelativeTolerance_ =
        std::max(relativeTolerance, NumericType(1e-12));
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

  void addRampStep(const NumericType duration, const NumericType startTemperatureK,
                   const NumericType endTemperatureK) {
    TemperatureStep step;
    step.duration = std::max(duration, NumericType(0));
    step.startTemperatureK = std::max(startTemperatureK, NumericType(1));
    step.endTemperatureK = std::max(endTemperatureK, NumericType(1));
    temperatureSchedule_.push_back(step);
  }

  void setTemperatureRamp(const std::vector<NumericType> &temperaturesK,
                          const std::vector<NumericType> &durations) {
    clearTemperatureSchedule();
    if (durations.empty())
      return;

    if (temperaturesK.size() == durations.size() + 1) {
      for (std::size_t i = 0; i < durations.size(); ++i) {
        addRampStep(durations[i], temperaturesK[i], temperaturesK[i + 1]);
      }
      return;
    }

    if (temperaturesK.size() == durations.size()) {
      for (std::size_t i = 0; i < durations.size(); ++i) {
        addIsothermalStep(durations[i], temperaturesK[i]);
      }
      return;
    }

    Logger::getInstance()
        .addWarning(
            "Anneal temperature ramp ignored: expected temperatures size "
            "to match stepDurations size (isothermal) or size+1 (ramps).")
        .print();
  }

  void setDiffusionMaterials(const std::vector<int> &materials) {
    diffusionMaterials_ = materials;
  }

  void setBlockingMaterials(const std::vector<int> &materials) {
    blockingMaterials_ = materials;
  }

  void enableDefectCoupling(const bool enable = true) {
    defectCouplingEnabled_ = enable;
  }

  void resetDefectInitialization() { defectFieldsInitialized_ = false; }

  void setDamageLabels(const std::string &damageLabel,
                       const std::string &damageLastImpLabel) {
    damageLabel_ = damageLabel;
    damageLastImpLabel_ = damageLastImpLabel;
  }

  void setDefectLabels(const std::string &interstitialLabel,
                       const std::string &vacancyLabel) {
    interstitialLabel_ = interstitialLabel;
    vacancyLabel_ = vacancyLabel;
  }

  void setDefectSourceWeights(const NumericType historyWeight,
                              const NumericType lastImpWeight) {
    defectFromDamageHistoryWeight_ = std::max(historyWeight, NumericType(0));
    defectFromDamageLastImpWeight_ = std::max(lastImpWeight, NumericType(0));
  }

  void setDefectPartition(const NumericType interstitialFraction,
                          const NumericType vacancyFraction) {
    defectToInterstitial_ = std::max(interstitialFraction, NumericType(0));
    defectToVacancy_ = std::max(vacancyFraction, NumericType(0));
  }

  void setDefectDiffusivities(const NumericType interstitialDiffusivity,
                              const NumericType vacancyDiffusivity) {
    interstitialDiffusivity_ = std::max(interstitialDiffusivity, NumericType(0));
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
  }

  void setDefectEnhancedDiffusion(const NumericType tedCoefficient,
                                  const NumericType normalization) {
    tedCoefficient_ = std::max(tedCoefficient, NumericType(0));
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

  NumericType getEffectiveDiffusionCoefficient() const {
    return getEffectiveDiffusionCoefficientAtTemperature(temperatureK_);
  }

  NumericType getEffectiveDiffusionCoefficientAtTemperature(
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
          .addWarning("Anneal material data '" + materialLabel_ + "' not found.")
          .print();
      return;
    }

    const auto diffusionMaterialsSet =
        std::unordered_set<int>(diffusionMaterials_.begin(),
                                diffusionMaterials_.end());
    const auto blockingMaterialsSet =
        std::unordered_set<int>(blockingMaterials_.begin(),
                                blockingMaterials_.end());

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

    std::vector<NumericType> *interstitial = nullptr;
    std::vector<NumericType> *vacancy = nullptr;
    std::vector<NumericType> *cluster = nullptr;
    if (defectCouplingEnabled_) {
      auto *cellData = &cellSet_->getCellGrid()->getCellData();
      const auto hasInterstitial =
          (cellData->getScalarData(interstitialLabel_, true) != nullptr);
      const auto hasVacancy =
          (cellData->getScalarData(vacancyLabel_, true) != nullptr);
      const auto hasCluster =
          (!defectClusteringEnabled_ ||
           cellData->getScalarData(clusterLabel_, true) != nullptr);

      if (!hasInterstitial) {
        cellSet_->addScalarData(interstitialLabel_, NumericType(0));
      }
      if (!hasVacancy) {
        cellSet_->addScalarData(vacancyLabel_, NumericType(0));
      }
      if (defectClusteringEnabled_ && !hasCluster) {
        cellSet_->addScalarData(clusterLabel_, NumericType(0));
      }

      // Retrieve pointers only after all potential addScalarData() calls.
      interstitial = cellSet_->getScalarData(interstitialLabel_);
      vacancy = cellSet_->getScalarData(vacancyLabel_);
      cluster =
          defectClusteringEnabled_ ? cellSet_->getScalarData(clusterLabel_) : nullptr;

      concentration = cellSet_->getScalarData(speciesLabel_);
      materials = cellSet_->getScalarData(materialLabel_);
      if (!interstitial || !vacancy || !concentration || !materials ||
          (defectClusteringEnabled_ && !cluster)) {
        Logger::getInstance()
            .addWarning(
                "Defect-coupled anneal field setup failed. Proceeding without defect coupling.")
            .print();
        defectCouplingEnabled_ = false;
        defectClusteringEnabled_ = false;
        interstitial = nullptr;
        vacancy = nullptr;
        cluster = nullptr;
      }

      auto damage = cellSet_->getScalarData(damageLabel_);
      auto damageLastImp = cellSet_->getScalarData(damageLastImpLabel_);
      if (!damage || !damageLastImp) {
        Logger::getInstance()
            .addWarning("Defect-coupled anneal requested but damage datasets '" +
                        damageLabel_ + "'/'" + damageLastImpLabel_ +
                        "' are missing. Proceeding without defect coupling.")
            .print();
        defectCouplingEnabled_ = false;
        defectClusteringEnabled_ = false;
        interstitial = nullptr;
        vacancy = nullptr;
        cluster = nullptr;
      } else if (!defectFieldsInitialized_) {
#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(interstitial->size()); ++i) {
          const auto source =
              defectFromDamageHistoryWeight_ * (*damage)[i] +
              defectFromDamageLastImpWeight_ * (*damageLastImp)[i];
          (*interstitial)[i] =
              std::max(NumericType(0), defectToInterstitial_ * source);
          (*vacancy)[i] = std::max(NumericType(0), defectToVacancy_ * source);
          if (cluster != nullptr) {
            const auto trapped =
                clusterInitFraction_ * std::max((*interstitial)[i], NumericType(0));
            (*cluster)[i] = trapped;
            (*interstitial)[i] =
                std::max((*interstitial)[i] - trapped, NumericType(0));
          }
        }
        defectFieldsInitialized_ = true;
      }
    }

    std::vector<NumericType> next(concentration->size(), NumericType(0));
    auto integrateSegment = [&](const NumericType segmentDuration,
                                const NumericType startTemperatureK,
                                const NumericType endTemperatureK) {
      if (segmentDuration <= NumericType(0))
        return;

      NumericType time = NumericType(0);
      while (time < segmentDuration) {
        const auto remaining = segmentDuration - time;
        NumericType localTemperatureK = startTemperatureK;
        if (useArrhenius_ &&
            std::abs(endTemperatureK - startTemperatureK) > NumericType(0)) {
          const auto progress =
              std::clamp(time / std::max(segmentDuration, NumericType(1e-30)),
                         NumericType(0), NumericType(1));
          localTemperatureK = startTemperatureK +
                              (endTemperatureK - startTemperatureK) * progress;
        }

        const auto diffusionCoefficient =
            getEffectiveDiffusionCoefficientAtTemperature(localTemperatureK);
        NumericType dt = timeStep_;
        if (mode_ == AnnealMode::Explicit) {
          if (diffusionCoefficient > NumericType(0)) {
            const auto maxStableDt =
                stabilityFactor_ * dx * dx /
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

        if (defectCouplingEnabled_ && interstitial && vacancy) {
          stepDiffusion(*interstitial, next, *materials, dx, dt,
                        interstitialDiffusivity_, isDiffusiveMaterial,
                        isBlockedMaterial);
          interstitial->swap(next);
          stepDiffusion(*vacancy, next, *materials, dx, dt, vacancyDiffusivity_,
                        isDiffusiveMaterial, isBlockedMaterial);
          vacancy->swap(next);

#pragma omp parallel for
          for (int i = 0; i < static_cast<int>(interstitial->size()); ++i) {
            const auto I = (*interstitial)[i];
            const auto V = (*vacancy)[i];
            const auto Ieq =
                defectEquilibriumEnabled_ ? interstitialEquilibrium_ : NumericType(0);
            const auto Veq =
                defectEquilibriumEnabled_ ? vacancyEquilibrium_ : NumericType(0);
            const auto recombination =
                recombinationRate_ * (I * V - Ieq * Veq);
            const auto sinkI =
                interstitialSinkRate_ * (I - (defectEquilibriumEnabled_ ? Ieq : NumericType(0)));
            const auto sinkV =
                vacancySinkRate_ * (V - (defectEquilibriumEnabled_ ? Veq : NumericType(0)));
            auto newI = I - dt * (recombination + sinkI);
            auto newV = V - dt * (recombination + sinkV);
            if (cluster != nullptr) {
              const auto C = (*cluster)[i];
              const auto capture = dt * (ikfi_ * std::max(newI, NumericType(0)) +
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

          stepDiffusionVariable(
              *concentration, next, *materials, dx, dt,
              [&](const int idx) {
                const auto defectTerm =
                    std::max(NumericType(0),
                             ((*interstitial)[idx] - (*vacancy)[idx]) /
                                 tedNormalization_);
                return diffusionCoefficient *
                       (NumericType(1) + tedCoefficient_ * defectTerm);
              },
              isDiffusiveMaterial, isBlockedMaterial);
        } else {
          stepDiffusion(*concentration, next, *materials, dx, dt,
                        diffusionCoefficient, isDiffusiveMaterial,
                        isBlockedMaterial);
        }

        concentration->swap(next);
        if (clampNonNegative_) {
#pragma omp parallel for
          for (int i = 0; i < static_cast<int>(concentration->size()); ++i) {
            if ((*concentration)[i] < NumericType(0))
              (*concentration)[i] = NumericType(0);
          }
        }
        time += dt;
      }
    };

    if (temperatureSchedule_.empty()) {
      integrateSegment(duration_, temperatureK_, temperatureK_);
    } else {
      for (const auto &step : temperatureSchedule_) {
        integrateSegment(step.duration, step.startTemperatureK,
                         step.endTemperatureK);
      }
    }
  }

private:
  template <typename IsDiffusiveFn, typename IsBlockedFn>
  void stepDiffusion(const std::vector<NumericType> &field,
                     std::vector<NumericType> &next,
                     const std::vector<NumericType> &materials, NumericType dx,
                     NumericType dt, NumericType diffusivity,
                     const IsDiffusiveFn &isDiffusiveMaterial,
                     const IsBlockedFn &isBlockedMaterial) {
    stepDiffusionVariable(
        field, next, materials, dx, dt,
        [&](const int) { return diffusivity; }, isDiffusiveMaterial,
        isBlockedMaterial);
  }

  template <typename DiffusivityAtIndexFn, typename IsDiffusiveFn,
            typename IsBlockedFn>
  void stepDiffusionVariable(const std::vector<NumericType> &field,
                             std::vector<NumericType> &next,
                             const std::vector<NumericType> &materials,
                             NumericType dx, NumericType dt,
                             const DiffusivityAtIndexFn &diffusivityAt,
                             const IsDiffusiveFn &isDiffusiveMaterial,
                             const IsBlockedFn &isBlockedMaterial) {
    if (mode_ == AnnealMode::Explicit) {
      const auto invDx2 = NumericType(1) / (dx * dx);
#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(field.size()); ++i) {
        const auto mat = static_cast<int>(materials[i]);
        const auto centerValue = field[i];
        if (!isDiffusiveMaterial(mat) || isBlockedMaterial(mat)) {
          next[i] = centerValue;
          continue;
        }

        NumericType laplacian = NumericType(0);
        int neighborCount = 0;
        const auto &neighbors = cellSet_->getNeighbors(i);
        for (const auto n : neighbors) {
          if (n < 0)
            continue;
          const auto neighborMat = static_cast<int>(materials[n]);
          if (isBlockedMaterial(neighborMat) || !isDiffusiveMaterial(neighborMat))
            continue;

          laplacian += (field[n] - centerValue);
          neighborCount++;
        }

        if (neighborCount == 0) {
          next[i] = centerValue;
          continue;
        }

        const auto Dlocal = std::max(diffusivityAt(i), NumericType(0));
        auto updated = centerValue + dt * Dlocal * invDx2 * laplacian;
        if (clampNonNegative_ && updated < NumericType(0))
          updated = NumericType(0);
        next[i] = updated;
      }
      return;
    }

    // Implicit backward-Euler solve:
    // (I - dt*D*L) u_new = u_old, solved by Gauss-Seidel on the sparse
    // neighborhood stencil.
    next = field;
    const auto invDx2 = NumericType(1) / (dx * dx);
    for (int iter = 0; iter < implicitMaxIterations_; ++iter) {
      NumericType maxDelta = NumericType(0);
      NumericType maxValue = NumericType(0);
      for (int i = 0; i < static_cast<int>(field.size()); ++i) {
        const auto mat = static_cast<int>(materials[i]);
        if (!isDiffusiveMaterial(mat) || isBlockedMaterial(mat)) {
          next[i] = field[i];
          maxValue = std::max(maxValue, std::abs(next[i]));
          continue;
        }

        int neighborCount = 0;
        NumericType neighborSum = NumericType(0);
        const auto &neighbors = cellSet_->getNeighbors(i);
        for (const auto n : neighbors) {
          if (n < 0)
            continue;
          const auto neighborMat = static_cast<int>(materials[n]);
          if (isBlockedMaterial(neighborMat) || !isDiffusiveMaterial(neighborMat))
            continue;
          neighborSum += next[n];
          neighborCount++;
        }

        if (neighborCount == 0) {
          next[i] = field[i];
          maxValue = std::max(maxValue, std::abs(next[i]));
          continue;
        }

        const auto Dlocal = std::max(diffusivityAt(i), NumericType(0));
        const auto alpha = dt * Dlocal * invDx2;
        const auto diag = NumericType(1) + alpha * neighborCount;
        const auto updated = (field[i] + alpha * neighborSum) / diag;
        const auto old = next[i];
        next[i] = clampNonNegative_ ? std::max(updated, NumericType(0)) : updated;
        maxDelta = std::max(maxDelta, std::abs(next[i] - old));
        maxValue = std::max(maxValue, std::abs(next[i]));
      }

      const auto rel = maxDelta / std::max(maxValue, NumericType(1e-30));
      if (rel < implicitRelativeTolerance_)
        break;
    }
  }

  SmartPointer<DenseCellSet<NumericType, D>> cellSet_;
  std::string speciesLabel_ = "concentration";
  std::string materialLabel_ = "Material";
  NumericType duration_ = NumericType(0);
  NumericType timeStep_ = NumericType(-1);
  NumericType stabilityFactor_ = NumericType(0.45);
  bool clampNonNegative_ = true;
  AnnealMode mode_ = AnnealMode::Explicit;
  int implicitMaxIterations_ = 400;
  NumericType implicitRelativeTolerance_ = NumericType(1e-6);

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
  std::string damageLastImpLabel_ = "Damage_LastImp";
  std::string interstitialLabel_ = "Interstitial";
  std::string vacancyLabel_ = "Vacancy";
  NumericType defectFromDamageHistoryWeight_ = NumericType(0);
  NumericType defectFromDamageLastImpWeight_ = NumericType(1);
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
  NumericType tedCoefficient_ = NumericType(0);
  NumericType tedNormalization_ = NumericType(1e20);
  bool defectClusteringEnabled_ = false;
  std::string clusterLabel_ = "ICluster";
  NumericType ikfi_ = NumericType(0);
  NumericType ikfc_ = NumericType(0);
  NumericType ikr_ = NumericType(0);
  NumericType clusterInitFraction_ = NumericType(0);
};

} // namespace viennacs
