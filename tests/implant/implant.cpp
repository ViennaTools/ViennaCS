#include <csAnneal.hpp>
#include <csDenseCellSet.hpp>
#include <csImplant.hpp>
#include <csImplantModel.hpp>

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMaterialMap.hpp>

#include <vcTestAsserts.hpp>

#include <cmath>

namespace ls = viennals;
namespace cs = viennacs;

using T = double;
constexpr int D = 2;

// Box-profile implant model: depth profile = 1 for depth in [0, maxDepth),
// lateral = 1 everywhere.
struct BoxImplantModel : public cs::ImplantModel<T, D> {
  T maxDepth_;
  explicit BoxImplantModel(T maxDepth) : maxDepth_(maxDepth) {}
  T getDepthProfile(T depth) override {
    return (depth >= 0 && depth < maxDepth_) ? T(1) : T(0);
  }
  T getLateralProfile(T /*offset*/, T /*depth*/) override { return T(1); }
  T getMaxDepth() override { return maxDepth_; }
  T getMaxLateralRange() override { return T(0); }
};

// 2D slab: substrate (material 1) from y=0 to y=subH, cover (material 0) above.
// Material 0 is the default voidMaterial in csImplant, so the beam passes
// through cover cells and enters the substrate.
cs::SmartPointer<cs::DenseCellSet<T, D>>
makeImplantCellSet(T gridDelta, T xExtent, T subH, T topSpace) {
  T bounds[2 * D] = {-0.5 * xExtent, 0.5 * xExtent, -gridDelta,
                     subH + topSpace};
  ls::BoundaryConditionEnum bc[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T origin[D] = {};
  T normal[D] = {};
  normal[D - 1] = 1.;

  auto makePlane = [&](T y) {
    origin[D - 1] = y;
    auto levelSet =
        ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bc, gridDelta);
    ls::MakeGeometry<T, D>(
        levelSet, ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
        .apply();
    return levelSet;
  };

  auto bottom = makePlane(0.);
  auto top = makePlane(subH);
  ls::BooleanOperation<T, D>(top, bottom, ls::BooleanOperationEnum::UNION)
      .apply();

  auto matMap = ls::SmartPointer<ls::MaterialMap>::New();
  matMap->insertNextMaterial(1);
  matMap->insertNextMaterial(1);

  std::vector<ls::SmartPointer<ls::Domain<T, D>>> levelSets = {bottom, top};

  auto cellSet = cs::SmartPointer<cs::DenseCellSet<T, D>>::New();
  cellSet->setCellSetPosition(true);
  cellSet->setCoverMaterial(
      0); // material 0 = void; beam passes through cover cells
  cellSet->fromLevelSets(levelSets, matMap, topSpace);
  return cellSet;
}

// --- Test 1: Implant deposits non-zero concentration into substrate cells ---
void testImplantDepositsConcentration() {
  const T gridDelta = 1.0;
  const T subH = 10.0;

  auto cellSet = makeImplantCellSet(gridDelta, 6.0, subH, 2.0);
  auto model = cs::SmartPointer<BoxImplantModel>::New(5.0);

  cs::Implant<T, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.apply();

  auto conc = cellSet->getScalarData("concentration");
  VC_TEST_ASSERT(conc != nullptr);

  auto mats = cellSet->getScalarData("Material");
  T total = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    VC_TEST_ASSERT((*conc)[i] >= 0.);
    if (static_cast<int>((*mats)[i]) == 1)
      total += (*conc)[i];
  }
  VC_TEST_ASSERT(total > 0.);
}

// --- Test 2: Mask material stops the beam; no concentration is deposited ---
void testMaskBlocksImplant() {
  const T gridDelta = 1.0;
  const T subH = 10.0;

  auto cellSet = makeImplantCellSet(gridDelta, 6.0, subH, 2.0);

  // Convert all substrate cells to mask material (3)
  auto mats = cellSet->getScalarData("Material");
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    if (static_cast<int>((*mats)[i]) == 1)
      (*mats)[i] = T(3);

  auto model = cs::SmartPointer<BoxImplantModel>::New(5.0);

  cs::Implant<T, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.setMaskMaterials(3);
  implant.apply();

  // Concentration field is always created; all values must be zero.
  auto conc = cellSet->getScalarData("concentration");
  VC_TEST_ASSERT(conc != nullptr);
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    VC_TEST_ASSERT((*conc)[i] == 0.);
}

// --- Test 3: Damage model creates and populates the "Damage" field ---
void testDamageFieldCreated() {
  const T gridDelta = 1.0;
  const T subH = 10.0;

  auto cellSet = makeImplantCellSet(gridDelta, 6.0, subH, 2.0);
  auto model = cs::SmartPointer<BoxImplantModel>::New(5.0);
  auto dmgModel = cs::SmartPointer<BoxImplantModel>::New(5.0);

  cs::Implant<T, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.setDamageModel(dmgModel);
  implant.apply();

  auto dmg = cellSet->getScalarData("Damage");
  VC_TEST_ASSERT(dmg != nullptr);

  auto mats = cellSet->getScalarData("Material");
  T total = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    if (static_cast<int>((*mats)[i]) == 1)
      total += (*dmg)[i];
  VC_TEST_ASSERT(total > 0.);
}

// --- Test 4: WaferDose control scales total concentration linearly with dose
// ---
void testDoseControlScalesConcentration() {
  const T gridDelta = 1.0;
  const T subH = 10.0;
  const T dose1 = 1e13;
  const T dose2 = 2e13;
  const T lengthUnit = 1e-7; // 1 nm in cm

  auto model = cs::SmartPointer<BoxImplantModel>::New(5.0);

  auto run = [&](T dose) -> T {
    auto cellSet = makeImplantCellSet(gridDelta, 6.0, subH, 2.0);
    cs::Implant<T, D> implant;
    implant.setCellSet(cellSet);
    implant.setImplantModel(model);
    implant.setDoseControl(cs::ImplantDoseControl::WaferDose);
    implant.setDose(dose);
    implant.setLengthUnitInCm(lengthUnit);
    implant.apply();
    auto conc = cellSet->getScalarData("concentration");
    T total = 0.;
    if (conc)
      for (auto v : *conc)
        total += v;
    return total;
  };

  T total1 = run(dose1);
  T total2 = run(dose2);

  VC_TEST_ASSERT(total1 > 0.);
  VC_TEST_ASSERT(total2 > 0.);
  // Doubling dose doubles concentration: dose appears as a scalar multiplier.
  VC_TEST_ASSERT_ISCLOSE(total2 / total1, 2.0, 1e-10);
}

// --- Test 5: Anneal spreads the implanted profile and conserves total dose ---
// Runs a full implant → anneal pipeline.  After annealing the peak
// concentration must decrease (diffusion spreads the box profile) while the
// total integrated dose is approximately conserved under reflective boundary
// conditions.
void testAnnealSpreadsImplantedProfile() {
  const T gridDelta = 1.0;
  const T subH = 20.0;

  auto cellSet = makeImplantCellSet(gridDelta, 6.0, subH, 2.0);
  auto model =
      cs::SmartPointer<BoxImplantModel>::New(4.0); // deposit in top 4 cells

  cs::Implant<T, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.apply();

  auto mats = cellSet->getScalarData("Material");
  auto conc = cellSet->getScalarData("concentration");
  VC_TEST_ASSERT(conc != nullptr);

  // Record pre-anneal peak and total in substrate
  T peakBefore = 0., sumBefore = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1)
      continue;
    peakBefore = std::max(peakBefore, (*conc)[i]);
    sumBefore += (*conc)[i];
  }
  VC_TEST_ASSERT(peakBefore > 0.);

  cs::Anneal<T, D> anneal;
  anneal.setCellSet(cellSet);
  anneal.setSpeciesLabel("concentration");
  anneal.setDiffusionMaterials({1});
  anneal.setBlockingMaterials(
      {0}); // cover material blocks diffusion out of substrate
  anneal.setDiffusionCoefficient(1.0);
  anneal.setMode(cs::AnnealMode::Explicit);
  anneal.setDuration(10.0);
  anneal.apply();

  T peakAfter = 0., sumAfter = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1)
      continue;
    VC_TEST_ASSERT((*conc)[i] >= 0.);
    peakAfter = std::max(peakAfter, (*conc)[i]);
    sumAfter += (*conc)[i];
  }

  // Diffusion must have spread the profile: peak decreases
  VC_TEST_ASSERT(peakAfter < peakBefore);
  // Total dose is conserved under reflective BCs (within 0.1 %)
  VC_TEST_ASSERT_ISCLOSE(sumBefore, sumAfter, 1e-3 * sumBefore);
}

int main() {
  testImplantDepositsConcentration();
  testMaskBlocksImplant();
  testDamageFieldCreated();
  testDoseControlScalesConcentration();
  testAnnealSpreadsImplantedProfile();
}
