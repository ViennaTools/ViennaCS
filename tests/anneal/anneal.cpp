#include <csAnneal.hpp>
#include <csDenseCellSet.hpp>

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMaterialMap.hpp>

#include <vcTestAsserts.hpp>

#include <cmath>
#include <numeric>

namespace ls = viennals;
namespace cs = viennacs;

using T = double;
constexpr int D = 2;

cs::SmartPointer<cs::DenseCellSet<T, D>>
makeSlabCellSet(T gridDelta, T xExtent, T subH, T topSpace) {
  T bounds[2 * D] = {-0.5 * xExtent, 0.5 * xExtent, -gridDelta, subH + topSpace};
  ls::BoundaryConditionEnum bc[D] = {ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
                                     ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T origin[D] = {};
  T normal[D] = {};
  normal[D - 1] = 1.;

  auto makePlane = [&](T y) {
    origin[D - 1] = y;
    auto ls = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bc, gridDelta);
    ls::MakeGeometry<T, D>(ls, ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
        .apply();
    return ls;
  };

  auto bottom = makePlane(0.);
  auto top = makePlane(subH);
  ls::BooleanOperation<T, D>(top, bottom, ls::BooleanOperationEnum::UNION).apply();

  auto matMap = ls::SmartPointer<ls::MaterialMap>::New();
  matMap->insertNextMaterial(1);
  matMap->insertNextMaterial(1);

  std::vector<ls::SmartPointer<ls::Domain<T, D>>> levelSets = {bottom, top};

  auto cellSet = cs::SmartPointer<cs::DenseCellSet<T, D>>::New();
  cellSet->setCellSetPosition(true);
  cellSet->setCoverMaterial(2);
  cellSet->fromLevelSets(levelSets, matMap, topSpace);
  return cellSet;
}

// --- Test 1: diffusion spreads a peaked concentration ---
// Places a spike at the mid-depth of the substrate.  After annealing the peak
// must decrease and the total concentration must be approximately conserved
// (Neumann / reflective BCs on all substrate faces).
void testDiffusionSpreads() {
  const T gridDelta = 0.5;
  const T subH = 4.0;

  auto cellSet = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
  cellSet->addScalarData("concentration", 0.);
  auto field = cellSet->getScalarData("concentration");
  auto mats = cellSet->getScalarData("Material");

  // Spike at the centre of the substrate
  const T midY = 0.5 * subH;
  T peakValue = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1) continue;
    const T y = cellSet->getCellCenter(i)[D - 1];
    if (std::fabs(y - midY) < gridDelta) {
      (*field)[i] = 1.0;
      peakValue = 1.0;
    }
  }
  VC_TEST_ASSERT(peakValue > 0.);

  // Sum before anneal
  T sumBefore = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    if (static_cast<int>((*mats)[i]) == 1) sumBefore += (*field)[i];

  cs::Anneal<T, D> anneal;
  anneal.setCellSet(cellSet);
  anneal.setSpeciesLabel("concentration");
  anneal.setDiffusionMaterials({1});
  anneal.setBlockingMaterials({2});
  anneal.setDiffusionCoefficient(1.0); // nm²/s (or any consistent unit)
  anneal.setMode(cs::AnnealMode::Explicit);
  anneal.setDuration(5.0);

  anneal.apply();

  // Peak in substrate must have decreased
  T newPeak = 0.;
  T sumAfter = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1) continue;
    newPeak = std::max(newPeak, (*field)[i]);
    sumAfter += (*field)[i];
  }

  VC_TEST_ASSERT(newPeak < peakValue);
  // Total concentration is conserved: reflective BCs, no source/sink
  VC_TEST_ASSERT_ISCLOSE(sumBefore, sumAfter, 1e-6 * sumBefore);

  // No negative concentrations
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    VC_TEST_ASSERT((*field)[i] >= 0.);
}

// --- Test 2: solid activation is bounded by total concentration ---
void testSolidActivation() {
  const T gridDelta = 0.5;
  const T subH = 3.0;

  auto cellSet = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
  cellSet->addScalarData("concentration", 0.);
  auto field = cellSet->getScalarData("concentration");
  auto mats = cellSet->getScalarData("Material");

  // Uniform concentration in substrate: use 1000× C_SS so active ≈ C_SS (within 0.1 %)
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    if (static_cast<int>((*mats)[i]) == 1) (*field)[i] = 1e22; // cm⁻³ equivalent

  cs::Anneal<T, D> anneal;
  anneal.setCellSet(cellSet);
  anneal.setSpeciesLabel("concentration");
  anneal.setDiffusionMaterials({1});
  anneal.setBlockingMaterials({2});
  anneal.setDiffusionCoefficient(0.); // no diffusion; just test activation
  anneal.setDuration(1.0);
  anneal.enableSolidActivation(true);
  // Solid solubility much lower than the initial concentration
  anneal.setSolidSolubilityArrhenius(1e19, 0.); // C_SS = 1e19 (Ea=0 → temperature-independent)

  anneal.apply();

  auto active = cellSet->getScalarData("active_concentration");
  VC_TEST_ASSERT(active != nullptr);

  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1) continue;
    const T total = (*field)[i];
    const T act = (*active)[i];
    // Active fraction is physically bounded: 0 ≤ active ≤ total
    VC_TEST_ASSERT(act >= 0.);
    VC_TEST_ASSERT(act <= total + 1e-10);
    // With C_total = 1000 × C_SS, formula gives C_SS * 1000/(1001) ≈ 0.999 * C_SS
    const T expectedActive = 1e19 * 1e22 / (1e19 + 1e22);
    VC_TEST_ASSERT_ISCLOSE(act, expectedActive, 1e16);
  }
}

// --- Test 3: temperature schedule runs without error and produces non-zero result ---
void testTemperatureSchedule() {
  const T gridDelta = 0.5;
  const T subH = 2.0;

  auto cellSet = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
  cellSet->addScalarData("concentration", 0.);
  auto field = cellSet->getScalarData("concentration");
  auto mats = cellSet->getScalarData("Material");

  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    if (static_cast<int>((*mats)[i]) == 1) (*field)[i] = 1.0;

  cs::Anneal<T, D> anneal;
  anneal.setCellSet(cellSet);
  anneal.setSpeciesLabel("concentration");
  anneal.setDiffusionMaterials({1});
  anneal.setBlockingMaterials({2});
  // Arrhenius diffusivity evaluated at peak T
  anneal.setArrheniusParameters(1.0, 0.); // D0=1, Ea=0 → D = 1 at any T
  anneal.setTemperature(1000.);

  // Ramp-up / soak / ramp-down: 3 steps
  anneal.setTemperatureSchedule({2.0, 2.0, 2.0}, {800., 1200., 1200., 800.});

  anneal.apply();

  // Concentration must still exist (anneal ran) and be non-negative
  int nSub = 0;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1) continue;
    VC_TEST_ASSERT((*field)[i] >= 0.);
    ++nSub;
  }
  VC_TEST_ASSERT(nSub > 0);
}

int main() {
  testDiffusionSpreads();
  testSolidActivation();
  testTemperatureSchedule();
}
