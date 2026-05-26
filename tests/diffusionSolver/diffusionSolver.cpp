#include <csDenseCellSet.hpp>
#include <csDiffusionSolver.hpp>
#include <csBoundaryDiffusionSolver.hpp>

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMaterialMap.hpp>

#include <vcTestAsserts.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace ls = viennals;
namespace cs = viennacs;

using T = double;
constexpr int D = 2;

// Substrate (material 1) from y=0 to subH; cover (material 2) above.
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
  cellSet->buildNeighborhood();
  return cellSet;
}

cs::SmartPointer<cs::DenseCellSet<T, D>>
makeEmbeddedPlaneCellSet(T gridDelta, T interfaceY, T topDepth) {
  T bounds[2 * D] = {-1., 1., -gridDelta, interfaceY + topDepth + gridDelta};
  ls::BoundaryConditionEnum bc[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T origin[D] = {};
  T normal[D] = {};
  origin[D - 1] = interfaceY;
  normal[D - 1] = 1.;

  auto surface = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bc, gridDelta);
  ls::MakeGeometry<T, D>(
      surface, ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
      .apply();

  auto matMap = ls::SmartPointer<ls::MaterialMap>::New();
  matMap->insertNextMaterial(1);

  auto cellSet = cs::SmartPointer<cs::DenseCellSet<T, D>>::New();
  cellSet->setCellSetPosition(true);
  cellSet->setCoverMaterial(2);
  cellSet->enableEmbeddedBoundaries();
  cellSet->fromLevelSets({surface}, matMap, interfaceY + topDepth);
  cellSet->buildNeighborhood();
  return cellSet;
}

// --- Test 1: explicit Euler diffuses a spike in the substrate ---
// Mark topmost substrate cells as "source" (material 3, fixed to 1).
// After N steps the remainder of the substrate (material 1) must have
// received concentration, and source cells must stay at 1.
void testExplicitWithInternalSource() {
  const T gridDelta = 0.5;
  const T subH = 4.0;

  auto cellSet = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
  cellSet->addScalarData("field", 0.);
  auto field = cellSet->getScalarData("field");
  auto mats  = cellSet->getScalarData("Material");

  // Tag the topmost substrate layer as source (material 3) and set it to 1.
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1) continue;
    const T y = cellSet->getCellCenter(i)[D - 1];
    if (std::fabs(y - (subH - 0.5 * gridDelta)) < 0.5 * gridDelta) {
      (*mats)[i]  = 3.;
      (*field)[i] = 1.;
    }
  }

  const T D_coeff  = 1.0;
  const T dt       = 0.4 * gridDelta * gridDelta / (2. * D * D_coeff);

  cs::DiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::Explicit);
  solver.setDiffusiveMaterials({1});
  solver.setSourceMaterials({3}); // Dirichlet: fixed at their initial value

  for (int s = 0; s < 20; ++s)
    solver.step(*field, *mats, gridDelta, dt, D_coeff);

  // Source cells must hold Dirichlet value
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) == 3)
      VC_TEST_ASSERT_ISCLOSE((*field)[i], 1.0, 1e-12);
  }

  // Concentration must have diffused into material-1 cells
  T sumConc = 0.;
  int nSub  = 0;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) == 1) {
      sumConc += (*field)[i];
      ++nSub;
    }
  }
  VC_TEST_ASSERT(nSub > 0);
  VC_TEST_ASSERT(sumConc > 0.);

  // No cell should go negative
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    VC_TEST_ASSERT((*field)[i] >= 0.);
}

// --- Test 2: GaussSeidel reaches the same steady state as explicit ---
void testGaussSeidelSteadyState() {
  const T gridDelta = 0.5;
  const T subH = 2.0;
  const T D_coeff = 1.0;

  // Run to steady state with a given mode and return mean substrate concentration.
  auto runToSteadyState = [&](cs::DiffusionSolverMode mode) -> T {
    auto cellSet = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
    cellSet->addScalarData("field", 0.);
    auto field = cellSet->getScalarData("field");
    auto mats  = cellSet->getScalarData("Material");

    // Spike at mid-depth to drive diffusion
    for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
      if (static_cast<int>((*mats)[i]) == 1 &&
          std::fabs(cellSet->getCellCenter(i)[D - 1] - 0.5 * subH) < gridDelta)
        (*field)[i] = 1.;
    }

    cs::DiffusionSolver<T, D> solver;
    solver.setCellSet(cellSet);
    solver.setMode(mode);
    solver.setDiffusiveMaterials({1});

    const T dt = (mode == cs::DiffusionSolverMode::Explicit)
                     ? 0.4 * gridDelta * gridDelta / (2. * D * D_coeff)
                     : 5.0;
    const T totalTime = 20.0;
    T time = 0.;
    while (time < totalTime) {
      T step = std::min(dt, totalTime - time);
      solver.step(*field, *mats, gridDelta, step, D_coeff);
      time += step;
    }

    T sum = 0.;
    int n  = 0;
    for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
      if (static_cast<int>((*mats)[i]) == 1) { sum += (*field)[i]; ++n; }
    }
    return n > 0 ? sum / n : 0.;
  };

  const T meanExplicit = runToSteadyState(cs::DiffusionSolverMode::Explicit);
  const T meanGS       = runToSteadyState(cs::DiffusionSolverMode::GaussSeidel);

  VC_TEST_ASSERT(meanExplicit > 0.);
  VC_TEST_ASSERT(meanGS > 0.);
  // Both modes conserve mass; means must agree within 5 %
  VC_TEST_ASSERT_ISCLOSE(meanExplicit, meanGS, 0.05 * meanExplicit);
}

// --- Test 3: blocked cells are never updated ---
void testBlockedMaterial() {
  const T gridDelta = 0.5;
  const T subH = 2.0;
  const T D_coeff = 1.0;

  auto cellSet = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);

  // Add field FIRST so all subsequent getScalarData calls return stable pointers.
  cellSet->addScalarData("field", 0.);
  auto mats  = cellSet->getScalarData("Material");
  auto field = cellSet->getScalarData("field");

  // Mark the bottommost substrate layer as blocked (material 3)
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) == 1 &&
        cellSet->getCellCenter(i)[D - 1] < gridDelta)
      (*mats)[i] = 3.;
  }

  // Uniform concentration in the diffusive region
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    if (static_cast<int>((*mats)[i]) == 1) (*field)[i] = 1.;

  cs::DiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::Explicit);
  solver.setDiffusiveMaterials({1});
  solver.setBlockedMaterials({3});

  const T dt = 0.4 * gridDelta * gridDelta / (2. * D * D_coeff);
  for (int s = 0; s < 50; ++s)
    solver.step(*field, *mats, gridDelta, dt, D_coeff);

  // Blocked cells must remain at zero (they were initialised to 0 and never updated)
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) == 3)
      VC_TEST_ASSERT_ISCLOSE((*field)[i], 0., 1e-12);
  }
}

// --- Test 4: mass conservation under zero-flux (Neumann) BCs ---
// A Gaussian spike in a closed domain must preserve total concentration
// (sum over all diffusive cells) across all three solver modes.
void testMassConservation(cs::DiffusionSolverMode mode) {
  const T gridDelta = 0.5;
  const T subH = 8.0;
  const T D_coeff = 1.0;

  auto cellSet = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
  cellSet->addScalarData("field", 0.);
  auto field = cellSet->getScalarData("field");
  auto mats  = cellSet->getScalarData("Material");

  // Gaussian spike centred in the substrate
  const T centre = 0.5 * subH;
  const T sigma  = 0.5;
  T totalInit = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1) continue;
    const T y = cellSet->getCellCenter(i)[D - 1];
    (*field)[i] = std::exp(-0.5 * (y - centre) * (y - centre) / (sigma * sigma));
    totalInit += (*field)[i];
  }
  VC_TEST_ASSERT(totalInit > 0.);

  cs::DiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(mode);
  solver.setDiffusiveMaterials({1});

  const T dt = (mode == cs::DiffusionSolverMode::Explicit)
                   ? 0.4 * gridDelta * gridDelta / (2. * D * D_coeff)
                   : 2.0;
  const int nSteps = 40;
  for (int s = 0; s < nSteps; ++s)
    solver.step(*field, *mats, gridDelta, dt, D_coeff);

  T totalFinal = 0.;
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    if (static_cast<int>((*mats)[i]) == 1) totalFinal += (*field)[i];

  // Mass must be conserved to within 0.1 %
  VC_TEST_ASSERT_ISCLOSE(totalFinal, totalInit, 1e-3 * totalInit);
}

void testEmbeddedBoundaryDirichletUsesSubGridDistance() {
  const T gridDelta = 0.25;
  const T interfaceY = 0.10;
  const T diffusivity = 1.0;
  const T dt = 1.e-5;

  auto cellSet = makeEmbeddedPlaneCellSet(gridDelta, interfaceY, 0.75);
  auto mats = cellSet->getScalarData("Material");
  std::vector<T> field(cellSet->getNumberOfCells(), 0.);

  int targetCell = -1;
  T targetDistance = 0.;
  for (int cellId = 0; cellId < static_cast<int>(cellSet->getNumberOfCells());
       ++cellId) {
    const auto center = cellSet->getCellCenter(cellId);
    for (const auto pointId : cellSet->getEmbeddedBoundaryPointIds(cellId)) {
      const auto &point = cellSet->getEmbeddedBoundaryPoints()[pointId];
      if (point.levelSetIndex != 0 || point.axis != D - 1)
        continue;
      if (point.coordinate[D - 1] >= center[D - 1])
        continue;
      targetCell = cellId;
      targetDistance = center[D - 1] - point.coordinate[D - 1];
      break;
    }
    if (targetCell >= 0)
      break;
  }
  VC_TEST_ASSERT(targetCell >= 0)
  VC_TEST_ASSERT(targetDistance > 0.)
  VC_TEST_ASSERT(targetDistance < 0.5 * gridDelta)

  const int diffusiveMaterial = static_cast<int>((*mats)[targetCell]);
  cs::BoundaryDiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::Explicit);
  solver.setDiffusiveMaterials({diffusiveMaterial});
  solver.setBoundaryCondition(0,
                              cs::BoundaryCondition<T>::dirichlet(1.));

  solver.step(field, *mats, gridDelta, dt, diffusivity);

  const T expected = dt * 2. * diffusivity /
                     (targetDistance * (targetDistance + gridDelta));
  VC_TEST_ASSERT_ISCLOSE(field[targetCell], expected, 1.e-10)
}

void testEmbeddedBoundaryRobinUsesSubGridDistance() {
  const T gridDelta = 0.25;
  const T interfaceY = 0.10;
  const T diffusivity = 1.0;
  const T transfer = 2.0;
  const T exterior = 1.0;
  const T dt = 1.e-5;

  auto cellSet = makeEmbeddedPlaneCellSet(gridDelta, interfaceY, 0.75);
  auto mats = cellSet->getScalarData("Material");
  std::vector<T> field(cellSet->getNumberOfCells(), 0.);

  int targetCell = -1;
  T targetDistance = 0.;
  for (int cellId = 0; cellId < static_cast<int>(cellSet->getNumberOfCells());
       ++cellId) {
    const auto center = cellSet->getCellCenter(cellId);
    for (const auto pointId : cellSet->getEmbeddedBoundaryPointIds(cellId)) {
      const auto &point = cellSet->getEmbeddedBoundaryPoints()[pointId];
      if (point.levelSetIndex != 0 || point.axis != D - 1)
        continue;
      if (point.coordinate[D - 1] >= center[D - 1])
        continue;
      targetCell = cellId;
      targetDistance = center[D - 1] - point.coordinate[D - 1];
      break;
    }
    if (targetCell >= 0)
      break;
  }
  VC_TEST_ASSERT(targetCell >= 0)

  const int diffusiveMaterial = static_cast<int>((*mats)[targetCell]);
  cs::BoundaryDiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::Explicit);
  solver.setDiffusiveMaterials({diffusiveMaterial});
  solver.setBoundaryCondition(
      0, cs::BoundaryCondition<T>::robin(transfer, exterior));

  solver.step(field, *mats, gridDelta, dt, diffusivity);

  const T conductance = diffusivity / targetDistance;
  const T boundaryConstant = transfer * exterior / (conductance + transfer);
  const T expected = dt * 2. * diffusivity * boundaryConstant /
                     (targetDistance * (targetDistance + gridDelta));
  VC_TEST_ASSERT_ISCLOSE(field[targetCell], expected, 1.e-10)
}

// --- Test 7: explicit source field raises total mass ------------------------
// Verify that a uniform volumetric source S causes the total concentration
// in diffusive cells to increase by approximately nCells * dt * S after one
// step. Diffusion redistributes mass but does not create or destroy it; only
// the source changes the total. We compare with-source vs without-source runs
// starting from the same initial state.
void testSourceFieldExplicit() {
  const T gridDelta = 0.5;
  const T subH = 4.0;
  const T D_coeff = 1.0;
  const T C0 = 0.;
  const T S = 2.;
  const T dt = 0.01;

  auto sumMat1 = [](const std::vector<T> *f, const std::vector<T> *m,
                    int n) {
    T s = 0.;
    for (int i = 0; i < n; ++i)
      if (static_cast<int>((*m)[i]) == 1) s += (*f)[i];
    return s;
  };

  // Run without source
  auto csNoSrc = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
  csNoSrc->addScalarData("field", C0);
  auto fNoSrc = csNoSrc->getScalarData("field");
  auto mNoSrc = csNoSrc->getScalarData("Material");
  {
    cs::DiffusionSolver<T, D> solver;
    solver.setCellSet(csNoSrc);
    solver.setMode(cs::DiffusionSolverMode::Explicit);
    solver.setDiffusiveMaterials({1});
    solver.step(*fNoSrc, *mNoSrc, gridDelta, dt, D_coeff);
  }

  // Run with source
  auto csSrc = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
  csSrc->addScalarData("field", C0);
  auto fSrc = csSrc->getScalarData("field");
  auto mSrc = csSrc->getScalarData("Material");
  std::vector<T> sourceField(csSrc->getNumberOfCells(), S);
  {
    cs::DiffusionSolver<T, D> solver;
    solver.setCellSet(csSrc);
    solver.setMode(cs::DiffusionSolverMode::Explicit);
    solver.setDiffusiveMaterials({1});
    solver.setSourceField(&sourceField);
    solver.step(*fSrc, *mSrc, gridDelta, dt, D_coeff);
  }

  const int n = csSrc->getNumberOfCells();
  const T totalNoSrc = sumMat1(fNoSrc, mNoSrc, n);
  const T totalSrc   = sumMat1(fSrc,   mSrc,   n);

  // Source must have added positive mass
  VC_TEST_ASSERT(totalSrc > totalNoSrc);
  // The increase must be at least 50% of dt*S*nDiffusiveCells (accounting for
  // the fact that some active cells skip the source when they have no
  // contributing neighbours — a known conserved quantity is dt*S*nActive).
  int nDiffuse = 0;
  for (int i = 0; i < n; ++i)
    if (static_cast<int>((*mSrc)[i]) == 1) ++nDiffuse;
  VC_TEST_ASSERT(nDiffuse > 0);
  const T expectedIncrease = dt * S * nDiffuse;
  VC_TEST_ASSERT_ISCLOSE(totalSrc - totalNoSrc, expectedIncrease, 0.1 * expectedIncrease);
}

// --- Test 8: GaussSeidel source field raises concentration -----------------
void testSourceFieldGaussSeidel() {
  const T gridDelta = 0.5;
  const T subH = 2.0;
  const T D_coeff = 1.0;
  const T C0 = 1.;
  const T S = 2.;
  const T dt = 0.01;

  auto cellSet = makeSlabCellSet(gridDelta, 4.0, subH, 0.5);
  cellSet->addScalarData("field", 0.);
  auto field = cellSet->getScalarData("field");
  auto mats  = cellSet->getScalarData("Material");

  for (int i = 0; i < cellSet->getNumberOfCells(); ++i)
    if (static_cast<int>((*mats)[i]) == 1) (*field)[i] = C0;

  std::vector<T> sourceField(cellSet->getNumberOfCells(), S);

  cs::DiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::GaussSeidel);
  solver.setDiffusiveMaterials({1});
  solver.setSourceField(&sourceField);

  solver.step(*field, *mats, gridDelta, dt, D_coeff);

  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1) continue;
    VC_TEST_ASSERT((*field)[i] > C0);
  }
}

// --- Test 9: embedded-boundary source field --------------------------------
// Same uniform-field argument: with no boundary condition, the embedded
// boundary acts as a zero-flux wall, so the Laplacian is zero and after one
// explicit step c = C0 + dt*S exactly.
void testEmbeddedBoundarySourceField() {
  const T gridDelta = 0.25;
  const T interfaceY = 0.10;
  const T D_coeff = 1.0;
  const T C0 = 1.;
  const T S = 2.;
  const T dt = 1e-5;

  auto cellSet = makeEmbeddedPlaneCellSet(gridDelta, interfaceY, 0.75);
  auto mats = cellSet->getScalarData("Material");

  std::vector<T> field(cellSet->getNumberOfCells(), C0);
  std::vector<T> sourceField(cellSet->getNumberOfCells(), S);

  cs::BoundaryDiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::Explicit);
  solver.setDiffusiveMaterials({1});
  solver.setSourceField(&sourceField);

  solver.step(field, *mats, gridDelta, dt, D_coeff);

  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*mats)[i]) != 1) continue;
    VC_TEST_ASSERT_ISCLOSE(field[i], C0 + dt * S, 1e-10);
  }
}

// --- Test 10: getEffectiveMinBoundaryDistance < gridDelta/2 for sub-grid BC --
void testGetEffectiveMinBoundaryDistance() {
  const T gridDelta = 0.25;
  const T interfaceY = 0.10;

  auto cellSet = makeEmbeddedPlaneCellSet(gridDelta, interfaceY, 0.75);

  cs::BoundaryDiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);

  const T effectiveDx = solver.getEffectiveMinBoundaryDistance(gridDelta);
  VC_TEST_ASSERT(effectiveDx > 0.)
  VC_TEST_ASSERT(effectiveDx < 0.5 * gridDelta)
}

int main() {
  testExplicitWithInternalSource();
  testGaussSeidelSteadyState();
  testBlockedMaterial();
  testMassConservation(cs::DiffusionSolverMode::Explicit);
  testMassConservation(cs::DiffusionSolverMode::GaussSeidel);
  testEmbeddedBoundaryDirichletUsesSubGridDistance();
  testEmbeddedBoundaryRobinUsesSubGridDistance();
  testSourceFieldExplicit();
  testSourceFieldGaussSeidel();
  testEmbeddedBoundarySourceField();
  testGetEffectiveMinBoundaryDistance();
}
