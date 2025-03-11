#include <csDenseCellSet.hpp>
#include <csImplant.hpp>
#include <models/csImplantGaussian.hpp>
#include <models/csImplantPearson.hpp>

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <vcUtil.hpp>

namespace cs = viennacs;

template <typename T, int D>
auto makePlane(const T xExtent, const T yExtent, const T gridDelta) {
  T bounds[2 * D] = {0.};
  bounds[1] = xExtent;
  bounds[3] = yExtent;

  viennals::BoundaryConditionEnum boundaryConds[D] = {
      viennals::BoundaryConditionEnum::REFLECTIVE_BOUNDARY};
  boundaryConds[D - 1] = viennals::BoundaryConditionEnum::INFINITE_BOUNDARY;

  auto levelSet = viennals::SmartPointer<viennals::Domain<T, D>>::New(
      bounds, boundaryConds, gridDelta);

  T origin[D] = {0.};
  T normal[D] = {0.};
  normal[D - 1] = 1.;
  viennals::MakeGeometry<T, D>(
      levelSet,
      viennals::SmartPointer<viennals::Plane<T, D>>::New(origin, normal))
      .apply();
  return levelSet;
}

int main(int argc, char **argv) {
  constexpr int D = 3;
  using NumericType = double;

  // Parsing the parameters
  cs::util::Parameters params;
  if (argc > 1) {
    params.readConfigFile(argv[1]);
  } else {
    std::cout << "Usage: " << argv[0] << " <config file>" << std::endl;
    return 1;
  }

  auto ls = makePlane<NumericType, D>(10., 10., 0.5);
  std::vector<cs::SmartPointer<viennals::Domain<NumericType, D>>> levelSets;
  levelSets.push_back(ls);

  auto cellSet = cs::SmartPointer<cs::DenseCellSet<NumericType, D>>::New(
      levelSets, nullptr, -10, false);

  cellSet->writeVTU("test.vtu");

  auto model = cs::SmartPointer<viennacs::ImplantGaussian<NumericType, D>>::New(
      5., 1., 2., 1.);
  //  model = cs::SmartPointer<PearsonIVModel<NumericType, D>>::New(1.);
  // auto model = cs::SmartPointer<PearsonIVModel<NumericType, D>>::New(1.);
  cellSet->addScalarData("concentration", 0);
  cs::Implant<double, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.setImplantAngle(params.get("angle"));
  implant.apply();

  cellSet->writeVTU("final.vtu");
  return 0;
}
