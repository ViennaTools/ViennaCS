#include <csDenseCellSet.hpp>
#include <csImplant.hpp>

namespace cs = viennacs;

template <typename T, int D> auto makePlane(const T xExtent, const T yExtent,  const T gridDelta) {
  T bounds[2 * D] = {0.};
  if constexpr (D == 2) {
    bounds[0] = 0;
    bounds[1] = xExtent;
    bounds[2] = 0;
    bounds[3] = yExtent;
  } else {
    bounds[0] = 0;
    bounds[1] = xExtent;
    bounds[2] = 0;
    bounds[3] = yExtent;
    bounds[4] = -1.;
    bounds[5] = 1.;
  }

  viennals::BoundaryConditionEnum<D> boundaryConds[D] = {
      viennals::BoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY};
  boundaryConds[D - 1] = viennals::BoundaryConditionEnum<D>::INFINITE_BOUNDARY; // warum genau eine Infinite Boundary condition?

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

// overwrite ImplantModel to describe implementation logic
template <class NumericType, int D>
class GaussianModel : public cs::ImplantModel<NumericType, D> {
public:
  GaussianModel(const NumericType a) : a_(a) {};

  NumericType getDepthProfile(NumericType depth) override {
      NumericType stddev = 1;
      NumericType mean = 5;
      const double pi = 3.14159265358979323846;
      return (1.0 / (stddev * std::sqrt(2 * pi))) *
             std::exp(-0.5 * std::pow((depth - mean) / stddev, 2))*999;
  }

private:
  const NumericType a_ = 1.;
};


int main(int argc, char** argv) {
  constexpr int D = 2;
  using NumericType = double;

    // Parsing the parameters
    cs::util::Parameters params;
    if (argc > 1) {
        params.readConfigFile(argv[1]);
    } else {
        std::cout << "Usage: " << argv[0] << " <config file>" << std::endl;
        return 1;
    }

  auto levelSet = makePlane<NumericType, D>(params.get("xExtent"), params.get("yExtent"), params.get("gridDelta"));

  auto materialMap = cs::SmartPointer<viennals::MaterialMap>::New();
  materialMap->insertNextMaterial(0);
  auto cellSet = cs::SmartPointer<cs::DenseCellSet<NumericType, D>>::New();
  cellSet->setCoverMaterial(0); //do I need this line?
  cellSet->fromLevelSets({levelSet}, materialMap, params.get("depth"));
  auto concentration = cellSet->addScalarData("aaa", 0);
  auto model = cs::SmartPointer<GaussianModel<NumericType, D>>::New(1.);

  cs::Implant<double, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.apply();

  cellSet->writeVTU("final.vtu");

  return 0;
}
