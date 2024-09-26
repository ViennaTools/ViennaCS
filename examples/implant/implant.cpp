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
  NumericType getDepthProfile(NumericType depth, cs::util::Parameters params) override {
      NumericType mu = params.get("mu");
      NumericType sigma = params.get("sigma");
      const double pi = 3.14159265358979323846;
      return (1.0 / (sigma * std::sqrt(2 * pi))) *
             std::exp(-0.5 * std::pow((depth - mu) / sigma, 2))*999;
  }

private:
  const NumericType a_ = 1.;
};


template <class NumericType, int D>
class PearsonIVModel : public cs::ImplantModel<NumericType, D> {
public:
    PearsonIVModel(const NumericType a) : a_(a) {};
    NumericType getDepthProfile(NumericType depth, cs::util::Parameters params) override {
        NumericType mu = params.get("mu");
        NumericType sigma = params.get("sigma");
        NumericType gamma = params.get("gamma");
        NumericType beta = params.get("beta");
        depth  = depth-mu;
        NumericType A = 10*beta - 12*gamma*gamma - 18;
        NumericType a = -gamma*sigma*(beta + 3) / A;
        NumericType b_0 = -sigma*sigma*(4*beta-3*gamma*gamma) / A;
        NumericType b_1 = a;
        NumericType b_2 = -(2*beta - 3*gamma*gamma - 6) / A;
        NumericType m = 1/(2*b_2);
        return pow(abs(b_0+b_1*depth+b_2*depth * depth), m) * exp((-(b_1 / b_2 + 2*a) /
            sqrt(4*b_0*b_2 - b_1*b_1)) * atan((2*b_2*depth + b_1) / sqrt(4*b_0*b_2 - b_1*b_1)));
    }

private:
    const NumericType a_ = 1;
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
  auto concentration = cellSet->addScalarData("concentration", 0);

  auto model = cs::SmartPointer<GaussianModel<NumericType, D>>::New(1.);

  cs::Implant<double, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.setParameters(params);
  implant.apply();

  cellSet->writeVTU("final.vtu");

  return 0;
}
