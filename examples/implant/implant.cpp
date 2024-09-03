#include <csDenseCellSet.hpp>
#include <csImplant.hpp>

namespace cs = viennacs;

template <typename T, int D> auto makePlane(const T extent, const T gridDelta) {
  T bounds[2 * D] = {0.};
  if constexpr (D == 2) {
    bounds[0] = -extent / 2.;
    bounds[1] = extent / 2.;
    bounds[2] = -1;
    bounds[3] = 1;
  } else {
    bounds[0] = -extent / 2.;
    bounds[1] = extent / 2.;
    bounds[2] = -extent / 2.;
    bounds[3] = extent / 2.;
    bounds[4] = -1;
    bounds[5] = 1.;
  }

  viennals::BoundaryConditionEnum<D> boundaryConds[D] = {
      viennals::BoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY};
  boundaryConds[D - 1] = viennals::BoundaryConditionEnum<D>::INFINITE_BOUNDARY;

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

template <class NumericType, int D>
class MyImplantModel : public cs::ImplantModel<NumericType, D> {
public:
  MyImplantModel(const NumericType a) : a_(a) {};

  NumericType getDepthProfile(NumericType depth) override { return a_; }

private:
  const NumericType a_ = 1.;
};

int main() {

  constexpr int D = 2;
  using NumericType = double;

  auto levelSet = makePlane<NumericType, D>(10., 0.2);

  auto materialMap = cs::SmartPointer<viennals::MaterialMap>::New();
  materialMap->insertNextMaterial(1);
  auto cellSet = cs::SmartPointer<cs::DenseCellSet<NumericType, D>>::New();
  cellSet->setCoverMaterial(1);
  cellSet->fromLevelSets({levelSet}, materialMap, -5.);
  cellSet->writeVTU("initial.vtu");

  auto model = cs::SmartPointer<MyImplantModel<NumericType, D>>::New(1.);

  cs::Implant<double, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.apply();

  return 0;
}