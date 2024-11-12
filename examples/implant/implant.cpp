#include <csDenseCellSet.hpp>
#include <csImplant.hpp>
#include <psDomain.hpp>
#include <psProcess.hpp>
#include <geometries/psMakeTrench.hpp>
#include <models/psDirectionalEtching.hpp>
#include <models/psIsotropicProcess.hpp>

namespace cs = viennacs;
namespace ps = viennaps;

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



template <class T> class VelocityField : public viennals::VelocityField<T> {
public:
    VelocityField() {}

    /// Should return a scalar value for the velocity at coordinate
    /// for a point of material with the given normalVector.
    T getScalarVelocity(const viennals::Vec3D<T> & /*coordinate*/,
                        int /*material*/,
                        const viennals::Vec3D<T> & /*normalVector*/,
                        unsigned long /*pointId*/) override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dis(0, 1.);

        return dis(gen);
    }
};

template <class T> class MySurfaceModel : public ps::SurfaceModel<T> {
public:
    MySurfaceModel() {}

    ps::SmartPointer<std::vector<T>>
    calculateVelocities(ps::SmartPointer<viennals::PointData<T>> rates,
                        const std::vector<ps::Vec3D<T>> &coordinates,
                        const std::vector<T> &materialIds) override {

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dis(0, 1.);

        std::vector<T> velocities;

        for (unsigned i = 0; i < coordinates.size(); ++i) {
            if (ps::MaterialMap::isMaterial(materialIds[i], ps::Material::Si)) {
                velocities.push_back(dis(gen));
            } else {
                velocities.push_back(0);
            }
        }

        return ps::SmartPointer<std::vector<T>>::New(std::move(velocities));
    }
};

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
             std::exp(-0.5 * std::pow((depth - mu) / sigma, 2));
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

template <class NumericType, int D>
class DualPearsonModel : public cs::ImplantModel<NumericType, D> {
public:
    DualPearsonModel(const NumericType a) : a_(a) {};
    NumericType getDepthProfile(NumericType depth, cs::util::Parameters params) override {
        NumericType mu1 = params.get("mu1");
        NumericType mu2 = params.get("mu2");
        NumericType sigma1 = params.get("sigma1");
        NumericType sigma2 = params.get("sigma2");
        NumericType gamma1 = params.get("gamma1");
        NumericType gamma2 = params.get("gamma2");
        NumericType beta1 = params.get("beta1");
        NumericType beta2 = params.get("beta2");
        NumericType r = params.get("r");
        return (1-r)* singlePearsonIV(depth, mu1, sigma1, gamma1, beta1) + r * singlePearsonIV(depth, mu2, sigma2, gamma2, beta2);
    }

    NumericType singlePearsonIV(NumericType depth, NumericType mu, NumericType sigma, NumericType gamma, NumericType beta){
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

  // creating trench/ mask cover
  auto domain = ps::SmartPointer<ps::Domain<double, 2>>::New();
  ps::MakeTrench<double,2>(domain, params.get("gridDelta"), 10, 3, 0.5, 0.1, 0., 2., false, true,
                            ps::Material::Si).apply();
  ps::Process<double, 2> process;
  process.setDomain(domain);
  auto etchingModel = ps::SmartPointer<ps::DirectionalEtching<double, 2>>::New(
          ps::Vec3D<double>{0., -1., 0.}, 1., -0.1, ps::Material::Mask);
  process.setProcessModel(etchingModel);
  process.setProcessDuration(0.1);
  process.apply();

  domain->saveSurfaceMesh("initial.vtp");
  domain->generateCellSet(0.1, ps::Material::Si, true);
  auto cellSet = domain->getCellSet();
  auto model = cs::SmartPointer<GaussianModel<NumericType, D>>::New(1.);
  cellSet->addScalarData("concentration", 0);
  cs::Implant<double, D> implant;
  implant.setCellSet(cellSet);
  implant.setImplantModel(model);
  implant.setParameters(params);
  implant.apply(params);

  cellSet->writeVTU("final.vtu");
  return 0;
}
