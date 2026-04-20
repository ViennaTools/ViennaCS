#include <csDenseCellSet.hpp>
#include <csTracing.hpp>

#include "geometry.hpp"
#include "particle.hpp"

namespace cs = viennacs;

int main(int argc, char **argv) {
  constexpr int D = 3;
  using NumericType = double;

  // Parse the parameters
  cs::util::Parameters params;
  if (argc > 1) {
    params.readConfigFile(argv[1]);
  } else {
    // Try default config file
    params.readConfigFile("config.txt");
    if (params.m.empty()) {
      std::cout << "No configuration file provided!" << std::endl;
      std::cout << "Usage: " << argv[0] << " <config file>" << std::endl;
      return 1;
    }
  }

  auto levelSet = makeFin<NumericType, D>(params);

  // Create a cellSet
  auto cellSet = cs::SmartPointer<cs::DenseCellSet<NumericType, D>>::New();
  cellSet->fromLevelSets(levelSet, nullptr, -5.);
  cellSet->writeVTU("initial.vtu");

  auto particle = std::make_unique<PlasmaDamageIon<NumericType, D>>(
      params.get("ionEnergy"), params.get("meanFreePath"));

  cs::Tracing<NumericType, D> tracer;
  tracer.setCellSet(cellSet);
  tracer.setParticle(particle);
  tracer.setNumberOfRaysPerPoint(params.get("raysPerPoint"));
  tracer.apply();

  cellSet->writeVTU("final.vtu");
}
