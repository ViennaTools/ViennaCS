#include "pyWrap.hpp"

PYBIND11_MODULE(VIENNACS_MODULE_NAME, module) {
  module.doc() = "ViennaCS is a header-only C++ cell set library, which adds "
                 "the possibility of using volumetric representations on top "
                 "of existing level-set functionalities for surfaces. Combined "
                 "with ray tracing techniques, this enables the simulation of "
                 "particle scattering and ion implantation.";

  // set version string of python module
  module.attr("__version__") = versionString();
  module.attr("version") = versionString();

  // wrap omp_set_num_threads to control number of threads
  module.def("setNumThreads", &omp_set_num_threads);

  // Submodule for 2D
  auto m2 = module.def_submodule("d2", "2D bindings");
  m2.attr("__name__") = "viennacs.d2";
  m2.attr("__package__") = "viennacs";
  bindAPI<2>(m2);

  // Submodule for 3D
  auto m3 = module.def_submodule("d3", "3D bindings");
  m3.attr("__name__") = "viennacs.d3";
  m3.attr("__package__") = "viennacs";
  bindAPI<3>(m3);
}
