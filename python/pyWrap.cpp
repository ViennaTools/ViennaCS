#include "pyWrap.hpp"

PYBIND11_MODULE(VIENNACS_MODULE_NAME, module) {
  module.doc() = "ViennaCS cell-set framework and solver kernels.";

  module.attr("__version__") = versionString();
  module.attr("version") = versionString();

  module.def("setNumThreads", &omp_set_num_threads);

  py::enum_<ImplantDoseControl>(module, "ImplantDoseControl")
      .value("Off", ImplantDoseControl::Off)
      .value("WaferDose", ImplantDoseControl::WaferDose)
      .value("BeamDose", ImplantDoseControl::BeamDose);

  py::enum_<DiffusionSolverMode>(module, "DiffusionSolverMode")
      .value("Explicit", DiffusionSolverMode::Explicit)
      .value("GaussSeidel", DiffusionSolverMode::GaussSeidel)
      .value("Implicit", DiffusionSolverMode::GaussSeidel);

  auto m2 = module.def_submodule("d2", "2D bindings");
  m2.attr("__name__") = "viennacs.d2";
  m2.attr("__package__") = "viennacs";
  bindAPI<2>(m2);

  auto m3 = module.def_submodule("d3", "3D bindings");
  m3.attr("__name__") = "viennacs.d3";
  m3.attr("__package__") = "viennacs";
  bindAPI<3>(m3);
}
