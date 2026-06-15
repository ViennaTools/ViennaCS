#include "pyWrap.hpp"

PYBIND11_MODULE(VIENNACS_MODULE_NAME, module) {
  module.doc() = "ViennaCS cell-set framework and solver kernels.";

  module.attr("__version__") = versionString();
  module.attr("version") = versionString();

  module.def("setNumThreads", &omp_set_num_threads);

  // BoundaryConditionType and BoundaryCondition are dimension-independent;
  // register them here rather than inside bindAPI<D> to avoid pybind11's
  // global-registry "already registered" error on the second bindAPI call.
  py::enum_<BoundaryConditionType>(module, "BoundaryConditionType")
      .value("Neumann", BoundaryConditionType::Neumann)
      .value("Dirichlet", BoundaryConditionType::Dirichlet)
      .value("Robin", BoundaryConditionType::Robin)
      .export_values();

  py::class_<BoundaryCondition<T>>(module, "BoundaryCondition")
      .def(py::init<>())
      .def_readwrite("type", &BoundaryCondition<T>::type)
      .def_readwrite("value", &BoundaryCondition<T>::value)
      .def_readwrite("transferCoefficient",
                     &BoundaryCondition<T>::transferCoefficient)
      .def_static("dirichlet", &BoundaryCondition<T>::dirichlet,
                  py::arg("boundaryValue"),
                  "Create a Dirichlet condition with the given surface value.")
      .def_static("neumann", &BoundaryCondition<T>::neumann,
                  py::arg("outwardFlux") = T(0),
                  "Create a Neumann condition with the given outward flux "
                  "(default 0 = zero flux).")
      .def_static("robin", &BoundaryCondition<T>::robin,
                  py::arg("transfer"), py::arg("exteriorValue"),
                  "Create a Robin condition: flux = transfer*(c - exteriorValue).");

  py::enum_<ImplantDoseControl>(module, "ImplantDoseControl")
      .value("Off", ImplantDoseControl::Off)
      .value("WaferDose", ImplantDoseControl::WaferDose)
      .value("BeamDose", ImplantDoseControl::BeamDose);

  py::enum_<DiffusionSolverMode>(module, "DiffusionSolverMode")
      .value("Explicit", DiffusionSolverMode::Explicit)
      .value("GaussSeidel", DiffusionSolverMode::GaussSeidel);

  auto m2 = module.def_submodule("d2", "2D bindings");
  m2.attr("__name__") = "viennacs.d2";
  m2.attr("__package__") = "viennacs";
  bindAPI<2>(m2);

  auto m3 = module.def_submodule("d3", "3D bindings");
  m3.attr("__name__") = "viennacs.d3";
  m3.attr("__package__") = "viennacs";
  bindAPI<3>(m3);
}
