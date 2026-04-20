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
  module.def("getDefaultImplantTablePath", []() {
    py::module_ os = py::module_::import("os");
    py::module_ pkg = py::module_::import("viennacs");
    py::object moduleFile = pkg.attr("__file__");
    return py::str(os.attr("path").attr("join")(
        os.attr("path").attr("dirname")(moduleFile), "data", "implant",
        "implant_moments.csv"));
  });

  py::class_<constants::PearsonIVParameters<T>>(module, "PearsonIVParameters")
      .def(py::init<>())
      .def_readwrite("mu", &constants::PearsonIVParameters<T>::mu)
      .def_readwrite("sigma", &constants::PearsonIVParameters<T>::sigma)
      .def_readwrite("beta", &constants::PearsonIVParameters<T>::beta)
      .def_readwrite("gamma", &constants::PearsonIVParameters<T>::gamma);

  py::enum_<LateralStraggleModel>(module, "LateralStraggleModel")
      .value("Constant", LateralStraggleModel::Constant)
      .value("SentaurusS3", LateralStraggleModel::SentaurusS3)
      .value("Taurus", LateralStraggleModel::Taurus)
      .value("Dios", LateralStraggleModel::Dios);

  py::enum_<ImplantDoseControl>(module, "ImplantDoseControl")
      .value("Off", ImplantDoseControl::Off)
      .value("WaferDose", ImplantDoseControl::WaferDose)
      .value("BeamDose", ImplantDoseControl::BeamDose);

  py::class_<LateralStraggleParameters<T>>(module, "LateralStraggleParameters")
      .def(py::init<>())
      .def_readwrite("model", &LateralStraggleParameters<T>::model)
      .def_readwrite("mu", &LateralStraggleParameters<T>::mu)
      .def_readwrite("sigma", &LateralStraggleParameters<T>::sigma)
      .def_readwrite("scale", &LateralStraggleParameters<T>::scale)
      .def_readwrite("lv", &LateralStraggleParameters<T>::lv)
      .def_readwrite("deltaSigma",
                     &LateralStraggleParameters<T>::deltaSigma)
      .def_readwrite("referenceRange",
                     &LateralStraggleParameters<T>::referenceRange)
      .def_readwrite("p1", &LateralStraggleParameters<T>::p1)
      .def_readwrite("p2", &LateralStraggleParameters<T>::p2)
      .def_readwrite("p3", &LateralStraggleParameters<T>::p3)
      .def_readwrite("p4", &LateralStraggleParameters<T>::p4)
      .def_readwrite("p5", &LateralStraggleParameters<T>::p5);

  py::class_<tables::ImplantTableEntry<T>>(module, "ImplantTableEntry")
      .def(py::init<>())
      .def_readwrite("species", &tables::ImplantTableEntry<T>::species)
      .def_readwrite("material", &tables::ImplantTableEntry<T>::material)
      .def_readwrite("substrateType",
                     &tables::ImplantTableEntry<T>::substrateType)
      .def_readwrite("modelType", &tables::ImplantTableEntry<T>::modelType)
      .def_readwrite("energyKeV", &tables::ImplantTableEntry<T>::energyKeV)
      .def_readwrite("tiltDeg", &tables::ImplantTableEntry<T>::tiltDeg)
      .def_readwrite("rotationDeg", &tables::ImplantTableEntry<T>::rotationDeg)
      .def_readwrite("screenThickness",
                     &tables::ImplantTableEntry<T>::screenThickness)
      .def_readwrite("headFraction",
                     &tables::ImplantTableEntry<T>::headFraction)
      .def_readwrite("screenDecayLength",
                     &tables::ImplantTableEntry<T>::screenDecayLength)
      .def_readwrite("damageDecay", &tables::ImplantTableEntry<T>::damageDecay)
      .def_readwrite("tiltDecayDeg", &tables::ImplantTableEntry<T>::tiltDecayDeg)
      .def_readwrite("reshapeStrength",
                     &tables::ImplantTableEntry<T>::reshapeStrength)
      .def_readwrite("headParams", &tables::ImplantTableEntry<T>::headParams)
      .def_readwrite("headLateralMu",
                     &tables::ImplantTableEntry<T>::headLateralMu)
      .def_readwrite("headLateralSigma",
                     &tables::ImplantTableEntry<T>::headLateralSigma)
      .def_readwrite("headLateralModel",
                     &tables::ImplantTableEntry<T>::headLateralModel)
      .def_readwrite("headLateralScale",
                     &tables::ImplantTableEntry<T>::headLateralScale)
      .def_readwrite("headLateralLv",
                     &tables::ImplantTableEntry<T>::headLateralLv)
      .def_readwrite("headLateralDeltaSigma",
                     &tables::ImplantTableEntry<T>::headLateralDeltaSigma)
      .def_readwrite("headLateralP1",
                     &tables::ImplantTableEntry<T>::headLateralP1)
      .def_readwrite("headLateralP2",
                     &tables::ImplantTableEntry<T>::headLateralP2)
      .def_readwrite("headLateralP3",
                     &tables::ImplantTableEntry<T>::headLateralP3)
      .def_readwrite("headLateralP4",
                     &tables::ImplantTableEntry<T>::headLateralP4)
      .def_readwrite("headLateralP5",
                     &tables::ImplantTableEntry<T>::headLateralP5)
      .def_readwrite("tailParams", &tables::ImplantTableEntry<T>::tailParams)
      .def_readwrite("tailLateralMu",
                     &tables::ImplantTableEntry<T>::tailLateralMu)
      .def_readwrite("tailLateralSigma",
                     &tables::ImplantTableEntry<T>::tailLateralSigma)
      .def_readwrite("tailLateralModel",
                     &tables::ImplantTableEntry<T>::tailLateralModel)
      .def_readwrite("tailLateralScale",
                     &tables::ImplantTableEntry<T>::tailLateralScale)
      .def_readwrite("tailLateralLv",
                     &tables::ImplantTableEntry<T>::tailLateralLv)
      .def_readwrite("tailLateralDeltaSigma",
                     &tables::ImplantTableEntry<T>::tailLateralDeltaSigma)
      .def_readwrite("tailLateralP1",
                     &tables::ImplantTableEntry<T>::tailLateralP1)
      .def_readwrite("tailLateralP2",
                     &tables::ImplantTableEntry<T>::tailLateralP2)
      .def_readwrite("tailLateralP3",
                     &tables::ImplantTableEntry<T>::tailLateralP3)
      .def_readwrite("tailLateralP4",
                     &tables::ImplantTableEntry<T>::tailLateralP4)
      .def_readwrite("tailLateralP5",
                     &tables::ImplantTableEntry<T>::tailLateralP5);

  py::class_<tables::ImplantRecipe<T>>(module, "ImplantRecipe")
      .def(py::init<>())
      .def_readwrite("species", &tables::ImplantRecipe<T>::species)
      .def_readwrite("material", &tables::ImplantRecipe<T>::material)
      .def_readwrite("substrateType", &tables::ImplantRecipe<T>::substrateType)
      .def_readwrite("preferredModel", &tables::ImplantRecipe<T>::preferredModel)
      .def_readwrite("tableFileName", &tables::ImplantRecipe<T>::tableFileName)
      .def_readwrite("useTableLookup", &tables::ImplantRecipe<T>::useTableLookup)
      .def_readwrite("energyKeV", &tables::ImplantRecipe<T>::energyKeV)
      .def_readwrite("tiltDeg", &tables::ImplantRecipe<T>::tiltDeg)
      .def_readwrite("rotationDeg", &tables::ImplantRecipe<T>::rotationDeg)
      .def_readwrite("screenThickness", &tables::ImplantRecipe<T>::screenThickness)
      .def_readwrite("damageLevel", &tables::ImplantRecipe<T>::damageLevel)
      .def_readwrite("entry", &tables::ImplantRecipe<T>::entry);

  py::class_<tables::ImplantTable<T>>(module, "ImplantTable")
      .def(py::init<>())
      .def(py::init<const std::string &>())
      .def("load", &tables::ImplantTable<T>::load)
      .def("getEntries", &tables::ImplantTable<T>::getEntries)
      .def("lookup", &tables::ImplantTable<T>::lookup, py::arg("species"),
           py::arg("material"), py::arg("substrateType"),
           py::arg("energyKeV"), py::arg("tiltDeg"), py::arg("rotationDeg"),
           py::arg("screenThickness") = T(0),
           py::arg("preferredModel") = "auto");

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
