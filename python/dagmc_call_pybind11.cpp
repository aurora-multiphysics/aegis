#include <pybind11/pybind11.h>
#include "dagmc_call_lib/DataClass.hpp"
#include <string>

namespace py = pybind11;
//using your namespace for your project

// here dagmc_call is the python module name, must equal to modname in Cmakelists.txt
PYBIND11_MODULE(dagmc_call_wrap, m) {
    py::class_<DataClass>(m, "DataClass")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("message", &DataClass::message)
        .def("setMessage", &DataClass::setMessage)
        .def("__repr__", [](const DataClass& r)
            {
                return "DataClass(" + r.message() + ")";
            }
        );
}
