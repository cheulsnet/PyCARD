// CARDref_wrapper.cpp
#include "CARDref.hpp"

namespace py = pybind11;

// Binding code
PYBIND11_MODULE(CARDref_module, m) {
    
    m.def("CARDref", &CARDref, "CARDref function for Python");
}