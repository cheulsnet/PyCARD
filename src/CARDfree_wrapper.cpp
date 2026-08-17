// CARDfree_wrapper.cpp
#include "CARDfree.hpp"

namespace py = pybind11;

// Binding code
PYBIND11_MODULE(CARDfree_module, m) {
    
    m.def("CARDfree", &CARDfree, "CARDfree function for Python");
}