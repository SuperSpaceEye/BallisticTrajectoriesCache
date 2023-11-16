#include <pybind11/pybind11.h>
#include "py_bindings.h"

namespace py = pybind11;
namespace by = BallisticFunctions;

PYBIND11_MODULE(CannonBallisticFunctions, m) {
    m.def("time_in_air", &by::time_in_air);
    m.def("calculate_pitch", &py_try_pitch);

    m.def("make_dataset", &by::make_dataset);
}