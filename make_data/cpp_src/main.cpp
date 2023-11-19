#include <pybind11/pybind11.h>
#include "py_bindings.h"

namespace py = pybind11;
namespace by = BallisticFunctions;

PYBIND11_MODULE(CannonBallisticFunctions, m) {
    m.def("calculate_pitch_iterative", &py_calculate_pitch_iterative);
    m.def("calculate_pitch_Endal", &py_calculate_pitch_endals);

    m.def("make_dataset_iterative", &py_make_dataset_iterative);
    m.def("make_dataset_endals", &py_make_dataset_endals);
}