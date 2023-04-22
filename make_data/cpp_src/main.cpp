#include <pybind11/pybind11.h>
#include "ballistic_functions.h"

namespace py = pybind11;

PYBIND11_MODULE(CannonBallisticFunctions, m) {
    m.def("time_in_air", &time_in_air);
    m.def("rough_pitch_estimation", &py_rough_pitch_estimation);
    m.def("fine_pitch_estimation", &py_fine_pitch_estimation);
    m.def("calculate_pitch", &py_try_pitch);

    m.def("make_dataset", &make_dataset<calculate_pitch>);
}