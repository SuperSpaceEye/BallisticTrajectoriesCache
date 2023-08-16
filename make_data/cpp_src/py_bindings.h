//
// Created by spaceeye on 17.08.23.
//

#ifndef BALLISTICCACHE_PY_BINDINGS_H
#define BALLISTICCACHE_PY_BINDINGS_H

#include "ballistic_functions.h"
#include <pybind11/stl.h>

inline std::pair<std::array<float, 3>, std::array<float, 3>>
py_try_pitch(
        pybind11::tuple & cannon_t,
        pybind11::tuple & target_t,
        int32_t power, int32_t length, int32_t max_steps) {
    std::array<float, 3> cannon{
            cannon_t[0].cast<float>(),
            cannon_t[1].cast<float>(),
            cannon_t[2].cast<float>(),
    };
    std::array<float, 3> target{
            target_t[0].cast<float>(),
            target_t[1].cast<float>(),
            target_t[2].cast<float>(),
    };

    return BallisticFunctions::calculate_pitch(cannon, target, power, length, max_steps);
}

#endif //BALLISTICCACHE_PY_BINDINGS_H
