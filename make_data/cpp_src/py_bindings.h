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
        int32_t initial_speed, int32_t length, int amin = -30, int amax = 60, float gravity = 0.05, float drag = 0.99,
        float max_delta_t_error = 1, int32_t max_steps=100000,
        int num_iterations = 5, int num_elements = 20, bool check_impossible = true, bool lambertW = false) {
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

    return BallisticFunctions::calculate_pitch(cannon, target, initial_speed, length, amin, amax, gravity, drag, max_delta_t_error, max_steps, num_iterations, num_elements, check_impossible, lambertW);
}

#endif //BALLISTICCACHE_PY_BINDINGS_H
