//
// Created by spaceeye on 17.08.23.
//

#ifndef BALLISTICCACHE_PY_BINDINGS_H
#define BALLISTICCACHE_PY_BINDINGS_H

#include "ballistic_functions.h"
#include <pybind11/stl.h>
#include <iostream>

inline std::pair<std::array<ftype, 3>, std::array<ftype, 3>>
py_calculate_pitch_iterative(
        pybind11::tuple & cannon_t,
        pybind11::tuple & target_t,
        ftype initial_speed, ftype length, int amin = -30, int amax = 60, ftype gravity = 0.05, ftype drag = 0.99,
        ftype max_delta_t_error = 1, int32_t max_steps= 100000,
        int num_iterations = 5, int num_elements = 20, bool check_impossible = true, bool lambertW = false) {
    const std::array<ftype, 3> cannon{
            cannon_t[0].cast<ftype>(),
            cannon_t[1].cast<ftype>(),
            cannon_t[2].cast<ftype>(),
    };
    const std::array<ftype, 3> target{
            target_t[0].cast<ftype>(),
            target_t[1].cast<ftype>(),
            target_t[2].cast<ftype>(),
    };

    return BallisticFunctions::calculate_pitch(cannon, target, initial_speed, length, amin, amax, gravity, drag, max_delta_t_error, max_steps, num_iterations, num_elements, check_impossible, lambertW);
}

inline std::pair<std::array<ftype, 3>, std::array<ftype, 3>>
py_calculate_pitch_endals(
        pybind11::tuple & cannon_t,
        pybind11::tuple & target_t,
        ftype initial_speed, ftype length, ftype amin, ftype amax, ftype gravity, ftype drag,
        ftype mult_coeff, ftype acceptable_error_range,
        ftype starting_from_t, ftype max_mult_depth,
        ftype starting_multiplier, ftype starting_max_depth
        ) {
    const std::array<ftype, 3> cannon{
            cannon_t[0].cast<ftype>(),
            cannon_t[1].cast<ftype>(),
            cannon_t[2].cast<ftype>(),
    };
    const std::array<ftype, 3> target{
            target_t[0].cast<ftype>(),
            target_t[1].cast<ftype>(),
            target_t[2].cast<ftype>(),
    };

    return BallisticFunctions::calculate_pitch_Endal(cannon, target, initial_speed, length, gravity, drag, amin, amax, mult_coeff, acceptable_error_range, starting_from_t, max_mult_depth, starting_multiplier, starting_max_depth);
}

std::vector<thread_dataset_type> py_make_dataset_iterative(
        ftype initial_speed,
        ftype length,
        ftype amin,
        ftype amax,
        ftype gravity,
        ftype drag,
        ftype max_delta_t_error,
        int max_steps,
        int num_iterations,
        int num_elements,
        bool check_impossible,
        bool lambertW,

        int max_height_above,
        int max_height_below,
        int num_threads,
        bool verbose,
        int max_distance,
        ftype step,
        int impossible_cutoff
        ) {

    auto fn = [initial_speed, length, amin, amax, gravity, drag, max_delta_t_error, max_steps, num_iterations, num_elements, check_impossible, lambertW]
            (const std::array<ftype, 3> & cannon, const std::array<ftype, 3> & target){
        return BallisticFunctions::calculate_pitch(
                cannon, target, initial_speed, length, amin, amax, gravity, drag, max_delta_t_error, max_steps, num_iterations,
                num_elements, check_impossible, lambertW);
    };

    return BallisticFunctions::make_dataset(fn, length, max_height_above, max_height_below, num_threads, verbose, max_distance, step, impossible_cutoff);
}

std::vector<thread_dataset_type> py_make_dataset_endals(
        ftype initial_speed,
        ftype length,
        ftype amin,
        ftype amax,
        ftype gravity,
        ftype drag,

        ftype mult_coeff, ftype acceptable_error_range,
        ftype starting_from_t, ftype max_mult_depth,
        ftype starting_multiplier, ftype starting_max_depth,

        int max_height_above,
        int max_height_below,
        int num_threads,
        bool verbose,
        int max_distance,
        ftype step,
        int impossible_cutoff
) {

    auto fn = [initial_speed, length, amin, amax, gravity, drag, mult_coeff, acceptable_error_range, starting_from_t, max_mult_depth, starting_multiplier, starting_max_depth]
            (const std::array<ftype, 3> & cannon, const std::array<ftype, 3> & target){
        return BallisticFunctions::calculate_pitch_Endal(cannon, target, initial_speed, length, gravity, drag, amin, amax, mult_coeff, acceptable_error_range, starting_from_t, max_mult_depth, starting_multiplier, starting_max_depth);
    };

    return BallisticFunctions::make_dataset(fn, length, max_height_above, max_height_below, num_threads, verbose, max_distance, step, impossible_cutoff);
}

#endif //BALLISTICCACHE_PY_BINDINGS_H
