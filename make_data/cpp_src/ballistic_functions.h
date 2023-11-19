//
// Created by spaceeye on 14.04.23.
//

#ifndef CANNONBALLISTICCALCULATOR_BALLISTIC_FUNCTIONS_H
#define CANNONBALLISTICCALCULATOR_BALLISTIC_FUNCTIONS_H

#include <cstdint>
#include <vector>
#include <cmath>
#include <array>
#include <thread>
#include <future>
#include <chrono>

#include <iostream>
#include <ranges>

#include "linspace.h"
#include "LambertW.hpp"

#define DESTRUCTURE3(v1, v2, v3, fn) auto t = fn; v1 = t[0]; v2 = t[1]; v3 = t[2];

using ftype = float;
using item_type = std::tuple<int32_t, int32_t, ftype, ftype, ftype, ftype, ftype, ftype>;
using thread_dataset_type = std::vector<item_type>;
using pitch_fn_return_type = std::pair<std::array<ftype, 3>, std::array<ftype, 3>>;
#define PITCH_FN_ARGS std::array<ftype, 3>, std::array<ftype, 3>

namespace BallisticFunctions {

inline auto range(int start, int end, int step=1) {
    auto step_fun = [=](auto x) { return x * step + start; };
    auto end_pred = [=](auto x) { return step>0 ? x <= end: x > end; };
    auto range =
            std::views::iota(0)
            | std::views::transform(step_fun)
            | std::views::take_while(end_pred);
    return range;
}

inline std::vector<ftype> flinspace(ftype start, ftype stop, int num_elements, ftype min, ftype max) {
    return linspace<ftype>(std::max(start, min), std::min(stop, max), num_elements);
}

std::array<ftype, 3> get_root(const std::vector<std::array<ftype, 3>> & data,
                              bool from_end) {
    if (from_end) {
        for (int i = data.size()-2; i >= 0; i--) {
            if (data[i][0] > data[i+1][0]) {return data[i+1];}
        }
        return data[0];
    } else {
        for (int i = 1; i < data.size(); i++) {
            if (data[i-1][0] < data[i][0]) {return data[i-1];}
        }
        return data.back();
    }
}

//calculates time when it goes above and then below target y pos.
inline std::pair<double, double> time_in_air(ftype y0, ftype y, ftype Vy,
                                             ftype gravity = 0.05, ftype drag=0.99, int32_t max_steps = 100000) {
    int64_t t = 0;
    int64_t t_below = INT64_MAX;

    if (y0 < y) {
        ftype y0_p;
        while (t < max_steps) {
            y0_p = y0;
            y0 += Vy;
            Vy = drag * Vy - gravity;

            t+=1;

            if (y0>y) {
                t_below = t - 1;
                break;
            }

            //if projectile stopped ascending and didn't go above targets y pos.
            if (y0 - y0_p < 0) {return {-1, -1};}
        }
    }

    while (t < max_steps) {
        y0 += Vy;
        Vy = drag * Vy - gravity;

        t += 1;

        if (y0 <= y) {return {t_below, t};}
    }

    return {t_below, -1};
}

inline std::pair<double, double> time_in_air_lambertw(ftype y0, ftype y, ftype Vy,
                                                      ftype gravity = 0.05, ftype drag=0.99) {
    auto delta = (20 * Vy + 100) * (std::pow(drag, (20 * Vy + 100 + 0.20 * (y0 - y))));

    auto top0 = -utl::LambertW<0>(delta * std::log(drag));
    auto top1 = -utl::LambertW<-1>(delta * std::log(drag));

    auto n0 = (top0/std::log(drag)) + 20 * Vy + 100 + 0.20 * (y0 - y);
    auto n1 = (top1/std::log(drag)) + 20 * Vy + 100 + 0.20 * (y0 - y);

    return {n0, n1};
}

inline double rad(auto deg) {return deg * (M_PI / 180);}

inline std::pair<std::array<ftype, 3>, bool>
try_pitch(ftype pitch_to_try,
          ftype initial_speed,
          ftype length,
          ftype distance,
          const std::array<ftype, 3> &cannon,
          const std::array<ftype, 3> &target,
          ftype gravity = 0.05,
          ftype drag = 0.99,
          int32_t max_steps = 1000000,
          bool lambertW=false) {
    ftype tp_rad = rad(pitch_to_try);

    auto Vw = std::cos(tp_rad) * initial_speed;
    auto Vy = std::sin(tp_rad) * initial_speed;

    ftype x_coord_2d = length * std::cos(tp_rad);

    double horizontal_time_to_target;
    if (drag < 1.0) {
        if (Vw == 0) { return {{-1, -1, -1}, false}; }
        auto part = 1 - (distance - x_coord_2d) / (100 * Vw);
        if (part <= 0) { return {{-1, -1, -1}, false}; }
        horizontal_time_to_target = std::abs(std::log(part) / (std::log(drag) + 1e-200));
    } else {
        horizontal_time_to_target = distance / Vw;
    }

    ftype y_coord_end_of_barrel = cannon[1] + std::sin(tp_rad) * length;

    int64_t t_below, t_above;
    if (lambertW) {
        auto [_t_below, _t_above] = time_in_air_lambertw(y_coord_end_of_barrel, target[1], Vy, gravity, drag);
        t_below = _t_below; t_above = _t_above;
    } else {
        auto [_t_below, _t_above] = time_in_air(y_coord_end_of_barrel, target[1], Vy, gravity, drag, max_steps);
        t_below = _t_below; t_above = _t_above;
    }

    if (t_below < 0) { return {{-1, -1, -1}, false};}

    //if target is above cannon it may hit it on ascension
    auto delta_t = std::min(
            std::abs(horizontal_time_to_target-t_below),
            std::abs(horizontal_time_to_target-t_above)
            );

    return {{(ftype)delta_t, pitch_to_try, (ftype)(delta_t + horizontal_time_to_target)}, true};
}

template<typename A>
inline std::vector<std::array<ftype, 3>> try_pitches(
        A iter,
        ftype initial_speed,
        ftype length,
        ftype distance,
        const std::array<ftype, 3> &cannon,
        const std::array<ftype, 3> &target,
        ftype gravity = 0.05,
        ftype drag = 0.99,
        int32_t max_steps=1000000,
        bool lambertW = false
        ) {
    std::vector<std::array<ftype, 3>> delta_times{};
    delta_times.reserve(20);//preallocate some
    for (auto item: iter) {
        auto [items, is_successful] = try_pitch(item, initial_speed, length, distance, cannon, target, gravity, drag, max_steps, lambertW);
        if (!is_successful) { continue;}
        delta_times.emplace_back(items);
    }
    return delta_times;
}

inline std::array<ftype, 3> min_array(const std::vector<std::array<ftype, 3>> & vec) {
    return *std::min_element(vec.begin(), vec.end(),
                             [](const auto & a, const auto & b){return a[0] < b[0];});
}

inline std::pair<std::array<ftype, 3>, std::array<ftype, 3>> calculate_pitch(
        const std::array<ftype, 3> &cannon,
        const std::array<ftype, 3> &target,
        ftype initial_speed, ftype length, int amin = -30, int amax = 60, ftype gravity = 0.05, ftype drag = 0.99,
        ftype max_delta_t_error = 1, int32_t max_steps=100000,
        int num_iterations = 5, int num_elements = 20, bool check_impossible = true, bool lambertW = false) {
    auto Dx = cannon[0] - target[0];
    auto Dz = cannon[2] - target[2];
    auto distance = std::sqrt(Dx * Dx + Dz * Dz);

    auto delta_times = try_pitches(range(amax, amin  - 1, -1), initial_speed, length, distance, cannon, target, gravity, drag, max_steps, lambertW);
    if (delta_times.empty()) {return {{-1, -1, -1}, {-1, -1, -1}};}

    auto [dT1, p1, at1] = get_root(delta_times, false);
    auto [dT2, p2, at2] = get_root(delta_times, true);

    bool c1=true;
    bool c2= p1 != p2;
    bool same_res = p1 == p2;

    std::vector<std::array<ftype, 3>> dTs1, dTs2;

    for (int i = 0; i < num_iterations; i++) {
        if (c1) { dTs1 = try_pitches(flinspace(p1-std::pow(10.,-i), p1+std::pow(10.,-i), num_elements, amin, amax), initial_speed, length, distance, cannon, target, gravity, drag, max_steps, lambertW);}
        if (c2) { dTs2 = try_pitches(flinspace(p2-std::pow(10.,-i), p2+std::pow(10.,-i), num_elements, amin, amax), initial_speed, length, distance, cannon, target, gravity, drag, max_steps, lambertW);}

        if (c1 && dTs1.empty()) {c1 = false;}
        if (c2 && dTs2.empty()) {c2 = false;}

        if (!c1 && !c2) {return {{-1, -1, -1}, {-1, -1, -1}};}

        if (c1) { DESTRUCTURE3(dT1, p1, at1, min_array(dTs1))}
        if (c2) { DESTRUCTURE3(dT2, p2, at2, min_array(dTs2))}
    }

    if (same_res) {dT2 = dT1; p2 = p1; at2 = at1;}

    std::array<ftype, 3> r1 = {dT1, p1, at1}, r2 = {dT2, p2, at2};
    if (check_impossible && dT1 > max_delta_t_error) { r1 = {-1, -1, -1};}
    if (check_impossible && dT2 > max_delta_t_error) { r2 = {-1, -1, -1};}

    return std::pair<std::array<ftype, 3>, std::array<ftype, 3>> {r1, r2};
}

inline std::pair<std::array<ftype, 3>, std::array<ftype, 3>> solve_with_Endal(
        ftype x_r, ftype h, ftype v_m, ftype L, ftype g, ftype c_d,

        ftype mult_coeff = 0.25, ftype acceptable_error_range = 0.01,
        ftype starting_from_t = 0, ftype max_mult_depth = 8,
        ftype starting_multiplier = 2, ftype starting_max_depth = 32
        ) {

    const auto u = v_m/20.0;

    const auto A = (g * c_d) / (u * (1 - c_d));
    const auto C = (L/(u*x_r)) * ((g*c_d)/(1-c_d)) + h/x_r;
    const auto B = [g, c_d, x_r](ftype t){ return t * (g * c_d) / (1 - c_d) * (1 / x_r);};

    const auto a_r = [A, B, C](ftype t){
        const auto B_ = B(t);
        const auto part = -A*A + B_*B_ + C*C + 2*B_*C + 1;
        if (part < 0) {return -1e+101;}

        return 2 * std::atan(
                (std::sqrt(part) - 1)
                / (A + B_ + C)
                );
    };

    const auto x_r1 = [a_r, u, c_d, L](ftype t){
        const auto a_rr = a_r(t);
        if (a_rr < -1e+100) {return a_rr;}
        return (u * std::cos(a_rr)) / std::log(c_d) * (std::pow(c_d, t) - 1) + L * std::cos(a_rr);
    };

    const auto get_starting = [starting_from_t, x_r1, starting_multiplier, starting_max_depth](){
        const auto s_t = starting_from_t;
        ftype res = x_r1(s_t);

        ftype p = 1;
        ftype depth = 0;
        while (res < -1e+100) {
            res = x_r1(s_t + p);
            p *= starting_multiplier;

            depth += 1;
            if (depth > starting_max_depth) { throw std::runtime_error("Depth exceeds allowed"); }
        }
        return s_t + p;
    };

    ftype s_t = get_starting();

    ftype s_xr = x_r1(s_t);
    ftype s1_xr = x_r1(s_t+1);

    auto find_solution = [s_xr, max_mult_depth, x_r1, x_r, mult_coeff, acceptable_error_range, a_r](ftype s_t, bool inverse, bool comp_flag){
        ftype p_xr = s_xr;
        ftype c_t = 1 * (inverse ? -1 : 1);
        ftype mult = 2;
        int depth = 0;

        auto comp = (comp_flag ? [](ftype a, ftype b){return a > b;} : [](ftype a, ftype b){return a < b;});
        while (true) {
            if (depth > max_mult_depth) { return std::array<ftype, 3> {-1, -1, -1}; }

            auto n_t = c_t * mult;
            auto n_xr = x_r1(n_t + s_t);

            if (n_xr < x_r && comp(p_xr, n_xr)) { return std::array<ftype, 3> {-1, -1, -1}; }

            if (comp(n_xr, x_r)) {
                mult *= mult_coeff;
                depth++;
                continue;
            }

            if (std::abs(n_xr - x_r) <= acceptable_error_range) {
                if (n_t + s_t < 0) {return std::array<ftype, 3> {-1, -1, -1};}
                return std::array<ftype, 3> {(ftype)(std::abs(n_xr - x_r)), (ftype)a_r(n_t+s_t), (ftype)(n_t+s_t)};
            }

            s_t += n_t;
            p_xr = n_xr;
        }
    };

    std::tuple<ftype, bool, bool, bool> args1, args2;
    if (s_xr < x_r && s1_xr - s_xr >= 0) {args1 = {s_t, false, true,  true }; args2 = {s_t, false, false, true}; }
    if (s_xr < x_r && s1_xr - s_xr <  0) {args1 = {s_t, true,  true,  true }; args2 = {s_t, true,  false, true}; }
    if (s_xr > x_r && s1_xr - s_xr >= 0) {args1 = {s_t, true,  false, false}; args2 = {s_t, false, false, true}; }
    if (s_xr > x_r && s1_xr - s_xr <  0) {args1 = {s_t, false, false, false}; args2 = {s_t, true,  false, true}; }

    auto first_solution = find_solution(get<0>(args1), get<1>(args1), get<2>(args1));

    if (first_solution[0] >= 0 and get<3>(args1)) {
        s_t = first_solution[2] + (acceptable_error_range * 2 * (get<1>(args2) ? -1 : 1));
    }

    auto second_solution = find_solution(s_t, get<1>(args2), get<2>(args2));

    return {first_solution, second_solution};
}

inline std::pair<std::array<ftype, 3>, std::array<ftype, 3>> calculate_pitch_Endal(
        const std::array<ftype, 3> &cannon,
        const std::array<ftype, 3> &target,
        ftype initial_speed, ftype length, ftype gravity, ftype drag,
        ftype amin, ftype amax,
        ftype mult_coeff, ftype acceptable_error_range,
        ftype starting_from_t, ftype max_mult_depth,
        ftype starting_multiplier, ftype starting_max_depth
        ) {
    auto Dx = cannon[0] - target[0];
    auto Dz = cannon[2] - target[2];
    auto x_r = std::sqrt(Dx * Dx + Dz * Dz);

    auto h = target[1] - cannon[1];

    auto [s1, s2] = solve_with_Endal(x_r, h, initial_speed*20, length, gravity, drag, mult_coeff, acceptable_error_range, starting_from_t, max_mult_depth, starting_multiplier, starting_max_depth);

    if (s1[1] > rad(amax) || s1[1] < rad(amin)) {s1 = {-1, -1, -1};}
    if (s2[1] > rad(amax) || s2[1] < rad(amin)) {s2 = {-1, -1, -1};}

    return {s1, s2};
}

inline double
calculate_y_line(const std::function<pitch_fn_return_type(PITCH_FN_ARGS)>& pitch_fn,
        thread_dataset_type *dataset, double starting_x, int *points_simulated, int *y_done,
         int max_length,uint32_t impossible_cutoff, ftype step, int y, bool count_cutoff_at_the_start) {
    bool had_result = false;
    int cutoff_count = 0;
    bool dont_change_starting = false;
    if (starting_x == 0) { starting_x = 1; }
    for (double x = starting_x; x < max_length; x += step) {
        auto [res1, res2] = pitch_fn({0, 0, 0}, {(ftype)x, (ftype)y, 0});

//        if (res2[0] >= 0 && res1[0] < 0) {std::swap(res1, res2);}
//        if (res1[0] >= 0 && res2[0] >= 0) {
//            if (res1[2] > res2[2]) {
//                std::swap(res1, res2);
//            }
//        }

        if (res1[0] >= 0 || res2[0] >= 0) {
            dataset->push_back(
               {(int32_t) x, (int32_t) y,
                (ftype) res1[0], (ftype) res1[1], (ftype) res1[2],
                (ftype) res2[0], (ftype) res2[1], (ftype) res2[2]
                });
            had_result = true;
            dont_change_starting = true;
        } else {
            if (!dont_change_starting) {starting_x = x;}
        }
        (*points_simulated)++;

        if ((had_result || count_cutoff_at_the_start)
        && (res1[0] < 0 && res2[0] < 0)) {cutoff_count++;} else { cutoff_count = 0;}
        if (cutoff_count >= impossible_cutoff) { break;}
        }
    (*y_done)++;
    return starting_x;
}

auto make_dataset_thread(
        const std::function<pitch_fn_return_type(PITCH_FN_ARGS)>& pitch_fn,
        thread_dataset_type * dataset,
        int max_height_above,
        int max_height_below,
        int start_pos,
        int num_threads,
        int * points_simulated,
        int * y_done,
        int max_length,
        int length,
        uint32_t impossible_cutoff,
        ftype step = 1,

        uint8_t * done = nullptr
                ) {
    dataset->reserve(100000);

    if (start_pos == 0) {
        calculate_y_line(pitch_fn, dataset, 0, points_simulated, y_done, max_length, impossible_cutoff, step, 0, true);
    }

    double starting_x = length;
    for (int y = start_pos+1; y < max_height_above; y+=num_threads) {
        starting_x = calculate_y_line(pitch_fn, dataset, starting_x, points_simulated, y_done, max_length, impossible_cutoff, step, y, true);
    }

    starting_x = length;
    //y levels below cannon meanwhile can always be hit at some point, so just simulate until it hits a reachable point
    // and only then start calculating cutoff
    for (int y = start_pos-1; y > -max_height_below; y-=num_threads) {
        starting_x = calculate_y_line(pitch_fn, dataset, starting_x, points_simulated, y_done, max_length, impossible_cutoff, step, y, true);
    }

    *done = true;
}

auto make_dataset(
        const std::function<pitch_fn_return_type(PITCH_FN_ARGS)>& pitch_fn,
        ftype length,
        int max_height_above = 256,
        int max_height_below = 256,
        int num_threads=16,
        bool verbose=true,
        int max_length=600,
        float step=1,
        uint32_t impossible_cutoff = 50
        ) {
    using namespace std::chrono_literals;

    std::vector<thread_dataset_type> threads_result;
    std::vector<std::pair<int, int>> threads_progress;
    std::vector<std::thread> threads;
    std::vector<uint8_t> done;

    threads_progress.resize(num_threads);
    threads_result.resize(num_threads);
    threads.resize(num_threads);
    done.resize(num_threads, false);

    for (int i = 0; i < num_threads; i++) {
        threads[i] = std::thread(
                make_dataset_thread,
                pitch_fn,
                &threads_result[i],
                max_height_above,
                max_height_below, i,
                num_threads,
                &threads_progress[i].first,
                &threads_progress[i].second,
                max_length,
                length,
                impossible_cutoff,
                step,
                &done[i]
        );

        threads[i].detach();
    }

    bool not_ended = true;
    while (not_ended) {
        not_ended = false;
        std::this_thread::sleep_for(1s);

        for (auto & is_done: done) {
            if (!is_done) {
                not_ended += true;
                break;
            }
        }
        if (!verbose) { continue;}

        uint64_t y_levels_done = 0;
        uint64_t points_simulated = 0;

        for (const auto & [x_done, y_done]: threads_progress) {
            points_simulated += x_done;
            y_levels_done += y_done;
        }

        std::cout << "points calculated: " << points_simulated << " | y levels calculated: " << y_levels_done << "\n";
    }

    auto stop = std::chrono::high_resolution_clock::now();

    return threads_result;
}

}



#endif //CANNONBALLISTICCALCULATOR_BALLISTIC_FUNCTIONS_H
