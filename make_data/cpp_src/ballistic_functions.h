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

inline std::vector<float> flinspace(float start, float stop, int num_elements, float min, float max) {
    return linspace<float>(std::max(start, min), std::min(stop, max), num_elements);
}

std::array<float, 3> get_root(const std::vector<std::array<float, 3>> & data,
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
inline std::pair<int64_t, int64_t> time_in_air(float y0, float y, float Vy,
                                               float gravity = 0.05, float drag=0.99, int32_t max_steps = 100000) {
    int64_t t = 0;
    int64_t t_below = INT64_MAX;

    if (y0 < y) {
        float y0_p;
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

inline std::pair<int64_t, int64_t> time_in_air_lambertw(float y0, float y, float Vy,
                               float gravity = 0.05, float drag=0.99) {
    auto delta = (20 * Vy + 100) * (std::pow(drag, (20 * Vy + 100 + 0.20 * (y0 - y))));

    auto top0 = -utl::LambertW<0>(delta * std::log(drag));
    auto top1 = -utl::LambertW<-1>(delta * std::log(drag));

    auto n0 = (top0/std::log(drag)) + 20 * Vy + 100 + 0.20 * (y0 - y);
    auto n1 = (top1/std::log(drag)) + 20 * Vy + 100 + 0.20 * (y0 - y);

    return {n0, n1};
}

inline double rad(auto deg) {return deg * (M_PI / 180);}

inline std::pair<std::array<float, 3>, bool>
try_pitch(float pitch_to_try,
          int32_t initial_speed,
          int32_t length,
          float distance,
          const std::array<float, 3> &cannon,
          const std::array<float, 3> &target,
          float gravity = 0.05,
          float drag = 0.99,
          int32_t max_steps = 1000000,
          bool lambertW=false) {
    float tp_rad = rad(pitch_to_try);

    auto Vw = std::cos(tp_rad) * initial_speed;
    auto Vy = std::sin(tp_rad) * initial_speed;

    float x_coord_2d = length * std::cos(tp_rad);

    double horizontal_time_to_target;
    if (drag < 1.0) {
        if (Vw == 0) { return {{-1, -1, -1}, false}; }
        auto part = 1 - (distance - x_coord_2d) / (100 * Vw);
        if (part <= 0) { return {{-1, -1, -1}, false}; }
        horizontal_time_to_target = std::abs(std::log(part) / (std::log(drag) + 1e-200));
    } else {
        horizontal_time_to_target = distance / Vw;
    }

    float y_coord_end_of_barrel = cannon[1] + std::sin(tp_rad) * length;

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

    return {{(float)delta_t, pitch_to_try, (float)(delta_t + horizontal_time_to_target)}, true};
}

template<typename A>
inline std::vector<std::array<float, 3>> try_pitches(
        A iter,
        int32_t initial_speed,
        int32_t length,
        float distance,
        const std::array<float, 3> &cannon,
        const std::array<float, 3> &target,
        float gravity = 0.05,
        float drag = 0.99,
        int32_t max_steps=1000000,
        bool lambertW = false
        ) {
    std::vector<std::array<float, 3>> delta_times{};
    delta_times.reserve(20);//preallocate some
    for (auto item: iter) {
        auto [items, is_successful] = try_pitch(item, initial_speed, length, distance, cannon, target, gravity, drag, max_steps, lambertW);
        if (!is_successful) { continue;}
        delta_times.emplace_back(items);
    }
    return delta_times;
}

inline std::array<float, 3> min_array(const std::vector<std::array<float, 3>> & vec) {
    return *std::min_element(vec.begin(), vec.end(),
                             [](const auto & a, const auto & b){return a[0] < b[0];});
}

inline std::pair<std::array<float, 3>, std::array<float, 3>> calculate_pitch(
        const std::array<float, 3> &cannon,
        const std::array<float, 3> &target,
        int32_t initial_speed, int32_t length, int amin = -30, int amax = 60, float gravity = 0.05, float drag = 0.99,
        float max_delta_t_error = 1, int32_t max_steps=100000,
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

    std::vector<std::array<float, 3>> dTs1, dTs2;

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

    std::array<float, 3> r1 = {dT1, p1, at1}, r2 = {dT2, p2, at2};
    if (check_impossible && dT1 > max_delta_t_error) { r1 = {-1, -1, -1};}
    if (check_impossible && dT2 > max_delta_t_error) { r2 = {-1, -1, -1};}

    return std::pair<std::array<float, 3>, std::array<float, 3>> {r1, r2};
}

inline double
calculate_y_line(std::vector<std::tuple<int32_t, int32_t, float, float, float>> *dataset, int32_t charges, double starting_x, double length,
                 int *points_simulated, int *y_done, int max_length, int max_simulation_steps,
                 uint32_t impossible_cutoff, float max_delta_t_error, float step, int y,
                 bool count_cutoff_at_the_start, int amin = -30, int amax = 60, float gravity = 0.05, float drag = 0.99,
                 int num_iterations = 5, int num_elements = 20, bool check_impossible = true, bool lambertW = true) {
    bool had_result = false;
    int cutoff_count = 0;
    bool dont_change_starting = false;
    for (double x = starting_x; x < max_length; x += step) {
        auto [res1, res2] = calculate_pitch({0, 0, 0}, {(float)x, (float)y, 0}, charges, length, amin, amax, gravity, drag, max_delta_t_error, max_simulation_steps, num_iterations, num_elements, check_impossible, lambertW);
//        auto res = ((res1[0] < res2[0]) && (res1[0] >= 0)) ? res1 : res2;
        auto res = res2;
        if (res[0] >= 0) {
            dataset->push_back({(int32_t)x, (int32_t)y, (float)res[0], (float)res[1], (float)res[2]});
            had_result = true;
            dont_change_starting = true;
        } else {
            if (!dont_change_starting) {starting_x = x;}
        }
        (*points_simulated)++;

        if ((had_result || count_cutoff_at_the_start) && res[0] < 0) {cutoff_count++;} else { cutoff_count = 0;}
        if (cutoff_count >= impossible_cutoff) { break;}
        }
    (*y_done)++;
    return starting_x;
}

auto make_dataset_thread(
        std::vector<std::tuple<int32_t, int32_t, float, float, float>> * dataset,
                         int32_t charges,
                         int length,
                         int max_height_above,
                         int max_height_below,
                         int start_pos,
                         int num_threads,
                         int * points_simulated,
                         int * y_done,
                         int max_simulation_steps,
                         int max_length,
                         uint32_t impossible_cutoff,
                         float max_delta_t_error,
                         float step = 1,
                         int amin=-30, int amax=60, float gravity = 0.05, float drag = 0.99,
                         int num_iterations=5, int num_elements=20,
                         bool check_impossible=true,
                         bool lambertW = false,

                         uint8_t * done = nullptr
                                 ) {
    dataset->reserve(100000);

    if (start_pos == 0) {
    calculate_y_line(dataset, charges, length, length, points_simulated, y_done, max_length, max_simulation_steps,
                               impossible_cutoff, max_delta_t_error,
                               step, 0, true, amin, amax, gravity, drag, num_iterations, num_elements, check_impossible, lambertW);
    }

    double starting_x = length;
    for (int y = start_pos+1; y < max_height_above; y+=num_threads) {
        starting_x = calculate_y_line(dataset, charges, starting_x, length, points_simulated, y_done, max_length, max_simulation_steps,
                                   impossible_cutoff, max_delta_t_error,
                                   step, y, true, amin, amax, gravity, drag, num_iterations, num_elements, check_impossible, lambertW);
    }

    starting_x = length;
    //y levels below cannon meanwhile can always be hit at some point, so just simulate until it hits a reachable point
    // and only then start calculating cutoff
    for (int y = start_pos-1; y > -max_height_below; y-=num_threads) {
        starting_x = calculate_y_line(dataset, charges, starting_x, length, points_simulated, y_done, max_length, max_simulation_steps,
                                   impossible_cutoff, max_delta_t_error,
                                   step, y, true, amin, amax, gravity, drag, num_iterations, num_elements, check_impossible, lambertW);
    }

    *done = true;
}

auto make_dataset(
                  int32_t charges,
                  int length,
                  int max_height_above = 256,
                  int max_height_below = 256,
                  int num_threads=16,
                  bool verbose=true,
                  int max_steps=100000,
                  int max_length=600,
                  float step=1,
                  uint32_t impossible_cutoff = 50,
                  float max_delta_t_error = 1,
                  int amin=-30, int amax=60, float gravity=0.05, float drag = 0.99,
                  int num_iterations=5, int num_elements=20,
                  bool check_impossible=true,
                  bool lambertW = false
                  ) {
    using namespace std::chrono_literals;

    std::vector<std::vector<std::tuple<int32_t, int32_t, float, float, float>>> threads_result;
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
                &threads_result[i],
                charges, length, max_height_above,
                max_height_below, i, num_threads,
                &threads_progress[i].first,
                &threads_progress[i].second,
                max_steps,
                max_length,
                impossible_cutoff,
                max_delta_t_error,
                step,
                amin, amax, gravity, drag, num_iterations, num_elements, check_impossible, lambertW,
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
