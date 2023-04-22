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

#include <pybind11/stl.h>

#include "linspace.h"

//calculates time when it goes above and then below target y pos.
inline std::pair<int64_t, int64_t> time_in_air(float y0, float y, float Vy, int32_t max_steps= 100000) {
    int64_t t = 0;
    int64_t t_below = INT64_MAX;


    if (y0 < y) {
        float y0_p;
        while (t < max_steps) {
            y0_p = y0;
            y0 += Vy;
            Vy = 0.99 * Vy - 0.05;

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
        Vy = 0.99 * Vy - 0.05;

        t += 1;

        if (y0 <= y) {return {t_below, t};}
    }

    return {t_below, -1};
}

inline double rad(auto deg) {return deg * (M_PI / 180);}

inline std::pair<std::array<float, 3>, bool> try_pitch(
        const std::array<float, 3> & cannon,
        const std::array<float, 3> & target,
        float tried_pitch,
        float distance,
        int32_t initial_speed,
        int32_t length,
        int32_t max_steps,
        float delta_t_max_overshoot = 1) {
    auto tp_rad = rad(tried_pitch);

    auto Vw = std::cos(tp_rad) * initial_speed;
    auto Vy = std::sin(tp_rad) * initial_speed;

    auto x_coord_2d = length * std::cos(tp_rad);

    if (Vw == 0) {return {{-1, -1, -1}, false};}
    auto part = 1 - (distance - x_coord_2d) / (100 * Vw);
    if (part <= 0) { return {{-1, -1, -1}, false};}
    auto horizontal_time_to_target = std::abs(std::log(part) / (-0.010050335853501));

    auto y_coord_end_of_barrel = cannon[1] + std::sin(tp_rad) * length;

    auto [t_below, t_above] = time_in_air(y_coord_end_of_barrel, target[1], Vy, max_steps);

    if (t_above < 0) { return {{-1, -1, -1}, false};}

    //as vertical and horizontal time to targets are calculated separately, it falsely calculates that cannon can hit
    // a target very far away by just launching it horizontally. to allow some inaccuracy i allow a bit of overshoot
    if (t_above < horizontal_time_to_target - delta_t_max_overshoot) { return {{-1, -1, -1}, false};}


    //if target is above cannon it may hit it on ascension
    auto delta_t = std::min(
            std::abs(horizontal_time_to_target-t_below),
            std::abs(horizontal_time_to_target-t_above)
            );

    return {{(float)delta_t, tried_pitch, (float)(delta_t+horizontal_time_to_target)}, true};
}

inline std::vector<std::pair<float, float>>
rough_pitch_estimation(const std::array<float, 3> & cannon,
                       const std::array<float, 3> & target,
                       float distance,
                       int32_t initial_speed,
                       int32_t length,
                       int32_t max_steps,
                       float delta_t_max_overshoot = 1) {
    std::vector<std::pair<float, float>> delta_times;
    for (int tried_pitch = 60; tried_pitch >= -30; tried_pitch--) {
        auto [items, is_successful] = try_pitch(cannon, target, tried_pitch, distance, initial_speed, length, max_steps,
                                                delta_t_max_overshoot);
        if (!is_successful) { continue;}
        delta_times.emplace_back(items[0], items[1]);
    }
    return delta_times;
}

inline std::vector<std::pair<float, float>>
py_rough_pitch_estimation(
                       pybind11::tuple & cannon_t,
                       pybind11::tuple & target_t,
                       float distance,
                       int32_t initial_speed,
                       int32_t length,
                       int32_t max_steps,
                       float delta_t_max_overshoot = 1) {
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

    return rough_pitch_estimation(cannon, target, distance, initial_speed, length, max_steps, delta_t_max_overshoot);
}

inline std::vector<std::array<float, 3>>
fine_pitch_estimation(const std::array<float, 3> & cannon,
                      const std::array<float, 3> & target,
                      float distance,
                      int32_t initial_speed,
                      int32_t length,
                      float pitch,
                      int32_t num_refined = 20,
                      int32_t max_steps = 100000,
                      float delta_t_max_overshoot = 1) {
    std::vector<std::array<float, 3>> delta_times;
    auto pitches = linspace<float>(pitch-1, pitch+1, num_refined);
    for (auto & tried_pitch: pitches) {
        auto [items, is_successful] = try_pitch(cannon, target, tried_pitch, distance, initial_speed, length, max_steps,
                                                delta_t_max_overshoot);
        if (!is_successful) { continue;}
        delta_times.push_back(items);
    }
    return delta_times;
}

inline std::vector<std::array<float, 3>>
py_fine_pitch_estimation(
                      pybind11::tuple & cannon_t,
                      pybind11::tuple & target_t,
                      float distance,
                      int32_t initial_speed,
                      int32_t length,
                      float pitch,
                      int32_t num_refined = 20,
                      int32_t max_steps = 100000,
                      float delta_t_max_overshoot = 1) {
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

    return fine_pitch_estimation(cannon, target, distance, initial_speed, length, pitch, num_refined, max_steps, delta_t_max_overshoot);
}

inline std::array<float, 3>
            calculate_pitch(const std::array<float, 3> & cannon,
                            const std::array<float, 3> & target,
                            int32_t power, int32_t length,
                            int32_t max_steps,
                            float delta_t_max_overshoot = 1) {
    auto Dx = cannon[0] - target[0];
    auto Dz = cannon[2] - target[2];
    auto distance = std::sqrt(Dx * Dx + Dz * Dz);
    auto initial_speed = power;

    auto delta_times1 = rough_pitch_estimation(cannon, target, distance, initial_speed, length, max_steps, delta_t_max_overshoot);

    if (delta_times1.empty()) {return std::array<float, 3>{-1, -1, -1};}

    auto min_pair = std::min_element(delta_times1.begin(), delta_times1.end(),
                                        [](const auto & a, const auto & b){ return a.first < b.first;})[0];
    auto delta_times2 = fine_pitch_estimation(cannon, target, distance, initial_speed, length, min_pair.second, 20, max_steps, delta_t_max_overshoot);
    if (delta_times2.empty()) {return std::array<float, 3>{-1, -1, -1};}

    auto min_arr = std::min_element(delta_times2.begin(), delta_times2.end(),
                                    [](const auto & a, const auto & b){return a[0] < b[0];})[0];

    if (min_arr[1] > 60 || min_arr[1] < -30) {return {-1, -1, -1};}

    return min_arr;
}

inline std::array<float, 3>
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

    return calculate_pitch(cannon, target, power, length, max_steps);
}

template<auto pitch_fn>
inline void calculate_y_line(std::vector<std::array<std::array<float, 3>, 2>> *dataset, int32_t charges, int barrel_length,
                      int *points_simulated, int *y_done, int max_simulation_steps, int max_length,
                      uint32_t impossible_cutoff, float delta_t_max_overshoot, float step, int y,
                      bool count_cutoff_at_the_start) {
    bool had_result = false;
    int cutoff_count = 0;
    for (float x = barrel_length; x < max_length; x += step) {
        auto res = pitch_fn({0, 0, 0}, {x, y, 0}, charges, barrel_length, max_simulation_steps, delta_t_max_overshoot);
        if (res[0] >= 0) {
            dataset->push_back(std::array<std::array<float, 3>, 2>{
                    std::array<float, 3>{(float) x, (float) y, 0},
                    res
            });
            had_result = true;
        }
        (*points_simulated)++;

        if ((had_result || count_cutoff_at_the_start) && res[0] < 0) {cutoff_count++;} else { cutoff_count = 0;}
        if (cutoff_count >= impossible_cutoff) { break;}
        }
    (*y_done)++;
}

template<auto pitch_fn>
auto make_dataset_thread(
                         std::vector<std::array<std::array<float, 3>, 2>> * dataset,
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
                         float delta_t_max_overshoot,
                         float step = 1,

                         uint8_t * done = nullptr
                                 ) {
    //y level above cannot may not be possible to reach at all, so to prevent simulating whole like cutoff starts early
    dataset->reserve(100000);
    for (int y = start_pos; y < max_height_above; y+=num_threads) {
        calculate_y_line<pitch_fn>(dataset, charges, length, points_simulated, y_done, max_simulation_steps, max_length,
                                   impossible_cutoff, delta_t_max_overshoot,
                                   step, y, true);
    }

    //y levels below cannon meanwhile can always be hit at some point, so just simulate until it hits a reachable point
    // and only then start calculating cutoff
    for (int y = -start_pos-1; y > -max_height_below; y-=num_threads) {
        calculate_y_line<pitch_fn>(dataset, charges, length, points_simulated, y_done, max_simulation_steps, max_length,
                                   impossible_cutoff, delta_t_max_overshoot,
                                   step, y, false);
    }

    *done = true;
}

template <auto pitch_fn>
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
                  float delta_t_max_overshoot = 1
                  ) {
    using namespace std::chrono_literals;

    std::vector<std::vector<std::array<std::array<float, 3>, 2>>> threads_result;
    std::vector<std::pair<int, int>> threads_progress;
    std::vector<std::thread> threads;
    std::vector<uint8_t> done;

    threads_progress.resize(num_threads);
    threads_result.resize(num_threads);
    threads.resize(num_threads);
    done.resize(num_threads, false);

    for (int i = 0; i < num_threads; i++) {
        threads[i] = std::thread(
                make_dataset_thread<pitch_fn>,
                &threads_result[i],
                charges, length, max_height_above,
                max_height_below, i, num_threads,
                &threads_progress[i].first,
                &threads_progress[i].second,
                max_steps,
                max_length,
                impossible_cutoff,
                delta_t_max_overshoot,
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




#endif //CANNONBALLISTICCALCULATOR_BALLISTIC_FUNCTIONS_H
