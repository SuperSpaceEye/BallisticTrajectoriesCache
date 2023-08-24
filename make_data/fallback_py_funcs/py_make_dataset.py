import multiprocessing
import time
from .BallisticsCalculator import calculate_pitch
from multiprocessing import Process

def calculate_y_line(dataset, y, charges, barrel_length, points_simulated, y_done, max_length,
                     max_simulation_steps=100000, impossible_cutoff=50, delta_t_max_overshoot=1, step=1,
                     count_cutoff_at_the_start=False, amin=-30, amax=60, gravity=0.05, num_iterations=5,
                     num_elements=20, check_impossible=True):
    had_result = False
    cutoff_count = 0
    x = barrel_length
    while x < max_length:
        res1, res2 = calculate_pitch((0, 0, 0), (x, y, 0), charges, barrel_length, amin, amax, gravity,
                                     delta_t_max_overshoot, max_simulation_steps, num_iterations, num_elements,
                                     check_impossible)
        res = res1 if res2[0] > res1[0] >= 0 else res2
        if res[0] >= 0:
            dataset.append(((x, y, 0), res))
            had_result = True
        points_simulated[0]+=1

        if (had_result or count_cutoff_at_the_start) and res[0] < 0: cutoff_count+=1;
        else: cutoff_count = 0

        if cutoff_count >= impossible_cutoff: break

        x += step
    y_done[0]+=1

def make_dataset_thread(dataset, charges, length, max_height_above, max_height_below,
                        start_pos, num_threads, progress,
                        max_simulation_steps, max_length, impossible_cutoff,
                        delta_t_max_overshoot, step, done,
                        amin=-30, amax=60, gravity=0.05, num_iterations=5, num_elements=20,
                        check_impossible=True):
    for y in range(start_pos, max_height_above, num_threads):
        calculate_y_line(dataset, y, charges, length, progress[0], progress[1], max_length, max_simulation_steps,
                         impossible_cutoff, delta_t_max_overshoot, step, True, amin, amax, gravity, num_elements,
                         num_iterations, check_impossible)
        print(f"points calculated: {progress[0][0]} | y levels calculated: {progress[1][0]}")

    for y in range(start_pos-1, -max_height_below, -num_threads):
        calculate_y_line(dataset, y, charges, length, progress[0], progress[1], max_length, max_simulation_steps,
                         impossible_cutoff, delta_t_max_overshoot, step, False, amin, amax, gravity, num_elements,
                         num_iterations, check_impossible)
        print(f"points calculated: {progress[0][0]} | y levels calculated: {progress[1][0]}")

    done[0] = True

def make_dataset(charges, length, max_height_above=256, max_height_below=256,
                 num_threads=16, verbose=True, max_steps=100000, max_length=600,
                 step=1, impossible_cutoff=50, delta_t_max_overshoot=1,
                 amin=-30, amax=60, gravity=0.05, num_iterations=5, num_elements=20,
                 check_impossible=True):
    # m = multiprocessing.Manager()
    # dataset = [m.list() for i in range(num_threads)]
    # threads_progress = [m.list([[0], [0]]) for i in range(num_threads)]
    # done = [m.list([False]) for i in range(num_threads)]

    dataset = []
    done = [False]

    make_dataset_thread(dataset, charges, length, max_height_above, max_height_below, 0, 1, [[0], [0]],
                        max_steps, max_length, impossible_cutoff, delta_t_max_overshoot, step, done,
                        amin, amax, gravity, num_iterations, num_elements, check_impossible)

    return [dataset, ]

# def make_dataset(charges, length, max_height_above=256, max_height_below=256,
#                  num_threads=16, verbose=True, max_steps=100000, max_length=600,
#                  step=1, impossible_cutoff=50, delta_t_max_overshoot=1):
#     m = multiprocessing.Manager()
#     dataset = [m.list() for i in range(num_threads)]
#     threads_progress = [m.list([[0], [0]]) for i in range(num_threads)]
#     done = [m.list([False]) for i in range(num_threads)]
#     threads = []
#
#     for i in range(num_threads):
#         p = Process(target=make_dataset_thread,
#                     args=(dataset[i], charges, length,
#                           max_height_above, max_height_below, i, num_threads,
#                           threads_progress[i],
#                           max_steps, max_length, impossible_cutoff, delta_t_max_overshoot, step, done[i]))
#         p.start()
#         threads.append(p)
#     not_ended = True
#     while not_ended:
#         not_ended = False
#         time.sleep(1)
#
#         for is_done in done:
#             if not is_done[0]:
#                 not_ended += True
#                 break
#         if not verbose: continue
#
#         y_levels_done = 0
#         points_simulated = 0
#
#         for (x_done, y_done), t in zip(threads_progress, threads):
#             points_simulated += x_done[0]
#             y_levels_done += y_done[0]
#
#         print(f"points calculated: {points_simulated} | y levels calculated: {y_levels_done}")
#
#     return dataset