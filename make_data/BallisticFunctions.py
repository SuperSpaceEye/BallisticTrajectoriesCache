try:
    from CannonBallisticFunctions import make_dataset_iterative, make_dataset_endals
except:
    raise NotImplementedError("Not implemented yet")
    # from fallback_py_funcs.py_make_dataset import make_dataset

_make_dataset_iterative = make_dataset_iterative
_make_dataset_endals = make_dataset_endals
f = float
i = int
b = bool


def make_dataset_iterative(initial_speed, length, amin=-30, amax=60, gravity=0.05, drag=0.99, max_height_above=256,
                           max_height_below=256, num_threads=16, verbose=True, max_steps=100000, max_distance=600, step=1,
                           impossible_cutoff=50, max_delta_t_error=1, num_iterations=5, num_elements=20,
                           check_impossible=True, lambertW=False):
    return _make_dataset_iterative(f(initial_speed), f(length), f(amin), f(amax), f(gravity), f(drag), f(max_delta_t_error),
                                   i(max_steps), i(num_iterations), i(num_elements), b(check_impossible), b(lambertW),
                                   i(max_height_above), i(max_height_below), i(num_threads), b(verbose), i(max_distance),
                                   f(step), i(impossible_cutoff))


def make_dataset_endals(initial_speed, length, amin=-30, amax=60, gravity=0.05, drag=0.99, mult_coeff=0.25,
                        acceptable_error_range=0.01, starting_from_t=0, max_mult_depth=8, starting_multiplier=2,
                        starting_max_depth=32, max_height_above=256, max_height_below=256, num_threads=16,
                        verbose=True, max_distance=600, step=1, impossible_cutoff=50):
    return _make_dataset_endals(
        f(initial_speed), f(length), f(amin), f(amax), f(gravity), f(drag),
        f(mult_coeff), f(acceptable_error_range), f(starting_from_t), f(max_mult_depth), f(starting_multiplier), f(starting_max_depth),
        i(max_height_above), i(max_height_below), i(num_threads), b(verbose), i(max_distance), f(step), i(impossible_cutoff))
