try:
    a = 1 / 0
    from CannonBallisticFunctions import make_dataset
except:
    from fallback_py_funcs.py_make_dataset import make_dataset
import pickle

if __name__ == "__main__":
    # arguments:
    # charges, cannon length from base, max height above cannon (including 0),
    # max height below cannon, num_threads, verbose, max simulation steps,
    # max distance from cannon, x step, stop line after n impossible,
    # max delta_t overshoot time
    res_n = make_dataset(8, 32, 16, 16, 16, True, 10000, 200, 1, 50, 1)
    with open(f"../data", mode="wb") as file:
        pickle.dump(res_n, file)

import display_data