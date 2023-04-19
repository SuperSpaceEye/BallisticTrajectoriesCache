from CannonBallisticFunctions import make_dataset
import pickle


if __name__ == "__main__":
    # arguments:
    # charges, cannon length from base, max height above cannon (including 0),
    # max height below cannon, num_threads, verbose, max simulation steps,
    # max distance from cannon, x step, stop line after n impossible
    res_n = make_dataset(2, 8, 256, 256, 16, True, 100000, 600, 1, 100)
    with open("../data", mode="wb") as file:
        pickle.dump(res_n, file)

import display_data