from BallisticFunctions import make_dataset_iterative, make_dataset_endals
from interpolate_data import transform_data
import pickle

if __name__ == "__main__":
    # res = make_dataset_iterative(initial_speed=16, length=32, amin=-90, amax=90, gravity=0.05, drag=0.99,
    #                              max_height_above=1024, max_height_below=1024, num_threads=15, verbose=True, max_steps=1000000,
    #                              max_distance=7000, step=1, impossible_cutoff=1000, max_delta_t_error=1.5,
    #                              check_impossible=True, lambertW=True)

    res = make_dataset_endals(
        initial_speed=10, length=32, amin=-90, amax=90, gravity=0.05, drag=0.99,
        max_height_above=1024, max_height_below=1024, num_threads=15, verbose=True,
        max_distance=1000, step=1, impossible_cutoff=700000
    )

    res = transform_data(res)
    with open(f"../data_endal", mode="wb") as file:
        pickle.dump(res, file)

import display_data
