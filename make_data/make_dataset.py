from BallisticFunctions import make_dataset
from interpolate_data import transform_data
import pickle

if __name__ == "__main__":
    res = make_dataset(4, 32, 256, 256, 15, True, 1000000, 7000, 1, 7000, 1.5, check_impossible=True, lambertW=True,
                       # gravity=0.05*0.8, drag=1.0
                       amin=-90, amax=90
                       )
    res = transform_data(res)
    with open(f"../data_endal", mode="wb") as file:
        pickle.dump(res, file)

import display_data
