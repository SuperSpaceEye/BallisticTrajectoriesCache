from BallisticFunctions import make_dataset
import pickle

if __name__ == "__main__":
    res_n = make_dataset(110, 32, 1024, 1024, 15, True, 100000, 11000, 1, 7000, 1, check_impossible=True,
                         # gravity=0.05*0.8, drag=1.0
                         )
    with open(f"../data", mode="wb") as file:
        pickle.dump(res_n, file)

import display_data