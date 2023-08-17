from BallisticFunctions import make_dataset
import pickle

if __name__ == "__main__":
    res_n = make_dataset(8, 32, 256, 256, 1, True, 100000, 200, 1, 50, 1, check_impossible=False)
    with open(f"../data", mode="wb") as file:
        pickle.dump(res_n, file)

import display_data