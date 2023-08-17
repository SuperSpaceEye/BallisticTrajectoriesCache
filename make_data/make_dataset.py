from BallisticFunctions import make_dataset
import pickle

if __name__ == "__main__":
    res_n = make_dataset(8, 32, 256, 256, 16, True, 100000, 600, 1, 500, 1, check_impossible=True)
    with open(f"../data", mode="wb") as file:
        pickle.dump(res_n, file)

import display_data