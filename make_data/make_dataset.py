from BallisticFunctions import make_dataset
import pickle

if __name__ == "__main__":
    res_n = make_dataset(16, 32, 1024, 256, 15, True, 1000000, 10000, 1, 7000, 1.5, check_impossible=True, lambertW=True,
                         # gravity=0.05*0.8, drag=1.0
                         )
    with open(f"../data", mode="wb") as file:
        pickle.dump(res_n, file)

import display_data