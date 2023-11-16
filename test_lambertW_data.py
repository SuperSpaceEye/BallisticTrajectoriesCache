from make_data.interpolate_data import transform_data
import pickle
import gc
import tqdm
import numpy as np


print("Loading test data")
with open("data", mode="rb") as file:
    data_to_test = pickle.load(file)
print("Transforming test data")
data_to_test = transform_data(data_to_test, do_interpolate=False)

print("Loading validation data")
with open("data_n", mode="rb") as file:
    validation_data = pickle.load(file)
print("Transforming validation data")
validation_data = transform_data(validation_data, do_interpolate=False)

count = 0
no_validate = 0

pitch_error = np.array((), dtype=float)
delta_t_error = np.array((), dtype=float)
time_in_air_error = np.array((), dtype=float)

it = tqdm.tqdm(data_to_test)
for y in it:
    try:
        test_line = data_to_test[y]
        validate_line = validation_data[y]
    except:
        continue

    validate_map = {}
    for pos, item in validate_line:
        validate_map[pos] = item
    validate_line = None

    for pos, item in test_line:
        try:
            val_item = validate_map[pos]
        except:
            no_validate += 1
            continue

        pitch_error = np.append(pitch_error, (abs(val_item[0] - item[0])))
        delta_t_error = np.append(delta_t_error, abs(val_item[1] - item[1]))
        time_in_air_error = np.append(time_in_air_error, abs(val_item[2] - item[2]))

        count += 1

    validate_map = None
    gc.collect()

    it.set_postfix(
        {
            "count":count,
             "missed": no_validate,

             "avg_pitch_error":pitch_error.mean(),
             "avg_delta_t_error":delta_t_error.mean(),
             "avg_time_in_air_error":time_in_air_error.mean(),

            "median_pitch_error":np.median(pitch_error),
            "median_delta_t_error":np.median(delta_t_error),
            "median_time_in_air_error":np.median(time_in_air_error)
         })

print(f"Total error\npitch: {pitch_error.sum()} | delta_t {delta_t_error.sum()} | time_in_air {time_in_air_error.sum()}")
print(f"Avg error\npitch: {pitch_error.mean()} | delta_t {delta_t_error.mean()} | time_in_air {time_in_air_error.mean()}")
print(f"Median error\npitch: {np.median(pitch_error)} | delta_t {np.median(delta_t_error)} | time_in_air {np.median(time_in_air_error)}")
print(f"Missed validation {no_validate}")
print(f"Count: {count}")