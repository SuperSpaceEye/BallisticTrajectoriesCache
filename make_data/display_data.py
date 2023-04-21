import matplotlib.pyplot as plt
import pickle
import interpolate_data
import math

with open("../data2", mode="rb") as file:
    data = pickle.load(file)

data = interpolate_data.transform_data(data)

x_axis = []
y_axis = []
delta_t = []
pitch = []
airtime = []

# for thread_result in data_n:
for key in data:
    line = data[key]
    for item in line:
        # if item[1][0] > 5:
        #     continue

        x_axis.append(item[0][0])
        y_axis.append(item[0][1])

        # delta_t.append(item[1][0])
        delta_t.append(math.log(item[1][0], 2))
        pitch.append(item[1][1])
        airtime.append(item[1][2])

fig, ax = plt.subplots(3, figsize=(10,10))

sc1 = ax[0].scatter(x_axis, y_axis, c=delta_t)
sc2 = ax[1].scatter(x_axis, y_axis, c=pitch)
sc3 = ax[2].scatter(x_axis, y_axis, c=airtime)

# plt.scatter(x_axis, y_axis, c=z_axis)
plt.colorbar(sc1, ax=ax[0])
plt.colorbar(sc2, ax=ax[1])
plt.colorbar(sc3, ax=ax[2])

plt.show()