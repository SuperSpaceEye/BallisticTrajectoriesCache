import matplotlib
matplotlib.use("GTK3Agg")

import matplotlib.pyplot as plt
import pickle
import numpy as np
print("Loading data")
with open("../data", mode="rb") as file:
    data = pickle.load(file)

x_axis = []
y_axis = []
delta_t = []
pitch = []
airtime = []
accuracy = []

print("Repacking data")
for key in data:
    line = data[key]

    x_axis = np.append(x_axis, line[:, 0])
    y_axis += [int(key)] * len(line)

    delta_t  = np.append(delta_t, line[:, 1])
    pitch    = np.append(pitch, line[:, 2])
    airtime  = np.append(airtime, line[:, 3])
    accuracy = np.append(accuracy, (1 - line[:, 1]/(line[:, 3]+1e-200)))
x_axis = x_axis.astype(int)

print("Finished repacking")

fig, ax = plt.subplots(2, 2, figsize=(15,10))

ax[0, 0].title.set_text("delta_t")
ax[0, 1].title.set_text("pitch")
ax[1, 0].title.set_text("airtime")
ax[1, 1].title.set_text("accuracy")

sc1 = ax[0, 0].scatter(x_axis, y_axis, c=delta_t)
sc2 = ax[0, 1].scatter(x_axis, y_axis, c=pitch)
sc3 = ax[1, 0].scatter(x_axis, y_axis, c=airtime)
sc4 = ax[1, 1].scatter(x_axis, y_axis, c=accuracy)

# plt.scatter(x_axis, y_axis, c=z_axis)
plt.colorbar(sc1, ax=ax[0, 0])
plt.colorbar(sc2, ax=ax[0, 1])
plt.colorbar(sc3, ax=ax[1, 0])
plt.colorbar(sc4, ax=ax[1, 1])

plt.show()