import matplotlib
matplotlib.use("GTK3Agg")

import matplotlib.pyplot as plt
import pickle
import interpolate_data

print("Loading data")
with open("../data_1", mode="rb") as file:
    data = pickle.load(file)

print("Transforming data")
data = interpolate_data.transform_data(data, do_interpolate=False)

x_axis = []
y_axis = []
delta_t = []
pitch = []
airtime = []
accuracy = []

print("Repacking data")
# for thread_result in data_n:
for key in data:
    line = data[key]
    for item in line:
        x_axis.append(item[0])
        y_axis.append(int(key))

        delta_t.append(item[1][0])
        # delta_t.append(math.log(item[1][0], 2))
        pitch.append(item[1][1])
        airtime.append(item[1][2])
        accuracy.append(1 - item[1][0]/item[1][2])
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