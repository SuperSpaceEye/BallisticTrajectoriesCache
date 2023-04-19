import matplotlib.pyplot as plt
import pickle
import interpolate_data

with open("../data", mode="rb") as file:
    data = pickle.load(file)

data = interpolate_data.transform_data(data)

x_axis = []
y_axis = []
z_axis = []

# for thread_result in data_n:
for key in data:
    line = data[key]
    for item in line:
        x_axis.append(item[0][0])
        y_axis.append(item[0][1])
        z_axis.append(item[1][1])

plt.scatter(x_axis, y_axis, c=z_axis)
plt.show()