import datashader as ds
import pandas as pd
import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import pickle
import interpolate_data

figure_width = 1200
figure_height = 600

import matplotlib
matplotlib.use("GTK3Agg")

def make_df(X, Y, Z):
    return pd.DataFrame(np.array([X, Y, Z]).transpose((1, 0)), columns=["x", "y", "z"])

def make_agg(df, width=500, height=500):
    cvs = ds.Canvas(plot_width=width, plot_height=height)
    return cvs.points(df, "x", "y", agg=ds.mean("z"))


cmap = cc.cm["fire"].copy()
cmap.set_bad(cmap.get_under())  # set the color for 0

print("Loading data")
with open("../data", mode="rb") as file:
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

min_y = min(y_axis)
min_x = min(x_axis)
max_y = max(y_axis)
max_x = max(x_axis)

y_range = max_y - min_y
x_range = max_x - min_x

print(f"Min y {min(y_axis)} Max y {max(y_axis)} Min x {min(x_axis)} Max x {max(x_axis)}")

delta_t_agg  = make_agg(make_df(x_axis, y_axis, delta_t), 1200, 600)
pitch_agg    = make_agg(make_df(x_axis, y_axis, pitch), 1200, 600)
airtime_agg  = make_agg(make_df(x_axis, y_axis, airtime), 1200, 600)
accuracy_agg = make_agg(make_df(x_axis, y_axis, accuracy), 1200, 600)

fig, ax = plt.subplots(2, 2, figsize=(15,10))

ax[0, 0].title.set_text("delta_t")
ax[0, 1].title.set_text("pitch")
ax[1, 0].title.set_text("airtime")
ax[1, 1].title.set_text("accuracy")

sc1 = ax[0, 0].imshow(ds.tf.set_background(ds.tf.shade(delta_t_agg, cmap=cc.bmy), "white").to_pil())
sc2 = ax[0, 1].imshow(ds.tf.set_background(ds.tf.shade(pitch_agg, cmap=cc.bmy), "white").to_pil())
sc3 = ax[1, 0].imshow(ds.tf.set_background(ds.tf.shade(airtime_agg, cmap=cc.bmy), "white").to_pil())
sc4 = ax[1, 1].imshow(ds.tf.set_background(ds.tf.shade(accuracy_agg, cmap=cc.bmy), "white").to_pil())

def make_x_labels(labels):
    return [str(0)] + list(map(lambda x: str(int(x)), np.linspace(min_x, max_x, len(labels)-1)))

def make_y_labels(labels):
    return [str(0)] + list(map(lambda x: str(int(x)), np.linspace(max_y, min_y, len(labels)-1)))

def make_labels_for_ax(ax):
    ax.set_xticklabels(make_x_labels(ax.get_xticklabels()))
    ax.set_yticklabels(make_y_labels(ax.get_yticklabels()))

make_labels_for_ax(ax[0,0])
make_labels_for_ax(ax[0,1])
make_labels_for_ax(ax[1,0])
make_labels_for_ax(ax[1,1])

plt.show()