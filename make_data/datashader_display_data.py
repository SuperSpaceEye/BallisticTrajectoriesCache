from tqdm import tqdm
import datashader as ds
import pandas as pd
import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import pickle

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
with open("../data_endal", mode="rb") as file:
    data = pickle.load(file)

x_axis = []
y_axis = []
delta_t = []
pitch = []
airtime = []
accuracy = []

print("Repacking data")

for key in tqdm(data.keys()):
    line = data[key]

    x_axis = np.append(x_axis, line[:, 0])
    y_axis += [int(key)] * len(line)

    delta_t  = np.append(delta_t, line[:, 1])
    pitch    = np.append(pitch, line[:, 2])
    airtime  = np.append(airtime, line[:, 3])
    accuracy = np.append(accuracy, (1 - line[:, 1]/(line[:, 3]+1e-200)))
x_axis = x_axis.astype(int)

del data

print("Finished repacking")

min_y = min(y_axis)
min_x = min(x_axis)
max_y = max(y_axis)
max_x = max(x_axis)

y_range = max_y - min_y
x_range = max_x - min_x

figure_width = min(figure_width, x_range)
figure_height = min(figure_height, y_range)

# figure_width = x_range
# figure_height = y_range

print(f"Min y {min(y_axis)} Max y {max(y_axis)} Min x {min(x_axis)} Max x {max(x_axis)}")
print(f"Min pitch {min(pitch)} Max pitch {max(pitch)} | Min airtime {min(airtime)} Max airtime {max(airtime)}")

# delta_t_agg  = make_agg(make_df(x_axis, y_axis, delta_t), figure_width, figure_height)
pitch_agg    = make_agg(make_df(x_axis, y_axis, pitch), figure_width, figure_height)
airtime_agg  = make_agg(make_df(x_axis, y_axis, airtime), figure_width, figure_height)
# accuracy_agg = make_agg(make_df(x_axis, y_axis, accuracy), figure_width, figure_height)

fig, ax = plt.subplots(2, 2, figsize=(15,10))

# ax[0, 0].title.set_text("delta_t")
ax[0, 1].title.set_text("pitch")
ax[1, 0].title.set_text("airtime")
# ax[1, 1].title.set_text("accuracy")

# sc1 = ax[0, 0].imshow(ds.tf.set_background(ds.tf.shade(delta_t_agg, cmap=cc.bmy), "white").to_pil())
sc2 = ax[0, 1].imshow(ds.tf.set_background(ds.tf.shade(pitch_agg, cmap=cc.bmy), "white").to_pil())
sc3 = ax[1, 0].imshow(ds.tf.set_background(ds.tf.shade(airtime_agg, cmap=cc.bmy), "white").to_pil())
# sc4 = ax[1, 1].imshow(ds.tf.set_background(ds.tf.shade(accuracy_agg, cmap=cc.bmy), "white").to_pil())

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
