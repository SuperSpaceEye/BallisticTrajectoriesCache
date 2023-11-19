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

h_x_axis = []
h_y_axis = []

l_x_axis = []
l_y_axis = []

h_pitch = []
h_airtime = []
l_pitch = []
l_airtime = []

# delta_t = []
# accuracy = []

print("Repacking data")

for key in tqdm(data.keys()):
    line = data[key]

    h_x_axis = np.append(h_x_axis, line[:, 0])
    h_y_axis += [int(key)] * len(line)

    h_pitch    = np.append(h_pitch, line[:, 2])
    h_airtime  = np.append(h_airtime, line[:, 3])



    filtered = line[line[:, 4] >= 0]
    l_x_axis = np.append(l_x_axis, filtered[:, 0])
    l_y_axis += [int(key)] * len(filtered)

    l_pitch    = np.append(l_pitch, filtered[:, 5])
    l_airtime  = np.append(l_airtime, filtered[:, 6])



    # delta_t  = np.append(delta_t, line[:, 1])
    # accuracy = np.append(accuracy, (1 - line[:, 1]/(line[:, 3]+1e-200)))
h_x_axis = h_x_axis.astype(int)

del data

print("Finished repacking")

min_y = min(h_y_axis)
min_x = min(h_x_axis)
max_y = max(h_y_axis)
max_x = max(h_x_axis)

h_y_range = max_y - min_y
h_x_range = max_x - min_x

min_y = min(l_y_axis)
min_x = min(l_x_axis)
max_y = max(l_y_axis)
max_x = max(l_x_axis)

l_y_range = max_y - min_y
l_x_range = max_x - min_x

h_figure_width = min(figure_width, h_x_range)
h_figure_height = min(figure_height, h_y_range)
l_figure_width = min(figure_width, h_x_range)
l_figure_height = min(figure_height, h_y_range)

# figure_width = x_range
# figure_height = y_range

print(f"Min y {min(h_y_axis)} Max y {max(h_y_axis)} Min x {min(h_x_axis)} Max x {max(h_x_axis)}")
print(f"Min high pitch {min(h_pitch)} Max high pitch {max(h_pitch)} "
      f"| Min high airtime {min(h_airtime)} Max high airtime {max(h_airtime)}")
print(f"Min low pitch {min(l_pitch)} Max low pitch {max(l_pitch)} "
      f"| Min low airtime {min(l_airtime)} Max low airtime {max(l_airtime)}")

h_pitch_agg    = make_agg(make_df(h_x_axis, h_y_axis, h_pitch), h_figure_width, h_figure_height)
h_airtime_agg  = make_agg(make_df(h_x_axis, h_y_axis, h_airtime), h_figure_width, h_figure_height)

l_pitch_agg    = make_agg(make_df(l_x_axis, l_y_axis, l_pitch), l_figure_width, l_figure_height)
l_airtime_agg  = make_agg(make_df(l_x_axis, l_y_axis, l_airtime), l_figure_width, l_figure_height)

fig, ax = plt.subplots(2, 2, figsize=(15,10))

ax[0, 0].title.set_text("high pitch")
ax[0, 1].title.set_text("high airtime")

ax[1, 0].title.set_text("low pitch")
ax[1, 1].title.set_text("low airtime")


sc1 = ax[0, 0].imshow(ds.tf.set_background(ds.tf.shade(h_pitch_agg, cmap=cc.bmy), "white").to_pil())
sc2 = ax[0, 1].imshow(ds.tf.set_background(ds.tf.shade(h_airtime_agg, cmap=cc.bmy), "white").to_pil())
sc3 = ax[1, 0].imshow(ds.tf.set_background(ds.tf.shade(l_pitch_agg, cmap=cc.bmy), "white").to_pil())
sc4 = ax[1, 1].imshow(ds.tf.set_background(ds.tf.shade(l_airtime_agg, cmap=cc.bmy), "white").to_pil())


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
