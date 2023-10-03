import datashader as ds
import pandas as pd
import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
from fast_histogram import histogram2d
import matplotlib.colors as colors
import pickle
import interpolate_data
import itertools as it

X = np.array([[x, y, np.random.randint(0, 50, 1)[0]] for x, y in it.product(range(50), range(50))])
# X = X.reshape((-1, 2))
# Z = np.random.random(1000000)

# cmap = cc.cm["fire"].copy()
# cmap.set_bad(cmap.get_under())  # set the color for 0
# bounds = [[X[:, 0].min(), X[:, 0].max()], [X[:, 1].min(), X[:, 1].max()]]
# h = histogram2d(X[:, 0], X[:, 1], range=bounds, bins=500)
#
# plt.imshow(h, norm=colors.LogNorm(vmin=1, vmax=h.max()), cmap=cmap)
# plt.axis('off')
# plt.colorbar()
# plt.show()

# df = pd.DataFrame(data=X, columns=["x", "y", "val"])  # create a DF from array
# cvs = ds.Canvas(plot_width=500, plot_height=500)  # auto range or provide the `bounds` argument
# agg = cvs.points(df, 'x', 'y', agg=ds.mean("val"))  # this is the histogram
# img = ds.tf.set_background(ds.tf.shade(agg), "white").to_pil()  # create a rasterized image
# plt.imshow(img)
# plt.axis('off')
# plt.show()


import matplotlib
matplotlib.use("GTK3Agg")

def make_df(X, Y, Z):
    return pd.DataFrame(np.array([X, Y, Z]).transpose((1, 0)), columns=["x", "y", "z"])

def make_agg(df, width=500, height=500):
    cvs = ds.Canvas(plot_width=width, plot_height=height)
    return cvs.points(df, "x", "y", agg=ds.mean("z"))


cmap = cc.cm["fire"].copy()
cmap.set_bad(cmap.get_under())  # set the color for 0

with open("../data", mode="rb") as file:
    data = pickle.load(file)

data = interpolate_data.transform_data(data)

x_axis = []
y_axis = []
delta_t = []
pitch = []
airtime = []
accuracy = []

# for thread_result in data_n:
for key in data:
    line = data[key]
    for item in line:
        x_axis.append(item[0][0])
        y_axis.append(item[0][1])

        delta_t.append(item[1][0])
        # delta_t.append(math.log(item[1][0], 2))
        pitch.append(item[1][1])
        airtime.append(item[1][2])
        accuracy.append(1 - item[1][0]/item[1][2])

delta_t_agg  = make_agg(make_df(x_axis, y_axis, delta_t), 1200, 600)
# pitch_agg    = make_agg(make_df(x_axis, y_axis, pitch))
# airtime_agg  = make_agg(make_df(x_axis, y_axis, airtime))
# accuracy_agg = make_agg(make_df(x_axis, y_axis, accuracy))


fig, ax = plt.subplots(2, 2, figsize=(15,10))

ax[0, 0].title.set_text("delta_t")
ax[0, 1].title.set_text("pitch")
ax[1, 0].title.set_text("airtime")
ax[1, 1].title.set_text("accuracy")

sc1 = ax[0, 0].imshow(ds.tf.set_background(ds.tf.shade(delta_t_agg, cmap=cc.bmy), "white").to_pil())

# sc1 = ax[0, 0].scatter(x_axis, y_axis, c=delta_t)
# sc2 = ax[0, 1].scatter(x_axis, y_axis, c=pitch)
# sc3 = ax[1, 0].scatter(x_axis, y_axis, c=airtime)
# sc4 = ax[1, 1].scatter(x_axis, y_axis, c=accuracy)

# plt.scatter(x_axis, y_axis, c=z_axis)

# plt.colorbar(sc1, ax=ax[0, 0])
# plt.colorbar(sc2, ax=ax[0, 1])
# plt.colorbar(sc3, ax=ax[1, 0])
# plt.colorbar(sc4, ax=ax[1, 1])

plt.show()