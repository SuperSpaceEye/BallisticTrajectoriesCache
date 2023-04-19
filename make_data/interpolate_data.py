import numpy as np


def interpolate_line(line):
    y = line[0][0][1]
    x_arr = [it[0][0] for it in line]
    min_x = min(x_arr)
    max_x = max(x_arr)

    length = int(max(abs(min_x), abs(max_x)) - min(abs(min_x), abs(max_x)))+1

    results = [it[1] for it in line]

    a1 = np.array([np.nan for _ in range(length)], dtype=float)
    a2 = np.array([np.nan for _ in range(length)], dtype=float)
    a3 = np.array([np.nan for _ in range(length)], dtype=float)

    for x, res in zip(x_arr, results):
        x = int(x - min_x)
        a1[x] = res[0]
        a2[x] = res[1]
        a3[x] = res[2]

    ok1 = ~np.isnan(a1)
    ok2 = ~np.isnan(a2)
    ok3 = ~np.isnan(a3)

    xp1 = ok1.ravel().nonzero()[0]
    xp2 = ok2.ravel().nonzero()[0]
    xp3 = ok3.ravel().nonzero()[0]

    fp1 = a1[~np.isnan(a1)]
    fp2 = a2[~np.isnan(a2)]
    fp3 = a3[~np.isnan(a3)]

    x1 = np.isnan(a1).ravel().nonzero()[0]
    x2 = np.isnan(a2).ravel().nonzero()[0]
    x3 = np.isnan(a3).ravel().nonzero()[0]

    a1[np.isnan(a1)] = np.interp(x1, xp1, fp1)
    a2[np.isnan(a2)] = np.interp(x2, xp2, fp2)
    a3[np.isnan(a3)] = np.interp(x3, xp3, fp3)

    newresults = [(i1, i2, i3) for i1, i2, i3 in zip(a1, a2, a3)]
    newinputs = [(x, y, 0) for x in range(int(x_arr[0]), int(x_arr[0])+length)]

    newline = [(it1, it2) for it1, it2 in zip(newinputs, newresults)]

    return newline


def transform_data(data):
    collapsed_data = []
    for thread_result in data:
        collapsed_data += thread_result

    lines = {}

    for item in collapsed_data:
        x, y, z = item[0]

        y = int(y)

        if y not in lines:
            lines[y] = []

        lines[y].append(item)

    for k in lines:
        lines[k] = interpolate_line(lines[k])

    return lines