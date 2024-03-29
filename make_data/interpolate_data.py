from tqdm import tqdm
import numpy as np
import gc

def interpolate_line(line):
    x_arr = [it[0] for it in line]
    min_x = min(x_arr)
    max_x = max(x_arr)

    length = int(max(abs(min_x), abs(max_x)) - min(abs(min_x), abs(max_x)))+1

    a1 = np.array([np.nan for _ in range(length)], dtype=float)
    a2 = np.array([np.nan for _ in range(length)], dtype=float)
    a3 = np.array([np.nan for _ in range(length)], dtype=float)
    a4 = np.array([np.nan for _ in range(length)], dtype=float)
    a5 = np.array([np.nan for _ in range(length)], dtype=float)
    a6 = np.array([np.nan for _ in range(length)], dtype=float)

    for x, res in zip(x_arr, line):
        x = int(x - min_x)
        a1[x] = res[1]
        a2[x] = res[2]
        a3[x] = res[3]
        a4[x] = res[4]
        a5[x] = res[5]
        a6[x] = res[6]

    ok1 = ~np.isnan(a1)
    ok2 = ~np.isnan(a2)
    ok3 = ~np.isnan(a3)
    ok4 = ~np.isnan(a4)
    ok5 = ~np.isnan(a5)
    ok6 = ~np.isnan(a6)

    xp1 = ok1.ravel().nonzero()[0]
    xp2 = ok2.ravel().nonzero()[0]
    xp3 = ok3.ravel().nonzero()[0]
    xp4 = ok4.ravel().nonzero()[0]
    xp5 = ok5.ravel().nonzero()[0]
    xp6 = ok6.ravel().nonzero()[0]

    fp1 = a1[~np.isnan(a1)]
    fp2 = a2[~np.isnan(a2)]
    fp3 = a3[~np.isnan(a3)]
    fp4 = a4[~np.isnan(a4)]
    fp5 = a5[~np.isnan(a5)]
    fp6 = a6[~np.isnan(a6)]

    x1 = np.isnan(a1).ravel().nonzero()[0]
    x2 = np.isnan(a2).ravel().nonzero()[0]
    x3 = np.isnan(a3).ravel().nonzero()[0]
    x4 = np.isnan(a4).ravel().nonzero()[0]
    x5 = np.isnan(a5).ravel().nonzero()[0]
    x6 = np.isnan(a6).ravel().nonzero()[0]

    a1[np.isnan(a1)] = np.interp(x1, xp1, fp1)
    a2[np.isnan(a2)] = np.interp(x2, xp2, fp2)
    a3[np.isnan(a3)] = np.interp(x3, xp3, fp3)
    a4[np.isnan(a4)] = np.interp(x4, xp4, fp4)
    a5[np.isnan(a5)] = np.interp(x5, xp5, fp5)
    a6[np.isnan(a6)] = np.interp(x6, xp6, fp6)

    newresults = [(i1, i2, i3, i4, i5, i6) for i1, i2, i3, i4, i5, i6 in zip(a1, a2, a3, a4, a5, a6)]
    newinputs = [x for x in range(int(x_arr[0]), int(x_arr[0])+length)]

    newline = [(it1, *it2) for it1, it2 in zip(newinputs, newresults)]

    return newline


def transform_data(data, do_interpolate=False):
    print("Collapsing data from threads")
    collapsed_data = []
    for thread_result in data:
        collapsed_data += thread_result
    del data
    gc.collect()

    print("Putting data into map")
    lines = {}
    for i, it in tqdm(enumerate(collapsed_data)):
        y = int(it[1])

        if y not in lines:
            lines[y] = []

        lines[y].append((it[0], it[2], it[3], it[4], it[5], it[6], it[7]))

    if do_interpolate:
        keys = list(lines.keys())
        keys.sort()
        print("Checking for interpolation")
        for k in keys:
            line = lines[k]
            start = line[0][0]
            is_cont = True
            for item in line:
                if start+1 < item[0]:
                    is_cont = False
                    print(f"Interpolating line {k}")
                    break
                start += 1
            if not is_cont:
                lines[k] = interpolate_line(lines[k])

    for k in lines:
        lines[k] = np.array(lines[k])

    return lines
