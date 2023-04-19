import math


def prepare_data(data):
    for k in data:
        line = data[k]
        for i, item in enumerate(line):
            coords = item[0]

            distance = math.sqrt(coords[0] * coords[0]
                                 + coords[2] * coords[2])

            line[i] = (int(distance), item[1])
    return data