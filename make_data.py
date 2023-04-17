import math

def make_data(data):
    collapsed_data = []
    for thread_result in data:
        collapsed_data += thread_result

    inputs = []
    outputs = []

    for item in collapsed_data:
        coords = item[0]

        distance = math.sqrt(coords[0] * coords[0]
                             + coords[2] * coords[2])

        height_difference = coords[1]

        inputs.append((int(height_difference), int(distance)))

        outputs.append(item[1])

    zipped_data = [(*inpt, oupt) for inpt, oupt in zip(inputs, outputs)]
    zipped_data.sort(key=lambda x: x[0])

    data = {}
    for item in zipped_data:
        if item[0] not in data:
            data[item[0]] = []
        data[item[0]].append(item[1:])

    return data


if __name__ == "__main__":
    with open("data", mode="rb") as file:
        import pickle
        data = pickle.load(file)
        make_data(data)