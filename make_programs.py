MODEM_CHANNEL = 100
HEADER_NAME = "bheader"
ROUNDING = 4
RUNTIME_COMBINE_ITEMS = True

import pickle
import os

for path in ["compiled/", "compiled/dispatcher", "compiled/shard"]:
    try:os.mkdir(path)
    except:pass

max_str_length = int(0.5 * 1000 * 1000)

recombine_str = "true" if RUNTIME_COMBINE_ITEMS else "false"

def smart_to_str(it):
    if it == int(it): return str(int(it))
    if 1 > it > 0:
        return str(it)[1:]
    return str(it)

def make_line(line_data, rounding=ROUNDING):
    str_line = "{"
    str_line += str(line_data[0][0]) + ","  # add info about distance offset to first item

    for item in line_data:
        for it in item:
            str_line += smart_to_str(round(it, rounding)) + ","

    str_line = str_line[:-1]
    str_line += "},"
    return str_line

print("Loading data")
with open("data", mode="rb") as file:
    data = pickle.load(file)
with open("shard.lua", mode="r") as file:
    shard_str = file.read()
with open("dispatcher.lua", mode="r") as file:
    dispatcher_str = file.read()

shard_str = shard_str.replace("REPLACE_THIS_WITH_MODEM_CHANNEL", str(MODEM_CHANNEL))
shard_str = shard_str.replace("REPLACE_THIS_WITH_HEADER_NAME", HEADER_NAME)

dispatcher_str = dispatcher_str.replace("REPLACE_THIS_WITH_HEADER_NAME", HEADER_NAME)
dispatcher_str = dispatcher_str.replace("REPLACE_THIS_WITH_MODEM_CHANNEL", str(MODEM_CHANNEL))

lines = {}

sorted_keys = list(data.keys())
sorted_keys.sort()

for key in sorted_keys:
    item = data[key]
    line = make_line(item)
    lines[key] = line

shard_data = ""
displacement = abs(list(lines.keys())[0])+1

min_y = min(lines.keys())
boundaries = []

def write_data():
    global file, min_y, shard_data
    boundaries.append([min_y+displacement, key])
    print(f"Writing compiled/{HEADER_NAME}{min_y+displacement}")
    with open(f"compiled/{HEADER_NAME}{min_y+displacement}", mode="w") as file:
        file.write("{"+shard_data[:-1]+"}")
    min_y = key
    shard_data = ""

print("Starting file creation")
for key in sorted_keys:
    line = lines[key]

    if len(shard_str) + len(line) + 2 > max_str_length:
        raise RuntimeError("Line length is bigger than max_str_length")

    if len(shard_str) + len(shard_data) + len(line) + 2 > max_str_length:
        write_data()
    shard_data += line

if len(shard_data) != 0:
    write_data()

needed_headers_str = "{"
for item in boundaries:
    needed_headers_str += f"\"{str(item[0])}\","
needed_headers_str += "}"

dispatcher_str = dispatcher_str.replace("REPLACE_THIS_WITH_NEEDED_HEADERS", needed_headers_str)
dispatcher_str = dispatcher_str.replace("REPLACE_THIS_WITH_Y_DISPLACEMENT", str(displacement))
dispatcher_str = dispatcher_str.replace("REPLACE_THIS_WITH_IF_RECOMBINE", recombine_str)

with open("compiled/dispatcher/startup.lua", mode="w") as file:
    file.write(dispatcher_str)

with open("compiled/shard/startup.lua", mode="w") as file:
    file.write(shard_str)