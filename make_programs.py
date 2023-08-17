import pickle
from make_data.interpolate_data import transform_data
from prepare_data import prepare_data

max_str_length = int(0.5*1024*1024)

def make_line(line_data, rounding=2):
    str_line = "{"
    for item in line_data:
        str_line+= "{" + ",".join([str(round(it, rounding)) for it in item[1]]) + "},"
    str_line += "},"
    return str_line


with open("accurate_data", mode="rb") as file:
    data = prepare_data(transform_data(pickle.load(file)))
with open("shard.lua", mode="r") as file:
    shard_str = file.read()
with open("dispatcher.lua", mode="r") as file:
    dispatcher_str = file.read()

lines = {}

sorted_keys = list(data.keys())
sorted_keys.sort()

for key in sorted_keys:
    item = data[key]
    line = make_line(item)
    lines[key] = line

shard_programs = []
new_shard = str(shard_str)
shard_data = ""
displacement = abs(list(lines.keys())[0])+1

min_y = list(lines.keys())[0]
boundaries = []

def make_shard(key):
    global min_y, new_shard, shard_data, displacement, boundaries
    boundaries.append([min_y, key])
    min_y = key

    new_shard = new_shard.replace("REPLACE_THIS_WITH_SHARD_NUM", str(len(shard_programs) + 1), 1)
    new_shard = new_shard.replace("REPLACE_THIS_WITH_SHARD_DISPLACEMENT", str(displacement), 1)
    new_shard = new_shard.replace("REPLACE_THIS_WITH_DATA", shard_data, 1)
    shard_programs.append(new_shard)
    new_shard = str(shard_str)
    shard_data = ""
 
    displacement = -key + 1

for key in sorted_keys:
    line = lines[key]

    if len(new_shard) + len(shard_data) + len(line) > max_str_length:
        make_shard(key)

    shard_data += line

make_shard(list(lines.keys())[-1]+1)

for item in boundaries:
    dispatcher_str = dispatcher_str.replace("REPLACE_THIS_WITH_NUM_SHARDS", str(len(boundaries)), 1)
    dispatcher_str = dispatcher_str.replace("REPLACE_THIS_WITH_BOUNDARIES", ",".join(
        ["{"+",".join([str(it) for it in item])+"}" for item in boundaries]
    ), 1)


for i, program in enumerate(shard_programs):
    with open(f"compiled/shard_{i + 1}.lua", mode="w") as file:
        file.write(program)

with open("compiled/dispatcher.lua", mode="w") as file:
    file.write(dispatcher_str)
