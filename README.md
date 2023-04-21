# WIP BallisticTrajectoriesCache

## What is it
This project allows getting needed pitch for [Create Big Cannons](https://www.curseforge.com/minecraft/mc-mods/create-big-cannons) in a constant time by precalculating every possible position that a cannon with given characteristics can hit. The result of these calculations then transformed into a lookup table divided between several data shards with lua programs for [cc: tweaked](https://www.curseforge.com/minecraft/mc-mods/cc-tweaked) that work with it.

## How to compile programs
0. 1. Install requirements.txt
   2. Install c++ compiler that supports at least c++17
   3. Install pybind11 for c++ https://github.com/pybind/pybind11
   4. In make_data run "python startup.py install" (It will compile and install it as a python package. if you don't want it, run "pip uninstall CannonBallisticFunctions" after finishing working with it. I will probably make it a pip package at some point.)
1. Generate data with make_dataset.py in make_data. You can view the result by viewing display_data.py
2. Create directory named "compiled" and run make_programs.py. It will create a dispatcher and shard programs. Default maximum size for a shard is 0.5 Mb, but you can configure it if needed
## How to use
After programs are compiled you now have dispatcher.lua and shard_{n}.lua files. dispatcher.lua is a main file that communicates with the shards. Modify dispatcher's wait_for_request with whatever communication logic you need.

To setup, connect dispatcher computer to shard computers by using networking cable and wired modem (data is not encrypted.

Dispatcher takes x/z distance to object with target_y-cannon_y as parameters. 