# @sashafiesta#1978 (Discord) : Original formulas and principles.
# @Malex#6461: Python adaptation, some changes and improvements on the original formulas. https://github.com/Malex21/CreateBigCannons-BallisticCalculator
# @SpaceEye#2191: Optimized even more.

from math import sin, cos, atan, sqrt, pi, radians, log
from numpy import linspace


#filter linspace
def flinspace(start, stop, num_elements, min, max):
    items = []
    for item in linspace(start, stop, num_elements):
        if item < min or item > max: continue
        items.append(item)
    return items


def get_root(d, from_end):
    if from_end:
        for i in reversed(range(0, len(d) - 1)):
            if d[i][0] > d[i + 1][0]: return d[i + 1]
        return d[0]
    else:
        for i in range(1, len(d)):
            if d[i - 1][0] < d[i][0]: return d[i - 1]
        return d[-1]


def time_in_air(y0, y, Vy, gravity=0.05, max_steps=100000):
    t = 0
    t_below = 999_999_999_999

    if y0 <= y:
        # If cannon is lower than a target, simulating the way, up to the targets level
        while t < max_steps:
            y0p = y0
            y0 += Vy
            Vy = 0.99 * Vy - gravity

            t += 1

            if y0 > y:  # Will break when the projectile gets higher than target
                t_below = t-1
                break

            # if projectile stopped ascending and didn't go above targets y pos
            if y0 - y0p < 0:
                return -1, -1

    while t < max_steps:
        y0 += Vy
        Vy = 0.99 * Vy - gravity

        t += 1

        # Returns only when projectile is at same level as target or lower
        if y0 <= y:
            return t_below, t
    return t_below, -1


def calculate_if_pitch_hits(tried_pitch, initial_speed, length, distance,
                            cannon, target, delta_t_max_overshoot=1,
                            gravity=0.05, max_steps=100000):
    tp_rad = radians(tried_pitch)

    Vw = cos(tp_rad) * initial_speed
    Vy = sin(tp_rad) * initial_speed

    x_coord_2d = length * cos(tp_rad)

    if Vw == 0: return None, False
    part = 1 - ((distance - x_coord_2d) / (100 * Vw))
    if part <= 0: return None, False
    horizontal_time_to_target = abs(log(part) / (-0.010050335853501)) # This is the air resistance formula, here the denominator is ln(0.99)

    y_coord_of_end_barrel = cannon[1] + sin(tp_rad) * length

    t_below, t_above = time_in_air(y_coord_of_end_barrel, target[1], Vy, gravity, max_steps)

    if t_below < 0: return None, False
    if t_above < horizontal_time_to_target - delta_t_max_overshoot: return None, False

    # if target is above cannon it may hit it on ascension
    delta_t = min(
        abs(horizontal_time_to_target - t_below),
        abs(horizontal_time_to_target - t_above)
    )

    return (delta_t, tried_pitch, delta_t + horizontal_time_to_target), True


def try_pitches(iter, *args):
    delta_times = []
    for try_pitch in iter:
        items, is_successful = calculate_if_pitch_hits(try_pitch, *args)
        if not is_successful: continue
        delta_times.append(items)
    return delta_times


def calculate_pitch(cannon, target, initial_speed, length, amin=-30, amax=60, gravity=0.05, delta_t_max_overshoot=1,
                    max_steps=100000, num_iterations=5, num_elements=20, check_impossible=True):
    Dx, Dz = cannon[0] - target[0], cannon[2] - target[2]
    distance = sqrt(Dx * Dx + Dz * Dz)

    delta_times = try_pitches(range(amax, amin-1, -1), initial_speed, length, distance, cannon, target,
                              delta_t_max_overshoot, gravity, max_steps)
    if len(delta_times) == 0: return (-1, -1, -1), (-1, -1, -1)

    dT1, p1, at1 = get_root(delta_times, False)
    dT2, p2, at2 = get_root(delta_times, True)

    c1 = True
    c2 = not p1 == p2  # calculate second
    same_res = p1 == p2  # if result is same

    for i in range(0, num_iterations):
        if c1: dTs1 = try_pitches(flinspace(p1 - 10**(-i), p1 + 10**(-i), num_elements, amin, amax), initial_speed, length, distance, cannon, target, delta_t_max_overshoot, gravity, max_steps)
        if c2: dTs2 = try_pitches(flinspace(p2 - 10**(-i), p2 + 10**(-i), num_elements, amin, amax), initial_speed, length, distance, cannon, target, delta_t_max_overshoot, gravity, max_steps)

        if c1 and len(dTs1) == 0: c1=False
        if c2 and len(dTs2) == 0: c2=False

        if not c1 and not c2: return (-1, -1, -1), (-1, -1, -1)

        if c1: dT1, p1, at1 = min(dTs1, key=lambda x: x[0])
        if c2: dT2, p2, at2 = min(dTs2, key=lambda x: x[0])

    if same_res: dT2, p2, at2 = dT1, p1, at1

    r1, r2 = (dT1, p1, at1), (dT2, p2, at2)
    if check_impossible and dT1 > delta_t_max_overshoot: r1 = (-1, -1, -1)
    if check_impossible and dT2 > delta_t_max_overshoot: r2 = (-1, -1, -1)

    return r1, r2


def calculate_yaw(Dx, Dz, direction):
    if Dx != 0:
        yaw = atan(Dz / Dx) * 180 / pi
    else:
        yaw = 90

    if Dx >= 0:
        yaw += 180

    directions = [90, 180, 270, 0]
    return (yaw + directions[direction]) % 360
