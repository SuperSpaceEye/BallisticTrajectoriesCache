from math import sin, cos, atan, sqrt, pi, radians, log
from numpy import linspace

def time_in_air(y0, y, Vy, max_steps=100000):
    """Find the air time of the projectile, using recursive sequence.
    It gives time in air by comparing the y-coordinate of the shell to the targets.
    The shell will hit target only when its y-coord is the same as the targets,
    so by comparing those we can find for how long the shell is airborne

    Args:
        y0 (float): y coordinate of the projectile
        y (float): y coordinate of the target
        Vy (float): vertical velocity of the projectile

    Returns:
        int: Airtime of projectile in ticks / "Error" if timeout
        :param max_steps:
    """
    t = 0
    t_below = 999_999_999_999

    if y0 <= y:
        # If cannon is lower than a target, simulating the way, up to the targets level
        while t < max_steps:
            y0p = y0
            y0 += Vy
            Vy = 0.99 * Vy - 0.05

            t += 1

            if y0 > y:  # Will break when the projectile gets higher than target
                t_below = t-1
                break

            # if projectile stopped ascending and didn't go above targets y pos
            if y0 - y0p < 0:
                return -1, -1

    while t < max_steps:
        y0 += Vy
        Vy = 0.99 * Vy - 0.05

        t += 1

        # Returns only when projectile is at same level as target or lower
        if y0 <= y:
            return t_below, t
    return t_below, -1


def calculate_if_pitch_hits(tried_pitch, initial_speed, length, distance, cannon, target, delta_t_max_overshoot=1,
                            max_steps=100000):
    # Bias that the cannon is probably gonna aim up instead of down
    # No use for now, useful for a later optimisation

    tp_rad = radians(tried_pitch)  # cuz triedPitch is in degrees

    Vw = cos(tp_rad) * initial_speed
    Vy = sin(tp_rad) * initial_speed

    x_coord_2d = length * cos(tp_rad)
    # This value is the horizontal distance between the mount and the tip
    # of the cannon, on the ballistic plane. By substracting this value from
    # the distance between the mount and the target, we get the horizontal distance
    # between the end of the barrel and the target which is what we want

    if Vw == 0: return None, False
    part = 1 - ((distance - x_coord_2d) / (100 * Vw))
    if part <= 0: return None, False
    horizontal_time_to_target = abs(log(part) / (-0.010050335853501))

    # This is the air resistance formula, here the denominator is ln(0.99)

    y_coord_of_end_barrel = cannon[1] + sin(tp_rad) * length

    t_below, t_above = time_in_air(y_coord_of_end_barrel, target[1], Vy, max_steps)

    if t_above < 0: return None, False
    if t_above < horizontal_time_to_target - delta_t_max_overshoot: return None, False

    delta_t = min(
        abs(horizontal_time_to_target - t_below),
        abs(horizontal_time_to_target - t_above)
    )
    # We calculate the difference between the time to target and airtime of the shell
    # The way this whole thing works is by comparing those values
    # We try to find the angle that corresponds the most to the timeToTarget

    """
    TimeToTarget is the time it takes for the shell to reach the target on the horizontal plane
    timeAir is the time it takes for the shell to reach the target depending on the given angle
    We try to find the airTime that corresponds the most to TimeToTarget, by taking the difference of TimeToTarget
    by every airTime possible (Bruteforcing every angle between -30 and 60 degrees)
    """

    return (delta_t, tried_pitch, delta_t + horizontal_time_to_target), True


def calculate_pitch(cannon, target, power, length, Dx, Dz, max_steps=100000, delta_t_max_overshoot=1):
    distance = sqrt(Dx * Dx + Dz * Dz)
    # Horizontal distance between target and mount

    pitch: float

    deltaTimes = []
    for triedPitch in range(60, -31, -1):
        items, is_successful = calculate_if_pitch_hits(triedPitch, power, length, distance, cannon, target,
                                                       delta_t_max_overshoot, max_steps)
        if not is_successful: continue
        deltaTimes.append((items[0], items[1]))

    if len(deltaTimes) == 0: return (-1, -1, -1)

    deltaTime, pitch = min(deltaTimes, key=lambda x: x[0])
    # Gives the minimum value depending on deltaT, the difference between airtime and TimeToTarget

    # We do the same thing, but near pitch, to get a more precise angle

    deltaTimes = []
    for triedPitch in linspace(pitch - 1, pitch + 1, 20):
        items, is_successful = calculate_if_pitch_hits(triedPitch, power, length, distance, cannon, target,
                                                       delta_t_max_overshoot, max_steps)
        if not is_successful: continue
        deltaTimes.append(items)

    if len(deltaTimes) == 0: return (-1, -1, -1)

    deltaTime, pitch, airtime = min(deltaTimes, key=lambda x: x[0])

    if pitch > 60 or pitch < -30: return (-1, -1, -1)

    return (deltaTime, pitch, airtime)


def calculate_yaw(Dx, Dz, direction):
    if Dx != 0:
        yaw = atan(Dz / Dx) * 180 / pi
    else:
        yaw = 90

    if Dx >= 0:
        yaw += 180

    directions = [90, 180, 270, 0]
    return (yaw + directions[direction]) % 360

def BallisticsToTarget(cannon, target, power, direction, R1, R2, length):
    """Function that calculates the elevation angle to hit the target with a cannon

    Args:
        cannon (tuple): Position of the cannon block held by the mount (x, y, z)
        target (tuple): Position of the target (x, y, z)
        power (int): Power of the cannon / Number of powder charges
        direction (str): Direction of the cannon (East, West...)
        R1 (int): Rotation speed of the yaw shaft
        R2 (int): Rotation speed of the pitch shaft
        length (int): Lenght of the cannon from mount to tip

    Returns:
        tuple: The yaw, pitch required and predicted airtime of the projectile
    """
    directions = ["north", "west", "south", "east"]
    if direction not in directions: return "Invalid direction"
    direction = directions.index(direction)

    Dx, Dz = (cannon[0] - target[0], cannon[2] - target[2])

    delta_t, pitch, airtime = calculate_pitch(cannon, target, power, length, Dx, Dz, 0, 0)
    yaw = calculate_yaw(Dx, Dz, direction)

    # Now, let's get the times we need to take in order to aim our cannon

    yawTime = yaw * 20 / (0.75 * R1)  # in TICKS
    pitchTime = pitch * 20 / (0.75 * R2)  # in TICKS
    fuzeTime = airtime + (delta_t / 2) - 10

    return (
        yaw,
        pitch,
        airtime,
        yawTime,
        pitchTime,
        fuzeTime,
    )


# print("For the cannon coordinates, please input the coordinates of the cannon mount.")
# cannonCoord = (
#     float(input("x coord of cannon : ")),
#     float(input("y coord of cannon : ")) + 2,
#     float(input("z coord of cannon : ")),
# )
# targetCoord = (
#     float(input("x coord of target : ")),
#     float(input("y coord of target : ")),
#     float(input("z coord of target : ")),
# )
#
# powderCharges = int(input("Number of powder charges (int) : "))
#
# directionOfCannon = input(
#     "What is the standart direction of the cannon ? (north, south, east, west) "
# )
#
# yawRPM = float(input("What is the RPM of the yaw axis ? "))
# pitchRPM = float(input("What is the RPM of the pitch axis ? "))
#
# cannonLenght = int(
#     input(
#         "What is the lenght of the cannon ? (From the block held by the mount to the tip of the cannon, the held block excluded) "
#     )
# )
#
# ballistic = BallisticsToTarget(
#     cannonCoord,
#     targetCoord,
#     powderCharges,
#     directionOfCannon,
#     yawRPM,
#     pitchRPM,
#     cannonLenght,
# )

# yawRPM = 4
# pitchRPM = 4
#
# ballistic = BallisticsToTarget((0, 0, 0), (300, 0, 0), 7, "north", 4, 4, 4)
#
# if type(ballistic) is tuple:
#     print(f"Yaw is {round(ballistic[0], 1)}")
#     print(f"Pitch is {round(ballistic[1], 1)}")
#     print(f"Airtime is {round(ballistic[2], 0)} ticks")
#     print(
#         f"With the yaw axis set at {yawRPM} rpm, the cannon must take {round(ballistic[3], 0)} ticks of turning the yaw axis."
#     )
#     print(
#         f"With the pitch axis set at {pitchRPM} rpm, the cannon must take {round(ballistic[4], 0)} ticks of turning the pitch axis."
#     )
#     print(
#         f"You must set the fuze time to {round(ballistic[5], 0)} ticks,\
#         which is {round(ballistic[5], 0) // 20} seconds and {round(ballistic[5], 0) % 20} ticks."
#     )
# else:
#     print(ballistic)

# CREDITS OF ORIGINAL FORMULAS : @sashafiesta#1978 on Discord

