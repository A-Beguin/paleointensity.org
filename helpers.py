import math
import helpers
import random

"""
Calculations that are used throughout the scripts are placed in this helpers script
"""

def cart2dir(x, y, z):
    """
    Function for converting from Cartesian coordinates to directions

    input: Cartesian coordinates x, y, z
    output: declination inclination and radius
    """
    R = math.sqrt(x**2 + y**2 + z**2)
    Dec = math.atan2(y, x) * 180 / math.pi
    Inc = math.asin(z / R) * 180 / math.pi

    if Dec < 0:
        Dec = Dec + 360
    elif Dec > 360:
        Dec = Dec - 360

    return [Dec, Inc, R]


def dir2cart(Dec, Inc, Mag):
    """
    Function for converting the direction into Cartesian coordinates

    input:  declination inclination and magnetization
    output: Cartesian coordinates x, y, z
    """
    if (Inc == -90) or (Inc == 90):
        x = 0.0
        y = 0.0
        z = math.sin(math.radians(Inc)) * Mag
    elif (Dec == -90) or (Dec == 90) or (Dec == 270):
        x = 0.0
        y = math.sin(math.radians(Dec)) * math.cos(math.radians(Inc)) * Mag
        z = math.sin(math.radians(Inc)) * Mag
    elif (Dec == -180) or (Dec == 180):
        x = math.cos(math.radians(Dec)) * math.cos(math.radians(Inc)) * Mag
        y = 0.0
        z = math.sin(math.radians(Inc)) * Mag
    else:
        x = math.cos(math.radians(Dec)) * math.cos(math.radians(Inc)) * Mag
        y = math.sin(math.radians(Dec)) * math.cos(math.radians(Inc)) * Mag
        z = math.sin(math.radians(Inc)) * Mag

    return [x, y, z]



def cross_product(a, b):
    """
    Function that gives the cross product of two lists of length 3

    input:  two lists, both of length 3
    output: the cross product, lists/vector of length 3
    """

    c = [a[1] * b[2] - a[2] * b[1],
         a[2] * b[0] - a[0] * b[2],
         a[0] * b[1] - a[1] * b[0]]

    return c



def dot_product(v1, v2):
    """
    Function that gives the dot product of two lists of any length

    input:  two lists, v1 and v2
    output: the dot product of the two vectors
    """
    multiply=[i * j for (i, j) in zip(v1, v2)]

    return sum(multiply)


def norm(a):
    """
    Function that returns the norm of a lists

    input:  lists of any length
    output: the norm of the lists
    """
    norm_is = math.sqrt(sum(list(comp**2 for comp in a)))
    return norm_is


def get_angle_diff(v1, v2):
    """
    Function that returns the angular difference in degrees between two lists

    input:  two lists, v1 and v2 of any length
    output: angular difference in degrees
    """

    angle = math.atan2(norm(cross_product(v1, v2)), dot_product(v1, v2))
    return math.degrees(angle)


def difference(v1, v2):
    """
    Function that returns the difference between two lists at component level
    sqrt( (x2-x1)^2 + (y2 - y1)^2 + (z2- z1)^2 )

    input:  two lists, v1 and v2 of any length
    output: difference, one number
    """

    diff = math.sqrt(sum([(j - i)**2 for (i, j) in zip(v1, v2)]))
    return diff


def list_min_num(l, num):
    """
    Function that returns a list for which a number is subtracted from the original list.

    input:  a list and a number
    output: list
    """
    result = []
    for i in range(len(l)):
        result.append(l[i] - num)
    return result


def list_div_num(l, num):
    """
    Function that returns a list for which the elements in the original list are divided by a number

    input:  a list and a number
    output: list
    """
    result = []
    for i in range(len(l)):
        result.append(l[i] / num)
    return result


def list_mult_num(l, num):
    """
    Function that returns a list for which the elements in the original list are multiplied with a number

    input:  a list of any length and a number
    output: list
    """
    result = []
    for i in range(len(l)):
        result.append(l[i] * num)
    return result


def transpose_list(l):
    """
    Function that returns the transpose of the input list

    input:  a list
    output: list
    """
    return [[row[i] for row in l] for i in range(len(l[0]))]


def list_min_list(l1, l2):
    """
    Function that subtracts the elements of one list from another list

    input:  two lists of any length
    output: list
    """
    result = []
    for i in range(len(l1)):
        result.append(l1[i] - l2[i])
    return result


def list_plus_list(l1, l2):
    """
    Function that adds the elements of one list to the elements from another list

    input:  two lists of any length
    output: list
    """
    result = []
    for i in range(len(l1)):
        result.append(l1[i] + l2[i])
    return result


def rand_num():
    """
    Function generates a random float with x: 0.0 <= x < 1.0

    input:  -
    output: float
    """

    return random.random()         # Random float:  0.0 <= x < 1.0
