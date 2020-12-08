from itertools import count, islice

from scipy.optimize import bisect


def locate_zeros(points):
    px = py = None
    for x, y in points:
        if px is not None and py * y < 0:
            yield px, x
        px = x
        py = y


def locate_function_zeros(func, start, step):
    xs = (start + step * i for i in count())
    points = ((x, func(x)) for x in xs)
    return locate_zeros(points)


def locate_nth_function_zero(func, start, step, num=0):
    if num < 0:
        step = -step
        num = -num - 1
    return next(islice(locate_function_zeros(func, start, step), num, num + 1), None)


def find_nth_function_zero(func, x0, xtol_coarse, xtol_fine, num=0):
    bracket = locate_nth_function_zero(func, x0, xtol_coarse, num)
    return bisect(func, *bracket, xtol=xtol_fine)
