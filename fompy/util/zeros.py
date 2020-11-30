from itertools import count, islice


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


def locate_first_function_zero(func, start, step):
    return next(locate_function_zeros(func, start, step), None)

def locate_nth_function_zero(func, start, step, num):
    return next(islice(locate_function_zeros(func, start, step), num, num+1), None)
