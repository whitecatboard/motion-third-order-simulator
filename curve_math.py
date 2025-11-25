import sys
import numpy as np
import math

def solve_third_order_newton(f, df, x0, err):
    x1 = x0 - (f(x0) / df(x0))
    error = abs(x1 - x0)
    prev_error = sys.float_info.max

    while (error > err) and (error != prev_error):
        x0 = x1
        x1 = x0 - (f(x0) / df(x0))
        prev_error = error
        error = abs(x1 - x0)
        if (error > prev_error):
            x1 = np.nan
            
    return x1

def solve_third_order_newton(a, b, c, d, x0, err):
    f = a * x0**3 + b * x0**2 + c * x0 + d
    df = 3 * a * x0**2 + 2 * b * x0 + c

    x1 = x0 - (f / df)
    error = abs(x1 - x0)
    prev_error = sys.float_info.max

    while (error > err) and (error != prev_error):
        x0 = x1

        f = a * x0**3 + b * x0**2 + c * x0 + d
        df = 3 * a * x0**2 + 2 * b * x0 + c

        x1 = x0 - (f / df)
        prev_error = error
        error = abs(x1 - x0)
        if (error > prev_error):
            x1 = np.nan
            
    return x1

def solve_second_order_newton(a, b, c, x0, err):
    f = a * x0**2 + b * x0 + c
    df = 2 * a * x0 + 2 * b

    x1 = x0 - (f / df)
    error = abs(x1 - x0)
    prev_error = sys.float_info.max

    while (error > err) and (error != prev_error):
        x0 = x1

        f = a * x0**2 + b * x0 + c
        df = 2 * a * x0 + 2 * b

        x1 = x0 - (f / df)
        prev_error = error
        error = abs(x1 - x0)
        if (error > prev_error):
            x1 = np.nan
            
    return x1

def solve_second_order_pos(a, b, c):
    discriminant = b ** 2 - 4 * a * c
    unknown = np.nan

    if (discriminant > 0):
        unknown = (-b + math.sqrt(discriminant)) / (2 * a)

        if (unknown < 0):
            unknown = np.nan

    elif (discriminant == 0):
        unknown = - (b / (2 * a))
        if (unknown < 0):
            unknown = np.nan
    else:
        unknown = np.nan
    
    return unknown



