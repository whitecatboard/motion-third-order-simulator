import sys
import numpy as np
import math

def s_to_steps(s: float, alpha: float, epsilon: float) -> int:
    # Compute real step value
    r_stp = s * alpha

    # Compute natural step value
    stp = math.ceil(r_stp)

    if (stp - r_stp < epsilon):
        return stp
    else:
        return stp - 1

def solve_third_order_newton(a, b, c, d, x0, err):
    it = 1

    f = a * x0**3 + b * x0**2 + c * x0 + d
    df = 3 * a * x0**2 + 2 * b * x0 + c

    x1 = x0 - (f / df)
    error = abs(x1 - x0)
    prev_error = sys.float_info.max

    while (error > err) and (error != prev_error):
        it = it + 1
        x0 = x1

        f = a * x0**3 + b * x0**2 + c * x0 + d
        df = 3 * a * x0**2 + 2 * b * x0 + c
        
        x1 = x0 - (f / df)
        prev_error = error
        error = abs(x1 - x0)
        if (error > prev_error):
            x1 = np.nan

    if (it > 2):
        print(it)

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

def solve_second_order_neg(a, b, c):
    discriminant = b ** 2 - 4 * a * c
    unknown = np.nan

    if (discriminant > 0):
        unknown = (-b - math.sqrt(discriminant)) / (2 * a)

        if (unknown < 0):
            unknown = np.nan

    elif (discriminant == 0):
        unknown = - (b / (2 * a))
        if (unknown < 0):
            unknown = np.nan
    else:
        unknown = np.nan
    
    return unknown


