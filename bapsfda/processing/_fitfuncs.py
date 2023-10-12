__all__ = ["supergaussian", "expfit"]

import numpy as np


def supergaussian(t, amp, mean, std, p=1):
    """supergaussian function of order p for curve fitting.
    p = 1 is a gaussian, p = 2 is a supergaussian, etc.
    p = 1/2 is a special case that works for fitting the temperature from the derivative
    of a langmuir sweep IV trace. In this case, Te = std / sqrt(2).

    Args:
        t (_type_): _description_
        amp (_type_): _description_
        mean (_type_): _description_
        std (_type_): _description_
        p (int, optional): _description_. Defaults to 1.

    Returns:
        _type_: _description_
    """
    if p == 1 / 2:
        return amp * np.exp(-abs(t - mean) / (np.sqrt(2) * std))
    else:
        return amp * np.exp(-((((t - mean) ** 2) / (2 * std**2)) ** p))


def expfit(t, tau, a, b, t0):
    """exponential function for curve fitting
    in the case of fitting an IV trace, tau is the electron temperature.

    Args:
        t (_type_): _description_
        amp (_type_): _description_
        tau (_type_): _description_
        offset (_type_): _description_

    Returns:
        _type_: _description_
    """
    return a * np.exp((t - t0) / tau) + b