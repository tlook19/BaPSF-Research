import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert # type: ignore
import pandas as pd # type: ignore
from scipy import constants as const # type: ignore

e = const.elementary_charge
m_e = const.electron_mass
eps0 = const.epsilon_0
c = const.speed_of_light


class Interferometer:
    def __init__(self, f_uwave, num_passes) -> None:
        calibration = 1.0 / (
            (num_passes / 4.0 / np.pi / f_uwave) * (e**2 / m_e / c / eps0)
        )
        self._inter_params = {
            "freq": f_uwave,
            "num_passes": num_passes,
            "cal_factor": calibration,
        }

    def _get_phase(x):
        hilb = hilbert(x)
        return np.unwrap(np.angle(hilb - np.mean(hilb)))

    def calculate_density(self, signal, ref_signal):
        dphase = self._get_phase(ref_signal) - self._get_phase(signal)
        return (dphase - dphase[0]) * self._inter_params["cal_factor"]
