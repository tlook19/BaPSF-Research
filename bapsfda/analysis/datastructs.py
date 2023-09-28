import numpy as np
from astropy import units as u
from .datafuncs import butt_low, sav_smooth
from abc import ABC, abstractmethod


class SignalArray(ABC):
    pass


class DaqArray(SignalArray):
    """
    Class for signal arrays from the LAPD DAQ SIS Crate system.
    The LAPD DAQ generates a 2D array of signals with dimensions (ny * nx * nshots, nt) where nt is the number of time steps (samples) in the signal.
    The DAQ run properties are the ADC clock frequency, the number of samples averaged, the number of shots per position, and the number of probe positions in the x and y directions.
    The LAPD DAQ loop acquires nshots signals at each probe position, then moves to the next probe position, starting from the top left (ymax, xmin)
    and moving to the bottom right (ymin, xmax) stepping along x first, then y.

    This class may need to be adjusted for the XYZ drive, or if a new DAQ system is used.
    Args:
        SignalArray (ABC): Abstract Base Class for signal arrays
    """

    def __init__(
        self,
        signal: np.ndarray,
        adc_clk: float,
        sample_avg: int,
        ny: int,
        nx: int,
        nshots: int,
    ) -> None:
        """Initialize a DAQ array with a raw signal array and the DAQ run properties.

        Args:
            signal (np.ndarray): A raw signal array from the DAQ
            adc_clk (float): DAQ ADC clock frequency in Hz (1e6 for LAPD DAQ)
            sample_avg (int): Number of samples averaged by the DAQ (base 2 number for LAPD DAQ)
            ny (int): Number of probe positions in the y direction
            nx (int): Number of probe positions in the x direction
            nshots (int): Number of shots per probe position
        """
        self._signal_array = signal * u.V
        self._time_array = np.arange(self._nt) * self._dt
        self._run_props = {
            "sample_freq": adc_clk / sample_avg * u.Hz,
            "ny": ny,
            "nx": nx,
            "nshots": nshots,
            "dt": sample_avg / adc_clk * u.s,
            "nt": self._signal_array.shape[-1],
        }

    @property
    def signal_array(self) -> np.ndarray:
        return self._signal_array

    @property
    def time_array(self) -> np.ndarray:
        return self._time_array

    @property
    def run_props(self) -> dict:
        return self._run_props

    def shape(self) -> tuple:
        return self._signal_array.shape

    def reshape(self) -> np.ndarray:
        return self._signal_array.reshape(
            self._run_props["ny"],
            self._run_props["nx"],
            self._run_props["nshots"],
            self._run_props["nt"],
        )

    def low_pass_filter(self, cutoff, order=4) -> np.ndarray:
        return butt_low(
            self._signal_array, cutoff, self._run_props["sample_freq"], order
        )


class CurrentDaqArray(DaqArray):
    def __init__(self, gain, resistance) -> None:
        """Take a raw DAQ signal array and convert it to a current array in Amps.

        Args:
            gain (float): gain of the measurement circuit to be divided out
            resistance (float): resistance in ohms of the current sense resistor

        Raises:
            Exception: _description_
        """
        if resistance.unit != u.ohm:
            raise Exception("Sense resistance must be in ohms.")
        super().__init__()
        self._run_props["gain"] = gain
        self._run_props["resistance"] = resistance
        self._signal_array = (self._signal_array / (gain * resistance)).to(u.A)


class VoltageDaqArray(DaqArray):
    def __init__(self, gain) -> None:
        """Take a raw DAQ signal array and convert it to a voltage array in Volts.

        Args:
            gain (float): gain of the measurement circuit to be divided out
        """
        super().__init__()
        self._run_props["gain"] = gain
        self._signal_array = self._signal_array / gain
