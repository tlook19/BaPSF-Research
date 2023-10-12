__all__ = ["SignalArray", "DaqArray", "CurrentDaqArray", "VoltageDaqArray"]

import numpy as np
from astropy import units as u
from ._datafuncs import butt_low, sav_smooth
from abc import ABC, abstractmethod
from typing import Optional
from dataclasses import dataclass


@dataclass
class RunProperties:
    shots_total: int
    trigger_delay: float = 0.0


@dataclass
class DriveProperties:
    shots_per_pos: int
    nx: int = 1
    ny: int = 1
    nz: int = 1


@dataclass
class BoardProperties:
    samples_per_shot: int
    adc_clk: float = 100e6
    sample_avg: int = 1


@dataclass
class ChannelProperties:
    gain: float = 1
    resistance: Optional[float] = None


class SignalArray(ABC):
    pass


class DataRun:
    """
    Class for an LAPD data run with SiS Crate Digitizer DAQ system.
    The LAPD data run is a collection of signal arrays where the each "Board" has a base ADC clock frequency, number of samples averaged, and number of samples per shot.
    Each board has 8 channels the each digitize a voltage signal. Each signal may have some gain or resistance associated with it.
    A typical run will take a number of shots at each probe position, then move to the next probe position,
    starting from the top left (ymax, xmin) and scanning from left to right, then top to bottom.

    A drive is a collection of channels (which may be on different boards) that are moved together.

    Run Properties that stay for the entire run are:
        - shots_total: total number of shots in the run
        - trigger_delay: Delay between LAPD 1kA discharge current trigger and the DAQ trigger in seconds(usually 0)

    Drive Props that are the same for the entire datarun are nx * ny * nz * shots_per_pos should equal shots_total: (assumes grad)
        - shots_per_pos: number of shots at each probe position
        - ny: number of probe positions in the y direction
        - nx: number of probe positions in the x direction
        - nz: number of probe positions in the z direction

    Board properties are:
        - adc_clk: ADC clock frequency in Hz (100e6 for LAPD DAQ)
        - sample_avg: Number of samples averaged by the DAQ (base 2 number for LAPD DAQ)
        - ns: number of samples per shot

    Channel properties are:
        - gain: gain of the measurement circuit to be divided out
        - resistance(optional): resistance in Ohms of the current sense resistor

    From a composition standpoint a datarun is a collection of boards, and each board is a collection of channels.

    In the future, it may be relevent that an experiment is a collection of data runs.

    Returns:
        _type_: _description_
    """


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
        adc_clk: float = 100e6,
        sample_avg: int = 1,
        ny: int = 1,
        nx: int = 1,
        nshots: int = 1,
        gain: float = 1,
        resistance: Optional[float] = None,
    ) -> None:
        """Initialize a DAQ array with a raw signal array and the DAQ run properties.

        Args:
            signal (np.ndarray): A raw signal array from the DAQ
            adc_clk (float): DAQ ADC clock frequency in Hz (100e6 for LAPD DAQ)
            sample_avg (int): Number of samples averaged by the DAQ (base 2 number for LAPD DAQ)
            ny (int): Number of probe positions in the y direction
            nx (int): Number of probe positions in the x direction
            nshots (int): Number of shots per probe position
            gain (float): Gain of the measurement circuit to be divided out
            resistance (float): Resistance in ohms of the current sense resistor
        """
        self._signal_array = signal / gain
        self._units = u.V
        if resistance is not None:
            self._signal_array = self._signal_array / resistance
            self._units = u.A
        self._sample_freq = adc_clk / sample_avg * u.Hz
        self._ny = ny
        self._nx = nx
        self._nshots = nshots
        self._dt = (self._sample_freq**-1).to(u.s)
        self._nt = self._signal_array.shape[-1]
        self._units = u.V
        self._time_array = (np.arange(self._nt) * self._dt).to(u.s)

    @property
    def signal_array(self) -> np.ndarray:
        return self._signal_array

    @property
    def time_array(self) -> np.ndarray:
        return self._time_array

    @property
    def sample_freq(self) -> u.Quantity:
        return self._sample_freq

    @property
    def dt(self) -> u.Quantity:
        return self._dt

    @property
    def nt(self) -> int:
        return self._nt

    @property
    def ny(self) -> int:
        return self._ny

    @property
    def nx(self) -> int:
        return self._nx

    @property
    def nshots(self) -> int:
        return self._nshots

    @property
    def units(self) -> u.Unit:
        return self._units

    @property
    def shape(self) -> tuple:
        return self._signal_array.shape
