__all__ = ["get_board_props", "SISRun"]

from astropy import units as u  # type: ignore
from typing import Optional
from dataclasses import dataclass, field
from bapsflib.lapd import File as lapdfile  # type: ignore
from abc import ABC, abstractmethod
import numpy as np


@dataclass(frozen=True)
class BoardProperties:
    board_num: int
    channel_nums: list[int]
    num_shots: int
    num_samples: int
    sample_frequeny: float
    time_step: float
    time_array: np.ndarray


@dataclass(frozen=True)
class DriveProperties:
    shots_per_pos: int
    ny: int
    nx: int
    nz: int


@dataclass(frozen=True)
class ChannelProperties:
    gain: float = 1
    resistance: Optional[float] = None
    units: type[u.UnitBase] = u.V


@dataclass
class SISCrateChannel:
    crate_index: tuple[int, int]
    board_props: Optional[BoardProperties] = None
    channel_props: Optional[ChannelProperties] = None
    drive_props: Optional[DriveProperties] = None


class DataRun(ABC):
    pass


class SISRun(DataRun):
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        with lapdfile(file_path) as f:
            self.board_props = get_board_props(f)
        self.channel_dict = {}
        self.num_shots = []
        for keybp, valuebp in self.board_props.items():
            bnum = valuebp.board_num
            self.num_shots.append(valuebp.num_shots)
            for cnum in valuebp.channel_nums:
                self.channel_dict["b_{bnum}_c_{cnum}"] = SISCrateChannel(
                    (bnum, cnum), valuebp
                )
        if all([_ == self.num_shots[0] for _ in self.num_shots]):
            self.num_shots = self.num_shots[0]
        else:
            raise ValueError("Number of shots varies across boards")
        self.drives = {}  # type: dict

    def config_drive(self, drive_name, shots_per_pos, ny, nx, nz):
        if shots_per_pos * ny * nx * nz != self.num_shots:
            raise ValueError(
                f"Number of shots ({self.num_shots}) does not match "
                f"shots_per_pos ({shots_per_pos}) * ny ({ny}) * nx ({nx}) * nz ({nz})"
            )
        self.drives[drive_name] = DriveProperties(shots_per_pos, ny, nx, nz)

    def set_chan_props(self, c_key, gain=1, resistance=None):
        if c_key not in self.channel_dict:
            raise ValueError(f"{c_key} not found in channel_dict")
        if resistance is not None:
            units = u.A
        else:
            units = u.V
        self.channel_dict[c_key].channel_props = ChannelProperties(
            gain, resistance, units
        )

    def set_drive_props(self, c_key, d_key):
        if d_key not in self.drives:
            raise ValueError(f"{d_key} not found in drives")
        if c_key not in self.channel_dict:
            raise ValueError(f"{c_key} not found in channel_dict")
        self.channel_dict[c_key].drive_props = self.drives[d_key]

    def check_config(self, bc_key):
        if self.channel_dict[bc_key].channel_props is None:
            raise ValueError(f"{bc_key} has no channel_props")
        if self.channel_dict[bc_key].drive_props is None:
            raise ValueError(f"{bc_key} has no drive_props")

    def extract_signal(self, bc_key):
        self.check_config(bc_key)
        with lapdfile(self.file_path) as f:
            bnum, cnum = self.channel_dict[bc_key].crate_index
            sig = f.get_data(bnum, cnum)["signal"]
        sig /= self.channel_dict[bc_key].channel_props.gain
        if self.channel_dict[bc_key].channel_props.resistance is not None:
            sig /= self.channel_dict[bc_key].channel_props.resistance
        return sig


def get_board_props(file: lapdfile) -> dict:
    """
    _summary_

    Parameters
    ----------
    file : lapdfile
        _description_

    Returns
    -------
    dict
        _description_
    """
    boards = {}
    digitizer = file.file_map.digitizers["SIS crate"]
    if digitizer is None:
        raise ValueError("No SIS Crate digitizer found in file")
    num_active_configs = len(digitizer.active_configs)
    if num_active_configs == 0:
        raise ValueError("No active SIS Crate config found in file")
    if num_active_configs > 1:
        raise ValueError(f"{num_active_configs} SIS Crate configs found in file")
    active_config = digitizer.active_configs[0]
    for adc in digitizer.configs[active_config]["adc"]:
        for slot in digitizer.configs[active_config][adc]:
            bnum, chans, props = slot
            nt = props["nt"]
            samp_freq = (
                props["clock rate"].to(u.Hz).value
                / props["sample average (hardware)"]
                / props["shot average (software)"]
            )
            dt = 1 / samp_freq
            time_array = np.arange(nt) * dt
            bprops = {
                "board_num": bnum,
                "channel_nums": [_ for _ in chans],
                "num_shots": props["nshotnum"],
                "num_samples": nt,
                "sample_frequeny": samp_freq,
                "time_step": dt,
                "time_array": time_array,
            }
            boards[f"board_{bnum}_props"] = BoardProperties(**bprops)
    return boards
