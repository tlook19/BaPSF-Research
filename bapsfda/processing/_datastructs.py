__all__ = ["SignalArray", "get_board_props"]

from astropy import units as u # type: ignore
from typing import Optional
from dataclasses import dataclass, field
from bapsflib.lapd import File as lapdfile # type: ignore


@dataclass(frozen=True)
class BoardProperties:
    board_num: int
    channel_nums: list[int]
    num_shots: int
    num_samples: int
    adc_clk: float
    sw_sample_avg: Optional[int] = None
    hw_sample_avg: Optional[int] = None


@dataclass(frozen=True)
class DriveProperties:
    shots_per_pos: int
    ny: int
    nx: int
    nz: int


@dataclass(frozen=True)
class ChannelProperties:
    channel_num: int
    gain: float = 1
    resistance: Optional[float] = None


@dataclass(frozen=True)
class SISCrateChannel:
    crate_index: tuple[int, int]
    board_props: BoardProperties
    channel_props: ChannelProperties
    drive_props: DriveProperties

@dataclass
class DataRun:
    channels: list[SISCrateChannel]


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
            bprops = {
                "board_num": bnum,
                "channel_nums": [_ for _ in chans],
                "num_shots": props["nshotnum"],
                "num_samples": props["nt"],
                "adc_clk": props["clock rate"].to(u.Hz).value,
                "sw_sample_avg": props["shot average (software)"],
                "hw_sample_avg": props["sample average (hardware)"],
            }
            boards[f"board_{bnum}_props"] = BoardProperties(**bprops)
    return boards