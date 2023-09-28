import numpy as np
from pycine.raw import read_frames
from ._datafuncs import extract_flucts, sav_smooth


class FastCam:
    def __init__(self, path):
        self.raw_img_gen, self.setup, self.bpp = read_frames(path, start_frame=1)
        self.img_array = np.array(list(self.raw_img_gen))
        self.number_of_frames = self.img_array.shape[0]
        self.frame_rate = self.setup.FrameRate
        self.post_trigger = self.setup.PostTrigger
        self.pre_trigger = self.number_of_frames - self.post_trigger
        self.time = (
            np.arange(self.number_of_frames) - self.pre_trigger
        ) / self.frame_rate
        self.time_ms = self.time * 1e3

    def extract_fluctuations(self):
        self.flucts = extract_flucts(
            self.img_array, sav_smooth(self.img_array, 25, axis=0)
        )
