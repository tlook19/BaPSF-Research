import numpy as np
import matplotlib.pyplot as plt

# from ..processing._datastructs import DaqArray
from ..processing._datafuncs import butt_low, sav_smooth
from ..processing._fitfuncs import sgaussian_half, gaussian, expfit
from astropy import units as u  # type: ignore
from astropy.units.quantity import Quantity  # type: ignore
from astropy.constants import e, k_B, m_p  # type: ignore
from scipy.optimize import curve_fit  # type: ignore
from ._interferometer import uwave_calib_factor

__all__ = ["LangmuirSweep"]


class LangmuirSweep:
    @u.quantity_input
    def __init__(
        self,
        t_start: u.ms,
        t_ramp: u.ms,
        t_period: u.ms,
        ramp_sym: float,
        nsweeps: int,
    ):
        self._sweep_params = {
            "t_start": t_start.to(u.ms),
            "t_ramp_up": (t_ramp * ramp_sym).to(u.ms),
            "t_period": t_period.to(u.ms),
            "nsweeps": nsweeps,
        }

    def slice_sweeps(self, isweep, vsweep, t_array):
        if isweep.shape != vsweep.shape:
            raise Exception("Current and voltage signal arrays must be the same shape.")
        if isweep.dt != vsweep.dt:
            raise Exception("Current and voltage signals must have the same time step.")
        if t_array.shape[-1] != isweep.shape[-1]:
            raise Exception(
                "Time array must have the same number of samples as the signal arrays."
            )
        sample_freq = (1 / isweep.dt).to(u.Hz).value
        isweep_filtered = isweep
        vsweep_filtered = vsweep
        # removing DC offset from current signal
        isweep_filtered = np.moveaxis(
            np.moveaxis(isweep_filtered, 1, 0)
            - np.average(isweep_filtered[:, -100:-1], axis=-1),
            0,
            1,
        )
        # take 50 samples off each side of ramp to avoid edge effects
        t_start_ind = (
            int(self._sweep_params["t_start"].to(u.s).value * sample_freq) + 50
        )
        t_ramp_ind = (
            int(self._sweep_params["t_ramp_up"].to(u.s).value * sample_freq) - 100
        )
        t_period_ind = int(self._sweep_params["t_period"].to(u.s).value * sample_freq)
        ramp_times = (
            np.arange(self._sweep_params["nsweeps"]) * self._sweep_params["t_period"]
            + self._sweep_params["t_start"]
        )
        v_slices = np.empty(
            (vsweep_filtered.shape[0], self._sweep_params["nsweeps"], t_ramp_ind)
        )
        i_slices = np.empty(
            (isweep_filtered.shape[0], self._sweep_params["nsweeps"], t_ramp_ind)
        )
        for i in range(vsweep_filtered.shape[0]):
            for j in range(self._sweep_params["nsweeps"]):
                k = t_start_ind + j * t_period_ind
                ramp_slice = slice(k, k + t_ramp_ind)
                vsort = vsweep_filtered[i, ramp_slice]
                isort = isweep_filtered[i, ramp_slice]
                sort_ind = np.argsort(vsort)
                v_slices[i, j] = np.take_along_axis(vsort, sort_ind, axis=0)
                i_slices[i, j] = np.take_along_axis(isort, sort_ind, axis=0)
                if i_slices[i, j, 0] > i_slices[i, j, -1]:
                    i_slices[i, j] *= -1  # flip if IV trace is inverted
        return v_slices, i_slices, ramp_times

    def _plot_sweep_debug(self, vslice, islice, i2, arg_vf, arg_vp):
        plt.plot(vslice, islice)
        plt.plot(vslice, i2)
        plt.axvline(vslice[arg_vf], color="k")
        plt.axvline(vslice[arg_vp], color="k")
        plt.show()

    def _find_isat_vf(self, vslice, islice, t_pre_ind):
        arg_vf = islice.shape[0] - np.argmin(np.abs(islice[::-1]))
        vf = vslice[arg_vf]
        isat = np.average(islice[:t_pre_ind])
        isat_std = np.std(islice[:t_pre_ind])
        # a, b = curve_fit(
        #     lambda x, m, b: m * x + b, vslice[: arg_vf // 2], islice[: arg_vf // 2]
        # )
        # iline = a[0] * vslice + a[1]
        return isat, vf, isat_std, arg_vf

    def _find_plasma_potential(self, vslice, ivgrad):
        """Find the plasma potential for a single sweep via the first derivative method. Report Te estimate from super-gaussian fit.

        Args:
            vslice (1d array): Voltage slice sorted monotonically increasing.
            islice (1d array): Associated current slice.

        Returns:
            _type_: _description_
        """
        a1, b1 = curve_fit(sgaussian_half, vslice, ivgrad)
        vp = a1[1]
        te = a1[2] / np.sqrt(2)
        return vp, te

    def analyze_sweeps(self, v_slices, i_slices):
        """
        Analyze sweeps and store the plasma parameters in self.plasma_params.
        plasma_params index corresponds to the following:
        0: Plasma potential
        1: Electron temperature

        TODO: Add error handling for bad fits

        Args:
            v_slices (_type_): _description_
            i_slices (_type_): _description_

        Returns:
        """
        v_avg = np.mean(v_slices.reshape(61, 10, 28, 500), axis=1)
        i_avg = np.mean(i_slices.reshape(61, 10, 28, 500), axis=1)
        sorted_ind = np.argsort(v_avg, axis=2)
        v_avg = np.take_along_axis(v_avg, sorted_ind, axis=2)
        i_avg = np.take_along_axis(i_avg, sorted_ind, axis=2)
        i_vsm = sav_smooth(i_avg, b=50, axis=-1)
        plasma_params = np.empty((v_avg.shape[0], v_avg.shape[1], 2))
        for i in range(v_avg.shape[0]):
            for j in range(v_avg.shape[1]):
                ivgrad = np.gradient(i_vsm[i, j], v_avg[i, j])
                vp, te = self._find_plasma_potential(v_avg[i, j], ivgrad)
                plasma_params[i, j, 0] = vp
                plasma_params[i, j, 1] = te
        return plasma_params

    def calculate_density(
        self, plasma_params: np.ndarray, probe_area: Quantity, ion_mass_factor: float
    ):
        """
        Calculate the plasma density from the plasma parameters and probe area.

        Args:
            plasma_params (_type_): _description_
            probe_area (_type_): _description_

        Returns:
            _type_: _description_
        """
        isat = plasma_params[:, :, 0].to(u.A)  # type: ignore
        temp = (plasma_params[:, :, 3] * e / k_B).to(u.J)
        m_i = (ion_mass_factor * m_p).to(u.kg)
        area = probe_area.to(u.m**2)

        electron_density = (isat) / (e * area) * np.sqrt(m_i / temp)
        return electron_density.to(u.cm**-3)

    # TODO: Clean this up
    def calibrate_density(
        self,
        isats: np.ndarray,
        tes: np.ndarray,
        xpos: np.ndarray,
        ion_mass_factor: float,
        int_array: np.ndarray,
    ):
        ucal = uwave_calib_factor(288e9, 2)
        x = xpos * u.m
        i = isats * u.A
        cs = np.sqrt((tes * u.eV).to(u.J) / (ion_mass_factor * m_p))
        nish = np.exp(1 / 2) * i / (e.si * cs)
        nlid = np.trapz(nish, x, axis=0)
        smint = sav_smooth(int_array * ucal, b=40, axis=-1)
        t_int = np.round(((np.arange(4096) - 409.6) * 2.442e-5) * 1e3, 2)
        t_int_ramp_ind = np.arange(self._sweep_params["nsweeps"]) * 41 + 410
        area = (nlid / smint[t_int_ramp_ind]) * u.m**2
        return area, t_int[t_int_ramp_ind]
