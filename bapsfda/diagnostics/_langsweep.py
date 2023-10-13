import numpy as np
import matplotlib.pyplot as plt
from ..processing._datastructs import CurrentDaqArray, VoltageDaqArray
from ..processing._datafuncs import butt_low, sav_smooth
from ..processing._fitfuncs import supergaussian, expfit
from astropy import units as u
from astropy.units.quantity import Quantity
from astropy.constants import e, k_B, m_p
from scipy.optimize import curve_fit

__all__ = ["LangmuirSweep"]


class LangmuirSweep:
    """
    Langmuir Sweep Class
    Initialize with sweeping parameters. This class can then be used to slice the sweeps and analyze them.
    """

    def __init__(
        self,
        t_start: Quantity,
        t_ramp: Quantity,
        t_period: Quantity,
        ramp_sym: float,
        nsweeps: int,
    ):
        """_summary_

        Args:
            t_start (Quantity): _description_
            t_ramp (Quantity): _description_
            t_period (Quantity): _description_
            ramp_sym (float): _description_
            nsweeps (int): _description_

        Raises:
            Exception: _description_
        """
        if t_start.unit != u.ms or t_ramp.unit != u.ms or t_period.unit != u.ms:
            raise Exception("Sweep t_ parameters must be in milliseconds.")
        self._sweep_params = {
            "t_start": t_start,
            "t_ramp_up": t_ramp * ramp_sym,
            "t_period": t_period,
            "nsweeps": nsweeps,
        }

    def slice_sweeps(self, isweep: CurrentDaqArray, vsweep: VoltageDaqArray):
        """
        Take arrays of current and voltage traces and slice them into individual sweeps
        according to the parameters set in the LangmuirSweep object. Returns arrays of the slices
        with dimensions (nshots, nsweeps, t_ramp_ind), and an array of times corresponding to the
        start of each sweep.

        Args:
            isweep (CurrentDaqArray): _description_
            vsweep (VoltageDaqArray): _description_

        Raises:
            Exception: _description_
            Exception: _description_
            Exception: _description_

        Returns:
            _type_: _description_
        """
        if isweep.shape != vsweep.shape:
            raise Exception("Current and voltage signal arrays must be the same shape.")
        if isweep.run_props["dt"] != vsweep.run_props["dt"]:
            raise Exception("Current and voltage signals must have the same time step.")
        if np.any(isweep.time_array != vsweep.time_array):
            raise Exception(
                "Current and voltage signals must have the same time array."
            )
        iu = isweep.signal_array.unit
        vu = vsweep.signal_array.unit
        isweep_filtered = butt_low(
            isweep.signal_array, 5e5, isweep.run_props["sample_freq"].value, order=4
        )
        vsweep_filtered = butt_low(
            vsweep.signal_array, 5e5, vsweep.run_props["sample_freq"].value, order=4
        )
        isweep_filtered = np.moveaxis(
            np.moveaxis(isweep_filtered, 1, 0)
            - np.average(isweep_filtered[:, -2000:-1], axis=-1),
            0,
            1,
        )
        t_pre_ind = int((10 * u.us).to(u.s) / isweep.run_props["dt"])
        t_start_ind = int(
            self._sweep_params["t_start"].to(u.s) / isweep.run_props["dt"] - t_pre_ind
        )  # start 10 us early
        t_ramp_ind = int(
            self._sweep_params["t_ramp_up"].to(u.s) / isweep.run_props["dt"] + t_pre_ind
        )  # add 10 us to compensate for the 10 us early start
        t_period_ind = int(
            self._sweep_params["t_period"].to(u.s) / isweep.run_props["dt"]
        )
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
                vsort = vsweep_filtered[i, k : k + t_ramp_ind]
                isort = isweep_filtered[i, k : k + t_ramp_ind]
                sort_ind = np.argsort(vsort)
                v_slices[i, j] = sav_smooth(
                    np.take_along_axis(vsort, sort_ind, axis=0), 25
                )
                i_slices[i, j] = sav_smooth(
                    np.take_along_axis(isort, sort_ind, axis=0), 25
                )
        return v_slices * vu, i_slices * iu, ramp_times, t_pre_ind

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

    def _find_plasma_potential(self, vslice, islice):
        """Find the plasma potential for a single sweep via the first derivative method. Report Te estimate from super-gaussian fit.

        Args:
            vslice (1d array): Voltage slice sorted monotonically increasing.
            islice (1d array): Associated current slice.

        Returns:
            _type_: _description_
        """
        ivgrad = np.gradient(islice, vslice)
        arg_vp = np.argmax(ivgrad)
        # if arg_vf > arg_vp:
        #     print(f"arg_vf = {arg_vf}")
        #     print(f"arg_vp = {arg_vp}")
        #     print(f"i = {i}, j = {j}")
        #     self._plot_sweep_debug(vslice, islice, ivgrad, arg_vf, arg_vp)
        #     raise Exception("Algorithm thinks plasma potential is less then vf.")
        vp = vslice[arg_vp]
        a1, b1 = curve_fit(  # trying gaussian fit of vp to compare to max
            supergaussian,
            vslice,
            ivgrad,
            p0=[ivgrad[arg_vp], vp, 1],
        )

        a2, b2 = curve_fit(
            lambda t, amp, mean, temp: amp * np.exp(-abs(t - mean) / (2 * temp)),
            vslice,
            ivgrad,
            p0=[ivgrad[arg_vp], vp, 1],
        )
        return vp, a2[2]

    # def _find_te(self, vslice, islice, iline, arg_vf, arg_vp):
    #     """Find Te by exponential fit to IV curve.

    #     Args:
    #         vslice (_type_): _description_
    #         islice (bool): _description_
    #     """
    #     ecurr = islice - iline
    #     a, b = curve_fit(
    #         lambda t, amp, temp, offset: amp * np.exp((t + offset) / temp),
    #         vslice[arg_vf // 2 : (arg_vp + arg_vf) // 2],
    #         ecurr[arg_vf // 2 : (arg_vp + arg_vf) // 2],
    #         p0=[1, 1, 0],
    #     )
    #     return a[1]

    # def calculate_density(self, isat, temp, probe_area):
    # def calibrate_probe_area_if(self, inter_trace):
    # def calibrate_probe_area_ts(self, ts_density, ts_time):

    def analyze_sweeps(self, v_slices, i_slices, t_pre_ind):
        """
        Analyze sweeps and store the plasma parameters in self.plasma_params.
        plasma_params index corresponds to the following:
        0: Ion saturation current
        1: Floating potential
        2: Plasma potential
        3: Electron temperature

        TODO: Add error handling for bad fits

        Args:
            v_slices (_type_): _description_
            i_slices (_type_): _description_

        Returns:
        """
        v_units = u.V
        i_units = u.A
        if type(v_slices) == Quantity:
            v_units = v_slices.unit
            v_slices = v_slices.value
        if type(i_slices) == Quantity:
            i_units = i_slices.unit
            i_slices = i_slices.value
        pp_units = np.array([i_units, v_units, v_units, v_units])
        plasma_params = np.empty((v_slices.shape[0], v_slices.shape[1], 4))
        for i in range(v_slices.shape[0]):
            for j in range(v_slices.shape[1]):
                (isat, vf, isat_std, arg_vf) = self._find_isat_vf(
                    v_slices[i, j], i_slices[i, j], t_pre_ind
                )
                plasma_params[i, j, 0] = isat * -1
                plasma_params[i, j, 1] = vf
                plasma_params[i, j, 2], te1 = self._find_plasma_potential(
                    v_slices[i, j, 2 * t_pre_ind :], i_slices[i, j, 2 * t_pre_ind :]
                )
                # te2 = self._find_te(
                #     v_slices[i, j], i_slices[i, j], iline, arg_vf, arg_vp
                # )
                plasma_params[i, j, 3] = te1
                # if abs(te1 - te2) / te1 > 0.1:
                #     print(
                #         f"greater then 10% difference in Te methods, TeG: {te1}, TeE: {te2}, relative error = {abs(te1 - te2) / te1}, shot: {i}, sweep: {j}"
                #     )
                # if isat_std / plasma_params[i, j, 0] > 0.1:
                #     print(
                #         f"greater then 10% std in isat, isat: {plasma_params[i, j, 0]}, std: {isat_std}, relative error = {isat_std / plasma_params[i, j, 0]}, shot: {i}, sweep: {j}"
                #     )
                # if np.any(np.sqrt(np.diag(isat_pcov)) > isat_std):
                #     print(
                #         f"isat_pcov error greater then isat_std, isat: {plasma_params[i, j, 0]}, std: {isat_std}, pcov: {np.sqrt(np.diag(isat_pcov))}, shot: {i}, sweep: {j}"
                #     )
        return plasma_params, pp_units

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
        isat = plasma_params[:, :, 0].to(u.A)
        temp = (plasma_params[:, :, 3] * e / k_B).to(u.J)
        m_i = (ion_mass_factor * m_p).to(u.kg)
        area = probe_area.to(u.m**2)

        electron_density = (isat) / (e * area) * np.sqrt(m_i / temp)
        return electron_density.to(u.cm**-3)
