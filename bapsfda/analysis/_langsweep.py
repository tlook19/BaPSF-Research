import numpy as np
from ._datastructs import CurrentDaqArray, VoltageDaqArray
from ._datafuncs import butt_low, sav_smooth
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
        isweep_filt = isweep.signal_array.value
        vsweep_filt = vsweep.signal_array.value
        # isweep_filt = butt_low(
        #     isweep.signal_array, 5e5, isweep.run_props["sample_freq"].value, order=4
        # )
        # vsweep_filt = butt_low(
        #     vsweep.signal_array, 5e5, vsweep.run_props["sample_freq"].value, order=4
        # )
        isweep_filt = np.moveaxis(
            np.moveaxis(isweep_filt, 1, 0)
            - np.average(isweep_filt[:, -2000:-1], axis=-1),
            0,
            1,
        )
        t_start_ind = int(
            self._sweep_params["t_start"].to(u.s) / isweep.run_props["dt"]
        )
        t_ramp_ind = int(
            self._sweep_params["t_ramp_up"].to(u.s) / isweep.run_props["dt"]
        )
        t_period_ind = int(
            self._sweep_params["t_period"].to(u.s) / isweep.run_props["dt"]
        )
        ramp_times = (
            np.arange(self._sweep_params["nsweeps"]) * self._sweep_params["t_period"]
            + self._sweep_params["t_start"]
        )
        v_slices = np.empty(
            (vsweep_filt.shape[0], self._sweep_params["nsweeps"], t_ramp_ind)
        )
        i_slices = np.empty(
            (isweep_filt.shape[0], self._sweep_params["nsweeps"], t_ramp_ind)
        )
        for i in range(vsweep_filt.shape[0]):
            for j in range(self._sweep_params["nsweeps"]):
                k = t_start_ind + j * t_period_ind
                vsort = vsweep_filt[i, k : k + t_ramp_ind]
                isort = isweep_filt[i, k : k + t_ramp_ind]
                sort_ind = np.argsort(vsort)
                v_slices[i, j] = sav_smooth(
                    np.take_along_axis(vsort, sort_ind, axis=0), 25
                )
                i_slices[i, j] = sav_smooth(
                    np.take_along_axis(isort, sort_ind, axis=0), 25
                )
        return v_slices * vu, i_slices * iu, ramp_times

    def __find_isat_vf(self, vslice, islice):
        arg_vf = np.argmin(np.abs(islice))
        vf = vslice[arg_vf]
        isat = np.average(islice[: arg_vf // 2])
        isat_std = np.std(islice[: arg_vf // 2])
        a, b = curve_fit(
            lambda x, m, b: m * x + b, vslice[: arg_vf // 2], islice[: arg_vf // 2]
        )
        iline = a[0] * vslice + a[1]
        return isat, vf, isat_std, iline, arg_vf, b

    def __find_plasma_potential(self, vslice, islice):
        """Find the plasma potential for a single sweep via the first derivative method. Report Te estimate from super-gaussian fit.

        Args:
            vslice (1d array): Voltage slice sorted monotically increasing.
            islice (1d array): Associated current slice.

        Returns:
            _type_: _description_
        """
        ivgrad = np.gradient(islice, vslice)
        arg_vp = np.argmax(ivgrad)
        vp = vslice[arg_vp]
        a1, b1 = curve_fit(  # trying gaussian fit of vp to compare to max
            lambda t, amp, mean, std: amp * np.exp(-((t - mean) ** 2) / (2 * std**2)),
            vslice,
            ivgrad,
            p0=[ivgrad[vp], vp, 1],
        )
        print(vp, a1[1], vp - a1[1])
        a2, b2 = curve_fit(
            lambda t, amp, mean, temp: amp * np.exp(-abs(t - mean) / (2 * temp)),
            vslice,
            ivgrad,
            p0=[ivgrad[vp], vp, 1],
        )
        return vp, arg_vp, a2[2]

    def __find_te(self, vslice, islice, iline, arg_vf, arg_vp):
        """Find Te by exponential fit to IV curve.

        Args:
            vslice (_type_): _description_
            islice (bool): _description_
        """
        ecurr = islice - iline
        a, b = curve_fit(
            lambda t, amp, temp, offset: amp * np.exp((t + offset) / temp),
            vslice[arg_vf // 2 : (arg_vp + arg_vf) // 2],
            ecurr[arg_vf // 2 : (arg_vp + arg_vf) // 2],
            p0=[1, 1, 0],
        )
        return a[1]

    # def calculate_density(self, isat, temp, probe_area):
    # def calibrate_probe_area_if(self, inter_trace):
    # def calibrate_probe_area_ts(self, ts_density, ts_time):

    def analyze_sweeps(self, v_slices, i_slices):
        """
        Analyze sweeps and store the plamsa parameters in self.plasma_params.
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
                (
                    plasma_params[i, j, 0],
                    plasma_params[i, j, 1],
                    isat_std,
                    iline,
                    arg_vf,
                    isat_pcov,
                ) = self.__find_isat_vf(v_slices[i, j], i_slices[i, j])
                plasma_params[i, j, 2], arg_vp, te1 = self.__find_plasma_potential(
                    v_slices[i, j], i_slices[i, j]
                )
                te2 = self.__find_te(
                    v_slices[i, j], i_slices[i, j], iline, arg_vf, arg_vp
                )
                plasma_params[i, j, 3] = (te1 + te2) / 2
                if abs(te1 - te2) / te1 > 0.1:
                    print(
                        f"greater then 10% difference in Te methods, Te1: {te1}, Te2: {te2}, relative error = {abs(te1 - te2) / te1}, shot: {i}, sweep: {j}"
                    )
                if isat_std / plasma_params[i, j, 0] > 0.1:
                    print(
                        f"greater then 10% std in isat, isat: {plasma_params[i, j, 0]}, std: {isat_std}, relative error = {isat_std / plasma_params[i, j, 0]}, shot: {i}, sweep: {j}"
                    )
                if np.sqrt(np.diag(isat_pcov)) > isat_std:
                    print(
                        f"isat_pcov error greater then isat_std, isat: {plasma_params[i, j, 0]}, std: {isat_std}, pcov: {np.sqrt(np.diag(isat_pcov))}, shot: {i}, sweep: {j}"
                    )
        return plasma_params * pp_units

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
