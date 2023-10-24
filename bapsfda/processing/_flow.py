import numpy as np
from astropy import units as u
from astropy.constants import e, m_p


def calc_ExB_velocity(E, B):
    v = ((-E * u.V / u.m) / (B * u.G).to(u.T)).to(u.m / u.s)
    return v


def calc_ion_sound_speed(T, mu):
    cs = np.sqrt((T * u.eV) / (mu * m_p)).to(u.m / u.s)
    return cs
