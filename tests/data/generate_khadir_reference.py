"""Independent NumPy/SciPy translation of the authors' supplementary code.

Run explicitly to regenerate the frozen Fig. 3 profiles; tests never run this.
Source: https://doi.org/10.6084/m9.figshare.7704494 (QuantimageDipole.zip).
The archive omits indexRead/finter. Linearly interpolate the supplied nAu.txt
permittivity versus photon energy, as epsmRead.m specifies, then take sqrt.
Units below are µm. No vectorial_refocusing implementation is imported.
"""

from pathlib import Path

import numpy as np
from scipy.special import spherical_jn, spherical_yn


def generate():
    wavelength, n, m, na, radius = 0.590, 1.5, -100.0, 0.9, 0.037
    k0 = 2 * np.pi / wavelength
    # The two nAu.txt rows bracketing 590 nm: energy (eV), n, k.
    energy = np.array([2.01, 2.13])
    index = np.array([0.21 + 3.272j, 0.29 + 2.863j])
    photon_energy = 1239.8419 / (wavelength * 1000)
    epsilon = np.interp(photon_energy, energy, (index**2).real) + 1j * np.interp(photon_energy, energy, (index**2).imag)
    n_particle = np.sqrt(epsilon)
    relative_index = n_particle / n
    v = n * k0 * radius
    w = relative_index * v
    psi_v = v * spherical_jn(1, v)
    psi_w = w * spherical_jn(1, w)
    xi_v = v * (spherical_jn(1, v) + 1j * spherical_yn(1, v))
    dpv = spherical_jn(1, v) + v * spherical_jn(1, v, derivative=True)
    dpw = spherical_jn(1, w) + w * spherical_jn(1, w, derivative=True)
    dxv = dpv + 1j * (spherical_yn(1, v) + v * spherical_yn(1, v, derivative=True))
    a1 = (relative_index * psi_w * dpv - psi_v * dpw) / (relative_index * psi_w * dxv - xi_v * dpw)
    alpha = 1j * 6 * np.pi * a1 / (n * k0**3)

    size = 801
    dk = n * k0 / 100
    k = (np.arange(size) - size // 2) * dk
    kx, ky = np.meshgrid(k, k)
    pupil = kx**2 + ky**2 <= (na * k0)**2
    kx, ky = np.where(pupil, kx, 0), np.where(pupil, ky, 0)
    gamma = np.sqrt((n * k0)**2 - kx**2 - ky**2)
    gamma_img = np.sqrt(k0**2 - (kx**2 + ky**2) / m**2)
    ko = np.stack((kx, ky, gamma)) / (n * k0)
    ki = np.stack((kx / m, ky / m, gamma_img)) / k0
    kr = np.hypot(kx, ky)
    denominator = np.where(kr == 0, 1, kr)
    s = np.stack((-ky / denominator, np.where(kr == 0, 1, kx / denominator), np.zeros_like(kr)))
    p = np.cross(ko, s, axis=0)
    pp = np.cross(ki, s, axis=0)
    q = np.array([k0**2 * alpha, 0, 0])[:, None, None]
    emitted = 1j / (8 * np.pi**2) / gamma * (q - np.sum(ko * q, axis=0) * ko)
    rotated = s * np.sum(emitted * s, axis=0) + pp * np.sum(emitted * p, axis=0)
    weighted = np.where(pupil, np.sqrt(gamma / gamma_img) * rotated, 0)
    z = np.array([-1.4, -0.98, -0.55, -0.02, 0.55, 0.94, 1.4])
    transmittance, opd_nm = [], []
    crop = slice(size // 2 - 100, size // 2 + 101)
    for distance in z:
        # Put the source at zQ=0 and move the objective to zo=-distance.
        field = dk**2 / m * np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(weighted * np.exp(-1j * gamma * distance), axes=(-2, -1)), axes=(-2, -1)), axes=(-2, -1))
        incident = np.sqrt(n) / m * np.exp(-1j * n * k0 * distance)
        normalized = field[:, size // 2, crop] / incident
        normalized[0] += 1
        transmittance.append(np.sum(abs(normalized)**2, axis=0))
        opd_nm.append(np.angle(normalized[0]) / k0 * 1000)
    spacing_object = 2 * np.pi / (size * dk)
    np.savez(Path(__file__).with_name('khadir_fig3.npz'),
             x_object=np.arange(-100, 101) * spacing_object, z=z,
             transmittance=transmittance, opd_nm=opd_nm,
             n_particle=n_particle, alpha=alpha, dk=dk)
    print('n_particle:', n_particle, 'alpha (µm³):', alpha)
    print('central T:', np.asarray(transmittance)[:, 100])
    print('central OPD (nm):', np.asarray(opd_nm)[:, 100])


if __name__ == '__main__':
    generate()
