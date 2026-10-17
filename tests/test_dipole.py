"""Interface physics and the seven published Fig. 3 theory profiles."""

from pathlib import Path
import unittest

import jax.numpy as jnp
import numpy as np

from vectorial_refocusing import dipole_spectrum, refocus_field, sphere_polarizability, wavevector_grid


DATA = Path(__file__).parent / 'data'


class DipoleTests(unittest.TestCase):
    def test_interface_transmission_and_evanescent_source(self):
        k0, n_source, n_obj = 2 * np.pi, 1.0, 1.5
        kx = np.array([[0.0, 0.4, 1.0, 1.2]]) * k0
        ky = np.zeros_like(kx)
        q = np.array([0.4 + 0.1j, 0.7 - 0.3j, 1.0])
        zq, z_field = -0.07, 0.3
        actual = np.asarray(dipole_spectrum(
            q, kx=jnp.asarray(kx), ky=jnp.asarray(ky), k0=k0,
            n_source=n_source, n_obj=n_obj, na=1.3,
            position=(0.0, 0.0, zq), z_field=z_field,
        ))
        # For ky=0 the paper's s=y, p=(-kz,0,kx)/(n*k0).
        gamma_b = np.sqrt((n_source * k0)**2 - kx**2 + 0j)
        gamma = np.sqrt((n_obj * k0)**2 - kx**2)
        pb = np.stack((-gamma_b, ky, kx)) / (n_source * k0)
        p = np.stack((-gamma, ky, kx)) / (n_obj * k0)
        s = np.zeros_like(p)
        s[1] = 1
        expected = (
            2 / (gamma_b + gamma) * q[1] * s
            + 2 * n_source * n_obj / (n_obj**2 * gamma_b + n_source**2 * gamma)
            * np.sum(q[:, None, None] * pb, axis=0) * p
        ) * (1j / (8 * np.pi**2)) * np.exp(1j * (gamma * z_field - gamma_b * zq))
        np.testing.assert_allclose(actual, expected, atol=1e-12)
        self.assertTrue(np.isfinite(actual).all())  # includes gamma_b=0
        self.assertGreater(np.linalg.norm(actual[:, :, -1]), 0)
        ko = np.stack((kx, ky, gamma)) / (n_obj * k0)
        np.testing.assert_allclose(np.sum(ko * actual, axis=0), 0, atol=1e-12)

    def test_equal_indices_give_homogeneous_dipole(self):
        k0, n = 2 * np.pi / 0.59, 1.5
        kx, ky, _, _ = wavevector_grid((9, 11), dx_image=1.0, dy_image=1.0, magnification=2.0)
        pupil = np.asarray(kx**2 + ky**2 <= (0.9 * k0)**2)
        kx_safe, ky_safe = np.where(pupil, kx, 0), np.where(pupil, ky, 0)
        gamma = np.sqrt((n * k0)**2 - kx_safe**2 - ky_safe**2)
        ko = np.stack((kx_safe, ky_safe, gamma)) / (n * k0)
        q = np.array([1 + 0.2j, -0.3j, 0.8])[:, None, None]
        expected = 1j / (8 * np.pi**2 * gamma) * (q - ko * np.sum(ko * q, axis=0))
        expected *= np.exp(1j * gamma * 0.6)
        actual = dipole_spectrum(q[:, 0, 0], kx=kx, ky=ky, k0=k0,
                                 n_source=n, n_obj=n, na=0.9, z_field=0.6)
        np.testing.assert_allclose(actual, np.where(pupil, expected, 0), atol=1e-12)

    def test_fig3_theory_profiles(self):
        reference = np.load(DATA / 'khadir_fig3.npz')
        published = np.loadtxt(DATA / 'khadir_fig3_published.csv', delimiter=',')
        k0, n, m, size = 2 * np.pi / 0.590, 1.5, 100.0, 801
        spacing = m * 2 * np.pi / (size * float(reference['dk']))
        sampling = dict(dx_image=spacing, dy_image=spacing, magnification=m)
        kx, ky, _, _ = wavevector_grid((size, size), **sampling)
        alpha = sphere_polarizability(radius=0.037, k0=k0, n_medium=n,
                                      n_particle=complex(reference['n_particle']))
        np.testing.assert_allclose(alpha, reference['alpha'], rtol=1e-12)
        # The field is sampled beyond the particle, while the focal plane
        # moves through it. This exercises the actual simulation use case.
        z_field = 2.0
        spectrum = dipole_spectrum(
            jnp.array([k0**2 * alpha, 0, 0]), kx=kx, ky=ky, k0=k0,
            n_source=n, n_obj=n, na=0.9, z_field=z_field,
        )
        crop = slice(size // 2 - 100, size // 2 + 101)
        for index, z in enumerate(reference['z']):
            with self.subTest(z=z):
                z_focus = -float(z)
                field = refocus_field(
                    spectrum, **sampling, k0=k0, n_obj=n, n_img=1.0,
                    na=0.9, z_field=z_field, z_focus=z_focus,
                )
                incident = -np.sqrt(n) / m * np.exp(1j * n * k0 * z_focus)
                total = np.asarray(field[:, size // 2, crop] / incident).copy()
                total[0] += 1
                profiles = (np.sum(abs(total)**2, axis=0), np.angle(total[0]) / k0 * 1000)
                for observable, (profile, key) in enumerate(zip(profiles, ('transmittance', 'opd_nm'))):
                    np.testing.assert_allclose(profile, reference[key][index], rtol=0, atol=1e-11 if observable == 0 else 1e-9)
                    rows = published[(published[:, 0] == observable) & np.isclose(published[:, 1], z)]
                    difference = np.interp(rows[:, 2], reference['x_object'], profile) - rows[:, 3]
                    # Published curves are outlined graphics with manually
                    # placed grid lines; these tolerances include that scale.
                    self.assertLess(np.max(abs(difference)), 0.02 if observable == 0 else 1.0)
                    self.assertLess(np.sqrt(np.mean(difference**2)), 0.007 if observable == 0 else 0.45)


if __name__ == '__main__':
    unittest.main()
