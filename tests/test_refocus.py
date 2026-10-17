"""Physical checks; run with JAX_ENABLE_X64=true python -m unittest discover -s tests."""

import unittest
from types import SimpleNamespace

import jax.numpy as jnp
import numpy as np

from vectorial_refocusing import field_to_spectrum, image_coordinates, refocus_field, wavevector_grid
from vectorial_refocusing.rotation import rotate_field


PARAMS = dict(
    n_obj=1.5, n_img=1.2, magnification=3.0, na=1.1,
    k0=2 * np.pi / 0.6, dx_image=0.8, dy_image=0.7,
)


def reference_geometry(shape, params):
    """Independent NumPy grid and the paper's s, p, p' polarization basis."""
    ny, nx = shape
    m, k0 = params['magnification'], params['k0']
    dx, dy = params['dx_image'], params['dy_image']
    ix = np.arange(nx) - nx // 2
    iy = np.arange(ny) - ny // 2
    dkx, dky = 2 * np.pi * m / (nx * dx), 2 * np.pi * m / (ny * dy)
    kx, ky = np.meshgrid(ix * dkx, iy * dky)
    pupil = kx**2 + ky**2 <= (k0 * params['na'])**2
    kx_safe, ky_safe = np.where(pupil, kx, 0), np.where(pupil, ky, 0)
    gamma = np.sqrt((params['n_obj'] * k0)**2 - kx_safe**2 - ky_safe**2)
    gamma_img = np.sqrt((params['n_img'] * k0)**2 - (kx_safe**2 + ky_safe**2) / m**2)
    ko = np.stack((kx_safe, ky_safe, gamma)) / (params['n_obj'] * k0)
    ki = np.stack((-kx_safe / m, -ky_safe / m, gamma_img)) / (params['n_img'] * k0)
    kr = np.hypot(kx_safe, ky_safe)
    denominator = np.where(kr == 0, 1, kr)
    s = np.stack((-ky_safe / denominator, np.where(kr == 0, 1, kx_safe / denominator), np.zeros(shape)))
    p = np.cross(ko, s, axis=0)
    pp = np.cross(ki, s, axis=0)
    return SimpleNamespace(
        kx=kx, ky=ky, x=ix * dx, y=iy * dy, dkx=dkx, dky=dky,
        pupil=pupil, gamma=gamma, gamma_img=gamma_img,
        ko=ko, ki=ki, s=s, p=p, pp=pp,
    )


def transverse_spectrum(shape, params):
    geometry = reference_geometry(shape, params)
    pupil, ko = geometry.pupil, geometry.ko
    rng = np.random.default_rng(17)
    field = rng.normal(size=(3, *shape)) + 1j * rng.normal(size=(3, *shape))
    field -= ko * np.sum(ko * field, axis=0)
    return np.where(pupil, field, 0), geometry


class RefocusTests(unittest.TestCase):
    def test_axial_plane_wave(self):
        for shape in ((6, 8), (7, 9)):
            with self.subTest(shape=shape):
                geometry = reference_geometry(shape, PARAMS)
                dkx, dky = geometry.dkx, geometry.dky
                amplitude = np.array([1 + 0.3j, -0.2 + 0.7j, 0])
                spectrum = np.zeros((3, *shape), complex)
                spectrum[:, shape[0] // 2, shape[1] // 2] = amplitude / (dkx * dky)
                actual = np.asarray(refocus_field(jnp.asarray(spectrum), **PARAMS))
                expected = -np.sqrt(PARAMS['n_obj'] / PARAMS['n_img']) / PARAMS['magnification'] * amplitude
                self.assertTrue(np.isfinite(actual).all())
                np.testing.assert_allclose(actual, np.broadcast_to(expected[:, None, None], actual.shape), atol=1e-12)

    def test_rotation_against_eq3(self):
        spectrum, geometry = transverse_spectrum((7, 9), PARAMS)
        actual = np.asarray(rotate_field(
            jnp.asarray(spectrum), n_obj=PARAMS['n_obj'], n_img=PARAMS['n_img'],
            magnification=PARAMS['magnification'], k0=PARAMS['k0'],
            kx=jnp.asarray(np.where(geometry.pupil, geometry.kx, 0)),
            ky=jnp.asarray(np.where(geometry.pupil, geometry.ky, 0)),
        ))
        expected = (
            geometry.s * np.sum(spectrum * geometry.s, axis=0)
            + geometry.pp * np.sum(spectrum * geometry.p, axis=0)
        )
        np.testing.assert_allclose(actual, expected, atol=1e-12)
        np.testing.assert_allclose(np.sum(geometry.ki * actual, axis=0), 0, atol=1e-12)
        np.testing.assert_allclose(np.sum(abs(actual)**2, axis=0), np.sum(abs(spectrum)**2, axis=0), atol=1e-12)

    def test_fft_against_direct_integral(self):
        for shape in ((6, 8), (7, 9), (6, 9)):
            with self.subTest(shape=shape):
                spectrum, geometry = transverse_spectrum(shape, PARAMS)
                z_image = 0.37
                rotated = (
                    geometry.s * np.sum(spectrum * geometry.s, axis=0)
                    + geometry.pp * np.sum(spectrum * geometry.p, axis=0)
                )
                m = PARAMS['magnification']
                weighted = (
                    -np.sqrt(geometry.gamma / geometry.gamma_img) / m
                    * rotated * np.exp(1j * geometry.gamma_img * z_image)
                )
                expected = np.empty_like(spectrum)
                for iy, y_image in enumerate(geometry.y):
                    for ix, x_image in enumerate(geometry.x):
                        kernel = np.exp(-1j * (geometry.kx * x_image + geometry.ky * y_image) / m)
                        expected[:, iy, ix] = (
                            geometry.dkx * geometry.dky
                            * np.sum(weighted * kernel, axis=(-2, -1))
                        )
                actual = refocus_field(jnp.asarray(spectrum), **PARAMS, z_image=z_image)
                np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-11)
                actual_x, actual_y = image_coordinates(shape, dx_image=PARAMS['dx_image'], dy_image=PARAMS['dy_image'])
                np.testing.assert_allclose(actual_x, geometry.x, atol=1e-12)
                np.testing.assert_allclose(actual_y, geometry.y, atol=1e-12)
                grid = wavevector_grid(shape, dx_image=PARAMS['dx_image'], dy_image=PARAMS['dy_image'], magnification=m)
                for actual_grid, expected_grid in zip(grid, (geometry.kx, geometry.ky, geometry.dkx, geometry.dky)):
                    np.testing.assert_allclose(actual_grid, expected_grid, atol=1e-12)

    def test_object_translation_inverts_image(self):
        for shape in ((8, 10), (7, 9)):
            with self.subTest(shape=shape):
                spectrum, geometry = transverse_spectrum(shape, PARAMS)
                kx, ky = geometry.kx, geometry.ky
                # Positive object x shift gives a negative image x shift.
                x0 = 2 * PARAMS['dx_image'] / PARAMS['magnification']
                y0 = -PARAMS['dy_image'] / PARAMS['magnification']
                shifted = spectrum * np.exp(-1j * (kx * x0 + ky * y0))
                original = np.asarray(refocus_field(jnp.asarray(spectrum), **PARAMS))
                actual = refocus_field(jnp.asarray(shifted), **PARAMS)
                expected = np.roll(original, shift=(1, -2), axis=(-2, -1))
                np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-11)

    def test_defocus_and_pupil(self):
        # A bin exactly on the pupil edge; the full grid also contains
        # evanescent object- and image-side wavevectors outside the pupil.
        params = dict(PARAMS, magnification=2.0, n_img=1.0, k0=2 * np.pi,
                      na=1.0, dx_image=0.25, dy_image=0.25)
        shape = (8, 8)
        geometry = reference_geometry(shape, params)
        pupil, gamma_img, p = geometry.pupil, geometry.gamma_img, geometry.p
        iy, ix = 4, 5
        spectrum = np.zeros((3, *shape), complex)
        spectrum[:, iy, ix] = p[:, iy, ix]
        focus = np.asarray(refocus_field(jnp.asarray(spectrum), **params))
        self.assertGreater(np.linalg.norm(focus), 0)
        for z_image in (-0.31, 0.31):
            actual = refocus_field(jnp.asarray(spectrum), **params, z_image=z_image)
            np.testing.assert_allclose(actual, focus * np.exp(1j * gamma_img[iy, ix] * z_image), rtol=1e-12, atol=1e-11)
        contaminated = np.where(pupil, spectrum, 2 + 3j)
        actual = np.asarray(refocus_field(jnp.asarray(contaminated), **params))
        self.assertTrue(np.isfinite(actual).all())
        np.testing.assert_allclose(actual, focus, atol=1e-12)
        empty = refocus_field(jnp.zeros_like(jnp.asarray(spectrum)), **params)
        np.testing.assert_array_equal(empty, np.zeros_like(spectrum))

    def test_exit_plane_reference_is_invariant(self):
        spectrum, geometry = transverse_spectrum((7, 9), PARAMS)
        z_field, z_focus, distance = 1.7, -0.23, 0.6
        original = refocus_field(
            jnp.asarray(spectrum), **PARAMS,
            z_field=z_field, z_focus=z_focus, z_image=0.31,
        )
        # Moving the sampling plane forward advances its outgoing spectrum.
        advanced = spectrum * np.exp(1j * geometry.gamma * distance)
        actual = refocus_field(
            jnp.asarray(advanced), **PARAMS,
            z_field=z_field + distance, z_focus=z_focus, z_image=0.31,
        )
        np.testing.assert_allclose(actual, original, rtol=1e-12, atol=1e-11)
        # A common shift of the object-side coordinate origin changes nothing.
        shifted_origin = refocus_field(
            jnp.asarray(spectrum), **PARAMS,
            z_field=z_field + 3, z_focus=z_focus + 3, z_image=0.31,
        )
        np.testing.assert_allclose(shifted_origin, original, rtol=1e-12, atol=1e-11)

    def test_sampled_field_to_continuous_spectrum(self):
        for shape in ((6, 8), (7, 9)):
            with self.subTest(shape=shape):
                spectrum, g = transverse_spectrum(shape, PARAMS)
                m = PARAMS['magnification']
                ex = np.exp(1j * np.outer(g.x / m, g.kx[0]))
                ey = np.exp(1j * np.outer(g.y / m, g.ky[:, 0]))
                # Synthesize a real-space plane by direct positive-k summation.
                field = g.dkx * g.dky * np.einsum('bq,cqp,ap->cba', ey, spectrum, ex)
                actual = field_to_spectrum(
                    jnp.asarray(field), dx_object=PARAMS['dx_image'] / m,
                    dy_object=PARAMS['dy_image'] / m,
                )
                np.testing.assert_allclose(actual, spectrum, rtol=1e-12, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
