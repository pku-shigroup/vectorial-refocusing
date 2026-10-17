"""Refocus a given transverse Gaussian angular spectrum; lengths are in µm."""

import jax.numpy as jnp

from vectorial_refocusing import image_coordinates, refocus_field, wavevector_grid


def main():
    shape = (64, 64)
    sampling = dict(dx_image=0.5, dy_image=0.5, magnification=10.0)
    optics = dict(k0=2 * jnp.pi / 0.6, n_obj=1.5, n_img=1.0, na=1.0)
    kx, ky, _, _ = wavevector_grid(shape, **sampling)

    # Construct one possible input: project x polarization onto the plane
    # perpendicular to each object wavevector. No source model is assumed.
    pupil = kx**2 + ky**2 <= (optics['k0'] * optics['na'])**2
    kx_safe = jnp.where(pupil, kx, 0.0)
    ky_safe = jnp.where(pupil, ky, 0.0)
    k_obj = optics['n_obj'] * optics['k0']
    kz_obj = jnp.sqrt(k_obj**2 - kx_safe**2 - ky_safe**2)
    k_hat = jnp.stack((kx_safe, ky_safe, kz_obj)) / k_obj
    polarization = jnp.array([1.0, 0.0, 0.0])[:, None, None]
    transverse = polarization - k_hat * k_hat[0]
    envelope = jnp.exp(-(kx**2 + ky**2) / (0.4 * optics['k0'])**2)
    spectrum = jnp.where(pupil, transverse * envelope, 0.0).astype(jnp.complex128)

    # Sample the outgoing field beyond its reference plane. The objective
    # may focus behind that plane; refocus_field restores the reference.
    z_field, z_focus = 2.0, 0.3
    spectrum = spectrum * jnp.exp(1j * kz_obj * z_field)

    x, y = image_coordinates(shape, dx_image=sampling['dx_image'], dy_image=sampling['dy_image'])
    print(f'Image grid: {shape}, x = [{float(x[0]):.1f}, {float(x[-1]):.1f}] µm')
    for z_image in (0.0, 20.0):
        field = refocus_field(
            spectrum, **sampling, **optics,
            z_field=z_field, z_focus=z_focus, z_image=z_image,
        )
        field_squared = jnp.sum(jnp.abs(field)**2, axis=0)
        print(f'z_image = {z_image:5.1f} µm: peak |E|² = {float(field_squared.max()):.6g}, '
              f'all finite = {bool(jnp.isfinite(field).all())}')


if __name__ == '__main__':
    main()
