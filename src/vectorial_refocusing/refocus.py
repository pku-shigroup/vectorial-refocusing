"""Vectorial imaging and image-side defocus using Khadir's angular spectrum."""

import jax.numpy as jnp

from .grid import wavevector_grid
from .rotation import kz, rotate_field


def pupil_function(kx, ky, k0, na):
    """Boolean circular pupil in object-side wavevector coordinates."""
    return kx**2 + ky**2 <= (k0 * na) ** 2


def refocus_field(
    spectrum,
    *,
    n_img,
    n_obj,
    magnification,
    na,
    k0,
    dx_image,
    dy_image,
    z_field=0.0,
    z_focus=0.0,
    z_image=0.0,
):
    """Map a continuous object angular spectrum to a defocused image field.

    Input and output are centered complex arrays of shape ``(3, ny, nx)``.
    The input is the transverse angular spectrum of the field sampled at
    z_field: E(x,y,z_field) = integral spectrum(k) exp(i(kx*x+ky*y)) d²k.
    It is not an unscaled FFT array. z_field and the object focal position
    z_focus share an object-side coordinate system. Their phase conversion
    is exp(i*kz_obj*(z_focus-z_field)). This extrapolates outgoing waves in
    the background medium; it does not solve the field inside a scatterer.

    Lengths are in µm and k0 is the vacuum wavenumber in rad/µm. Positive
    magnification m describes an inverted image: k_image = -k_object / m.
    z_image is measured along +z from the image focal plane; its complete
    propagation phase is retained. The time convention is exp(-i omega t).

    Assumes homogeneous, lossless, nonmagnetic media, positive indices and
    sampling intervals, and 0 < na < min(n_obj, magnification * n_img).
    Only propagating waves transmitted by the pupil contribute.
    """
    kx, ky, dkx, dky = wavevector_grid(
        spectrum.shape[-2:],
        dx_image=dx_image,
        dy_image=dy_image,
        magnification=magnification,
    )
    pupil = pupil_function(kx, ky, k0, na)

    # Mask before square roots: multiplying NaN by zero is NaN.
    kx = jnp.where(pupil, kx, 0.0)
    ky = jnp.where(pupil, ky, 0.0)
    spectrum = jnp.where(pupil, spectrum, 0.0)
    rotated = rotate_field(
        spectrum,
        n_img=n_img,
        n_obj=n_obj,
        magnification=magnification,
        k0=k0,
        kx=kx,
        ky=ky,
    )
    kz_obj = kz(n_obj, k0, kx, ky)
    kz_img = kz(n_img, k0, kx / magnification, ky / magnification)
    amplitude = -jnp.sqrt(kz_obj / kz_img) / magnification
    phase = jnp.exp(1j * (kz_obj * (z_focus - z_field) + kz_img * z_image))
    weighted_spectrum = amplitude * phase * rotated

    axes = (-2, -1)
    field = jnp.fft.fftshift(
        jnp.fft.fft2(
            jnp.fft.ifftshift(weighted_spectrum, axes=axes),
            axes=axes,
            norm="backward",
        ),
        axes=axes,
    )
    return dkx * dky * field
