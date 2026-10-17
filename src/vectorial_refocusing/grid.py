"""Centered image coordinates (µm) and object wavevectors (rad/µm)."""

import jax.numpy as jnp


def field_to_spectrum(field, *, dx_object, dy_object):
    """Convert a centered field at its sampling plane to continuous e(k).

    ``field`` has shape (3, ny, nx). Object-side pixel intervals are in µm.
    For refocus_field use dx_image=m*dx_object, dy_image=m*dy_object and
    pass the sampling plane's position as z_field.
    """
    axes = (-2, -1)
    spectrum = jnp.fft.fftshift(
        jnp.fft.fft2(
            jnp.fft.ifftshift(field, axes=axes), axes=axes, norm="backward"
        ),
        axes=axes,
    )
    return spectrum * dx_object * dy_object / (2 * jnp.pi)**2


def image_coordinates(shape, *, dx_image, dy_image):
    """Return 1D ``x, y`` coordinates for an image of shape ``(ny, nx)``."""
    ny, nx = shape
    x = (jnp.arange(nx) - nx // 2) * dx_image
    y = (jnp.arange(ny) - ny // 2) * dy_image
    return x, y


def wavevector_grid(shape, *, dx_image, dy_image, magnification):
    """Return centered object-side ``kx, ky, dkx, dky``.

    ``kx`` and ``ky`` have shape ``(ny, nx)``. Image-side transverse
    wavevectors are ``-kx / magnification, -ky / magnification``.
    """
    ny, nx = shape
    kx = 2 * jnp.pi * magnification * jnp.fft.fftshift(
        jnp.fft.fftfreq(nx, d=dx_image)
    )
    ky = 2 * jnp.pi * magnification * jnp.fft.fftshift(
        jnp.fft.fftfreq(ny, d=dy_image)
    )
    dkx = 2 * jnp.pi * magnification / (nx * dx_image)
    dky = 2 * jnp.pi * magnification / (ny * dy_image)
    kx, ky = jnp.meshgrid(kx, ky, indexing="xy")
    return kx, ky, dkx, dky
