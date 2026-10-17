"""Rodrigues rotation from object to inverted-image propagation direction.

Here s = (sx, sy, 0) is Khadir's rotation axis u, opposite to the paper's
s polarization. Simultaneously reversing s, p and p' leaves Eq. 3 unchanged.
Inputs must describe propagating waves; refocus_field masks the pupil first.
"""

import jax.numpy as jnp


def kz(n, k0, kx, ky):
    """Forward longitudinal wavevector (rad/µm).

    For evanescent waves, pass a complex index to select the square-root
    branch with nonnegative imaginary part in a lossless medium.
    """
    return jnp.sqrt((n * k0) ** 2 - kx**2 - ky**2)


def rotation_cos_sin(n_img, n_obj, magnification, k0, kx, ky):
    """Appendix B: cos(theta) = k_hat.k_hat', sin(theta) = |k_hat × k_hat'|.

    The image is inverted, so the rotation about (sx, sy, 0) is positive.
    These expressions avoid computing an angle only to take its sin/cos.
    """
    kx_img, ky_img = -kx / magnification, -ky / magnification
    k_hat = jnp.stack((kx, ky, kz(n_obj, k0, kx, ky))) / (n_obj * k0)
    k_hat_img = jnp.stack((
        kx_img, ky_img, kz(n_img, k0, kx_img, ky_img)
    )) / (n_img * k0)
    cosine = jnp.sum(k_hat * k_hat_img, axis=0)
    sine = jnp.linalg.norm(jnp.cross(k_hat, k_hat_img, axis=0), axis=0)
    return cosine, sine


def apply_rotation(ex, ey, ez, kx, ky, cosine, sine):
    """Apply the Rodrigues component expansion without building a matrix."""
    k_parallel = jnp.sqrt(kx**2 + ky**2)
    denominator = jnp.where(k_parallel == 0, 1.0, k_parallel)
    sx = ky / denominator
    sy = -kx / denominator

    # On axis sx = sy = 0, cosine = 1 and sine = 0: identity rotation.
    one_minus_cosine = 1 - cosine
    ex_rot = (
        ex * (cosine + sx**2 * one_minus_cosine)
        + ey * (sx * sy * one_minus_cosine)
        + ez * (sy * sine)
    )
    ey_rot = (
        ex * (sy * sx * one_minus_cosine)
        + ey * (cosine + sy**2 * one_minus_cosine)
        - ez * (sx * sine)
    )
    ez_rot = ex * (-sy * sine) + ey * (sx * sine) + ez * cosine
    return ex_rot, ey_rot, ez_rot


def rotate_field(spectrum, *, n_img, n_obj, magnification, k0, kx, ky):
    """Rotate a transverse angular spectrum of shape ``(3, ny, nx)``."""
    cosine, sine = rotation_cos_sin(n_img, n_obj, magnification, k0, kx, ky)
    return jnp.stack(apply_rotation(*spectrum, kx, ky, cosine, sine))
