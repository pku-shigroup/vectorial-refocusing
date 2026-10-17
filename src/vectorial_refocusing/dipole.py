"""Dipole angular spectrum transmitted through one planar interface."""

import jax.numpy as jnp

from .rotation import kz


def sphere_polarizability(*, radius, k0, n_medium, n_particle):
    """Electric-dipole Mie polarizability in µm³ (Khadir Appendix C).

    The induced source strength is Q = k0² * alpha * E_excitation.
    Uses outgoing spherical Hankel functions for exp(-i omega t).
    """
    m = n_particle / n_medium
    v = n_medium * k0 * radius
    w = m * v
    psi_v = jnp.sin(v) / v - jnp.cos(v)
    psi_w = jnp.sin(w) / w - jnp.cos(w)
    xi_v = psi_v - 1j * (jnp.cos(v) / v + jnp.sin(v))
    psi_v_prime = jnp.sin(v) - psi_v / v
    psi_w_prime = jnp.sin(w) - psi_w / w
    xi_v_prime = jnp.sin(v) - 1j * jnp.cos(v) - xi_v / v
    a1 = (m * psi_w * psi_v_prime - psi_v * psi_w_prime) / (
        m * psi_w * xi_v_prime - xi_v * psi_w_prime
    )
    return 1j * 6 * jnp.pi * a1 / (n_medium * k0**3)


def dipole_spectrum(
    q, *, kx, ky, k0, n_source, n_obj, na, position=(0.0, 0.0, 0.0), z_field=0.0
):
    """Field angular spectrum at z_field, following Khadir Eq. 14–16.

    The interface is z=0; the given dipole is in n_source at z<=0 and the
    collection medium n_obj occupies z>=0. z_field>=0 is in that medium.
    For equal indices the interface disappears and any source position
    below z_field is allowed. Lengths are µm, wavevectors rad/µm.

    q is the paper's source strength Q = k0²*mu/epsilon0, shape (3,).
    Output is a centered (3, ny, nx) continuous angular spectrum, zero
    outside na. Requires na<n_obj. Source-side evanescent waves use the
    outgoing square-root branch; they may transmit into propagating waves.
    The prescribed dipole moment is not corrected for interface feedback.
    """
    pupil = kx**2 + ky**2 <= (k0 * na)**2
    kx, ky = jnp.where(pupil, kx, 0.0), jnp.where(pupil, ky, 0.0)
    gamma_b = kz(n_source + 0j, k0, kx, ky)
    gamma = kz(n_obj, k0, kx, ky)
    kr = jnp.sqrt(kx**2 + ky**2)
    denominator = jnp.where(kr == 0, 1.0, kr)
    s = jnp.stack((
        -ky / denominator,
        jnp.where(kr == 0, 1.0, kx / denominator),
        jnp.zeros_like(kr),
    ))
    kb_hat = jnp.stack((kx, ky, gamma_b)) / (n_source * k0)
    k_hat = jnp.stack((kx, ky, gamma)) / (n_obj * k0)
    pb = jnp.cross(kb_hat, s, axis=0)
    p = jnp.cross(k_hat, s, axis=0)
    q = jnp.asarray(q)[:, None, None]

    # Combine t_s/gamma_b and t_p/gamma_b analytically, including the
    # critical-angle limit gamma_b=0. These are bilinear dot products.
    qs = 2 / (gamma_b + gamma) * jnp.sum(q * s, axis=0) * s
    qp = (
        2 * n_source * n_obj / (n_obj**2 * gamma_b + n_source**2 * gamma)
        * jnp.sum(q * pb, axis=0) * p
    )
    xq, yq, zq = position
    phase = jnp.exp(1j * (-kx * xq - ky * yq - gamma_b * zq + gamma * z_field))
    spectrum = 1j / (8 * jnp.pi**2) * (qs + qp) * phase
    return jnp.where(pupil, spectrum, 0.0)
