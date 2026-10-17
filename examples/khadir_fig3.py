"""Plot computed Fig. 3 profiles against independent and published references.

Run from the repository: JAX_ENABLE_X64=true uv run --with matplotlib examples/khadir_fig3.py
The optional matplotlib dependency is only needed for this plot.
"""

from pathlib import Path

import jax.numpy as jnp
import numpy as np

from vectorial_refocusing import dipole_spectrum, refocus_field, sphere_polarizability, wavevector_grid


ROOT = Path(__file__).resolve().parents[1]


def simulate():
    reference = np.load(ROOT / 'tests/data/khadir_fig3.npz')
    k0, n, m, size = 2 * np.pi / 0.590, 1.5, 100.0, 801
    spacing = m * 2 * np.pi / (size * float(reference['dk']))
    sampling = dict(dx_image=spacing, dy_image=spacing, magnification=m)
    kx, ky, _, _ = wavevector_grid((size, size), **sampling)
    alpha = sphere_polarizability(radius=0.037, k0=k0, n_medium=n,
                                  n_particle=complex(reference['n_particle']))
    z_field = 2.0
    spectrum = dipole_spectrum(
        jnp.array([k0**2 * alpha, 0, 0]), kx=kx, ky=ky, k0=k0,
        n_source=n, n_obj=n, na=0.9, z_field=z_field,
    )
    profiles = []
    for z in reference['z']:
        z_focus = -float(z)
        field = refocus_field(spectrum, **sampling, k0=k0, n_obj=n, n_img=1.0,
                              na=0.9, z_field=z_field, z_focus=z_focus)
        incident = -np.sqrt(n) / m * np.exp(1j * n * k0 * z_focus)
        total = np.asarray(field[:, 400, 300:501] / incident).copy()
        total[0] += 1
        profiles.append((np.sum(abs(total)**2, axis=0), np.angle(total[0]) / k0 * 1000))
    return reference, np.asarray(profiles)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    reference, profiles = simulate()
    published = np.loadtxt(ROOT / 'tests/data/khadir_fig3_published.csv', delimiter=',')
    fig, axes = plt.subplots(2, 7, figsize=(17, 5), sharex=True, sharey='row')
    for i, z in enumerate(reference['z']):
        for observable, key in enumerate(('transmittance', 'opd_nm')):
            ax = axes[observable, i]
            ax.plot(reference['x_object'], profiles[i, observable], label='This implementation', lw=1.5)
            ax.plot(reference['x_object'], reference[key][i], 'k--', label='Independent reference', lw=0.8)
            rows = published[(published[:, 0] == observable) & np.isclose(published[:, 1], z)]
            ax.plot(rows[::3, 2], rows[::3, 3], '.', color='tab:orange', ms=2.5, label='Published theory (digitized)')
            ax.set_xlim(-3, 3)
            ax.grid(alpha=0.2)
            if observable == 0:
                ax.set_title(f'zQ − z_focus = {z * 1000:+.0f} nm', fontsize=9)
            else:
                ax.set_xlabel('Object coordinate (µm)')
    axes[0, 0].set_ylabel('Transmittance')
    axes[1, 0].set_ylabel('OPD (nm)')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3)
    fig.suptitle('Khadir Fig. 3: 74 nm Au, λ = 590 nm, n = 1.5, m = 100, NA = 0.9')
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    output = ROOT / 'doc/fig3_comparison.png'
    fig.savefig(output, dpi=180)
    print(f'Saved {output}')
    for observable, key in enumerate(('transmittance', 'opd_nm')):
        reference_error = np.max(abs(profiles[:, observable] - reference[key]))
        errors = []
        for i, z in enumerate(reference['z']):
            rows = published[(published[:, 0] == observable) & np.isclose(published[:, 1], z)]
            errors.extend(np.interp(rows[:, 2], reference['x_object'], profiles[i, observable]) - rows[:, 3])
        print(f'{key}: independent max error={reference_error:.6g}, '
              f'published max error={np.max(abs(np.asarray(errors))):.6g}, '
              f'published RMS={np.sqrt(np.mean(np.asarray(errors)**2)):.6g}')


if __name__ == '__main__':
    main()
