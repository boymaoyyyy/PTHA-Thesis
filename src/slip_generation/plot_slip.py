# plot_slip.py
import argparse
from pathlib import Path
import traceback
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def plot_slip(npy_path: Path, out_dir: Path, subfault_size_km: float = 20.0, show: bool = True, log_scale: bool = False):
    slip = np.load(npy_path)

    # slip shape: (n_down, n_along)
    n_down, n_along = slip.shape
    length_km = n_along * subfault_size_km
    width_km = n_down * subfault_size_km

    plt.figure(figsize=(10, 4))
    # Choose normalization: linear or log
    norm = None
    if log_scale:
        # avoid zeros or negative values for LogNorm by using the smallest positive entry
        positive = slip[slip > 0]
        if positive.size == 0:
            raise ValueError("Cannot use log scale: slip array has no positive values")
        vmin = positive.min()
        vmax = slip.max()
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    im = plt.imshow(
        slip,
        extent=[0, length_km, 0, width_km],
        origin='lower',
        cmap='hot_r',
        aspect='auto',
        interpolation='none',
        norm=norm
    )

    plt.xlabel('Along-strike distance (km)', fontsize=12)
    plt.ylabel('Down-dip distance (km)', fontsize=12)
    title = f'Slip Distribution — {npy_path.stem}'
    plt.title(title, fontsize=13, pad=15)

    cbar = plt.colorbar(im, shrink=0.8, pad=0.02)
    cbar.set_label('Slip (m)', fontsize=12)

    # Grid (optional)
    plt.xticks(np.arange(0, length_km + 1, max(10, int(length_km // 6))))
    plt.yticks(np.arange(0, width_km + 1, max(5, int(width_km // 4))))
    plt.grid(color='white', linestyle=':', linewidth=0.5, alpha=0.7)

    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = '_log' if log_scale else ''
    out_path = out_dir / (npy_path.stem + f'{suffix}.png')
    plt.tight_layout()
    # Debug: announce save path and catch exceptions during save
    print(f"[plot_slip] Saving plot to: {out_path}")
    try:
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        print(f"[plot_slip] Save complete: {out_path}")
    except Exception as e:
        print(f"[plot_slip] ERROR saving plot: {e}")
        traceback.print_exc()
        raise
    if show:
        plt.show()
    plt.close()
    return out_path


def main():
    repo_root = Path(__file__).resolve().parents[2]
    default_npy = repo_root / 'output' / 'slip_samples' / 'slip_Mw8.5_sample0.npy'
    default_out = repo_root / 'output' / 'slip_samples'

    p = argparse.ArgumentParser(description='Plot slip field saved as .npy')
    p.add_argument('npy', nargs='?', type=Path, default=default_npy, help='Path to .npy slip file')
    p.add_argument('--out-dir', type=Path, default=default_out, help='Directory to save PNG')
    p.add_argument('--subfault-size', type=float, default=20.0, help='Subfault size in km')
    p.add_argument('--no-show', action='store_true', help='Do not display the figure interactively')
    p.add_argument('--log', action='store_true', help='Use logarithmic color scale (LogNorm)')
    args = p.parse_args()

    if not args.npy.exists():
        raise FileNotFoundError(f"Slip file not found: {args.npy}")
    out_path = plot_slip(
        args.npy,
        args.out_dir,
        subfault_size_km=args.subfault_size,
        show=not args.no_show,
        log_scale=args.log,
    )
    print(f"Saved plot to: {out_path}")


if __name__ == '__main__':
    main()