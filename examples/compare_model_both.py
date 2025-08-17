#!/usr/bin/env python3
import os
import sys
import subprocess
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent


try:
    import hvswdpy as hv  
except ModuleNotFoundError:
    sys.path.insert(0, str(ROOT / 'src'))
    import hvswdpy as hv  


def read_model(path):
    with open(path, 'r') as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    n = int(lines[0])
    vp = []
    vs = []
    rho = []
    th = []
    for i in range(1, n + 1):
        t, a, b, d = lines[i].split()
        if i < n:
            th.append(float(t))
        vp.append(float(a))
        vs.append(float(b))
        rho.append(float(d))
    return (np.array(vp, dtype=np.float64),
            np.array(vs, dtype=np.float64),
            np.array(rho, dtype=np.float64),
            np.array(th, dtype=np.float64))


def run_cli(model_file, nf=100, fmin=0.1, fmax=100.0, nmr=3, nml=3, prec=1.0, nks=0):
    # One run for dispersion (-ph) creates Rph.dat/Lph.dat; capture HV to HV.dat
    hv_orig = ROOT / 'bin' / 'hv_orig'
    cmd = [
        str(hv_orig),
        '-f', str(model_file),
        '-fmin', str(fmin),
        '-fmax', str(fmax),
        '-nf', str(nf),
        '-logsam',
        '-nmr', str(nmr),
        '-nml', str(nml),
        '-prec', str(prec),
        '-nks', str(nks),
        '-ph',
        '-hv',
    ]
    with open(SCRIPT_DIR / 'HV.dat', 'w') as hv_out:
        res = subprocess.run(cmd, cwd=str(SCRIPT_DIR), stdout=hv_out, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"CLI failed: {res.stderr}")
    if not (SCRIPT_DIR / 'Rph.dat').exists():
        raise FileNotFoundError('Rph.dat was not created')
    if not (SCRIPT_DIR / 'HV.dat').exists():
        raise FileNotFoundError('HV.dat was not created (stdout redirection)')


def read_cli_hv(path=SCRIPT_DIR / 'HV.dat'):
    arr = np.fromstring(open(path, 'r').read(), sep=' ')
    if arr.size % 2 != 0:
        raise ValueError('Unexpected HV.dat format')
    pairs = arr.reshape(-1, 2)
    return pairs[:, 0], pairs[:, 1]


def read_cli_dispersion(path):
    # Reads Rph.dat or Lph.dat
    with open(path, 'r') as f:
        content = f.read().strip().split()
    nf = int(content[0]); nm = int(content[1])
    sl = np.array(list(map(float, content[2:2 + nf * nm])))
    flags_raw = content[2 + nf * nm: 2 + nf * nm + nf * nm]
    valid = np.array([tok == 'T' for tok in flags_raw])
    sl = sl.reshape(nm, nf)
    valid = valid.reshape(nm, nf)
    return sl, valid


def main():
    
    model_file = SCRIPT_DIR / 'model.txt'
    nf = 100
    fmin, fmax = 0.1, 10.0
    nmr, nml = 20, 20
    prec = 0.0001
    nks = 0

    # Ensure results dir under examples/
    results_dir = SCRIPT_DIR / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)

    # Read model
    vp, vs, rho, th = read_model(model_file)

    
    freq = np.logspace(np.log10(fmin), np.log10(fmax), nf)

    
    run_cli(model_file, nf=nf, fmin=fmin, fmax=fmax, nmr=nmr, nml=nml, prec=prec, nks=nks)
    
    f_cli, hv_cli = read_cli_hv(SCRIPT_DIR / 'HV.dat')
    
    r_sl_cli, r_va_cli = read_cli_dispersion(SCRIPT_DIR / 'Rph.dat')
    l_sl_cli, l_va_cli = (None, None)
    if (SCRIPT_DIR / 'Lph.dat').exists():
        l_sl_cli, l_va_cli = read_cli_dispersion(SCRIPT_DIR / 'Lph.dat')


    hv_py, status_hv = hv.hv(
        frequencies_hz=freq,
        vp=vp, vs=vs, rho=rho, thickness=th,
        n_rayleigh_modes=nmr, n_love_modes=nml, precision_percent=prec, nks=nks)


    disp = hv.dispersion(
        frequencies_hz=freq,
        vp=vp, vs=vs, rho=rho, thickness=th,
        n_rayleigh_modes=nmr, n_love_modes=nml, precision_percent=prec)
    r_sl_py = disp.rayleigh_slowness
    r_va_py = disp.rayleigh_valid
    l_sl_py = disp.love_slowness
    l_va_py = disp.love_valid
    status_dp = disp.status

    if not np.allclose(f_cli, freq):
        hv_py_plot = np.interp(f_cli, freq, hv_py)
        f_plot = f_cli
    else:
        hv_py_plot = hv_py
        f_plot = freq

    # HV plot
    plt.figure(figsize=(9, 6))
    plt.semilogx(f_cli, hv_cli, 'k-', lw=2, label='CLI HV')
    plt.semilogx(f_plot, hv_py_plot, 'r--', lw=2, label='API HV')
    plt.xlim(1, 5)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('H/V amplitude')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.title('H/V Comparison (CLI vs API)')
    plt.tight_layout()
    hv_png = results_dir / 'compare_hv.png'
    plt.savefig(str(hv_png), dpi=150)

 
    fig, ax = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

    colors = ['C0', 'C1', 'C2', 'C3', 'C4']
    nm_ray = r_sl_cli.shape[0]
    for im in range(nm_ray):
        mask = r_va_cli[im]
        vel_cli = 1.0 / r_sl_cli[im, mask]
        ax[0].scatter(f_cli[mask], vel_cli, s=12, color=colors[im % len(colors)], label=f'M{im+1}')
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Frequency (Hz)')
    ax[0].set_ylabel('Phase velocity (m/s)')
    ax[0].set_title('Rayleigh velocity (CLI)')
    ax[0].grid(True, alpha=0.3)
    if nm_ray > 0:
        ax[0].legend(title='Modes', fontsize=8)

    nm_ray_py = r_sl_py.shape[1]
    for im in range(nm_ray_py):
        mask = (r_va_py[:, im] != 0)
        vel_py = 1.0 / r_sl_py[mask, im]
        ax[1].scatter(freq[mask], vel_py, s=12, color=colors[im % len(colors)], label=f'M{im+1}')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Frequency (Hz)')
    ax[1].set_title('Rayleigh velocity (API)')
    ax[1].grid(True, alpha=0.3)
    if nm_ray_py > 0:
        ax[1].legend(title='Modes', fontsize=8)

    plt.tight_layout()
    ray_png = results_dir / 'compare_rayleigh_dispersion.png'
    plt.savefig(str(ray_png), dpi=150)

    if l_sl_cli is not None:
        fig, ax = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)
        nm_lov = l_sl_cli.shape[0]
        for im in range(nm_lov):
            mask = l_va_cli[im]
            vel_cli = 1.0 / l_sl_cli[im, mask]
            ax[0].scatter(f_cli[mask], vel_cli, s=12, color=colors[im % len(colors)], label=f'M{im+1}')
        ax[0].set_xscale('log')
        ax[0].set_xlabel('Frequency (Hz)')
        ax[0].set_ylabel('Phase velocity (m/s)')
        ax[0].set_title('Love velocity (CLI)')
        ax[0].grid(True, alpha=0.3)
        if nm_lov > 0:
            ax[0].legend(title='Modes', fontsize=8)

        nm_lov_py = l_sl_py.shape[1]
        for im in range(nm_lov_py):
            mask = (l_va_py[:, im] != 0)
            vel_py = 1.0 / l_sl_py[mask, im]
            ax[1].scatter(freq[mask], vel_py, s=12, color=colors[im % len(colors)], label=f'M{im+1}')
        ax[1].set_xscale('log')
        ax[1].set_xlabel('Frequency (Hz)')
        ax[1].set_title('Love velocity (API)')
        ax[1].grid(True, alpha=0.3)
        if nm_lov_py > 0:
            ax[1].legend(title='Modes', fontsize=8)

        plt.tight_layout()
        lov_png = results_dir / 'compare_love_dispersion.png'
        plt.savefig(str(lov_png), dpi=150)
    else:
        lov_png = None

    py_valid = int(np.sum(r_va_py != 0))
    cli_valid = int(np.sum(r_va_cli))
    print('Rayleigh valid points: API =', py_valid, ', CLI =', cli_valid)
    print('Saved plots:')
    print(' -', hv_png)
    print(' -', ray_png)
    if lov_png:
        print(' -', lov_png)


if __name__ == '__main__':
    main()


