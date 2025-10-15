"""
Homogenization parameters are computed for different frequencies using the model
in `./sif/6801.sif`. This will produce `./dat` data folder where all the parameter
data is located. This script converts those parameters in the form that is
needed by the transient homogenization model.

Author: Eelis Takala (ERS)
Email: eelis.takala@gmail.com
Original Date: September 2025
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

def read_elmer_data():
    import pandas as pd
    import re

    names_file = '6801/dat/6801.dat.names'
    names = []
    with open(names_file) as f:
        for line in f:
            m = re.match(r"\s*\d+:\s*(.*)", line)
            if m:
                var = m.group(1).strip().replace(":", "").replace(" ", "_")
                names.append(var)

    data_file = "6801/dat/6801.dat"
    df = pd.read_csv(data_file, delim_whitespace=True, names=names, engine="python")
    df.to_csv("6801.csv", index=False)
    return df

def read_constitutive_values():
    df=read_elmer_data()
    df['nu11']=df['res_nu_11_component(5)_re']+1.j*df['res_nu_11_component(5)_im']
    df['nu22']=df['res_nu_22_component(5)_re']+1.j*df['res_nu_22_component(5)_im']
    df['sigma33']=df['res_sigma_33_component(5)_re']+1.j*df['res_sigma_33_component(5)_im']
    return df

def fit_response(omega, y_meas, n=3, sign_omega=+1, reg_lambda=1e-8, verbose=1):
    """
    Fit a response of the form:
        y(omega) ≈ y0 + alpha * e1^T (I + j*sign_omega*(omega/omega_scale)*S_raw)^(-1) e1

    Parameters
    ----------
    omega : array (rad/s)
        Frequency array in rad/s
    y_meas : array (complex)
        Measured response (sigma, nu, etc.)
    n : int
        Matrix size for the model (default 3)
    sign_omega : +1 or -1
        Convention for Fourier sign (+jω or -jω)
    reg_lambda : float
        Regularization weight on S_raw
    verbose : int
        0 = silent, 1 = summary, 2 = detailed

    Returns
    -------
    results : dict
        {
          "params": fitted parameter vector,
          "S_raw": fitted scaled matrix,
          "Sigma": physical Sigma matrix (=S_raw/omega_scale),
          "alpha": fitted amplitude,
          "y0": fitted offset,
          "y_fit": fitted response (same units as y_meas),
          "residual_norm": final residual norm
        }
    """

    # scale omega and measurement for conditioning
    omega_scale = np.max(np.abs(omega)) if np.max(np.abs(omega)) > 0 else 1.0
    omega_s = omega / omega_scale

    y_scale = np.max(np.abs(y_meas))
    if y_scale == 0:
        y_scale = 1.0
    y_meas_norm = y_meas / y_scale

    # e1 vector
    e1 = np.zeros((n,), dtype=complex)
    e1[0] = 1.0 + 0j

    # build symmetric tridiagonal S from params
    def build_S_from_params(params):
        diag = params[:n]
        off  = params[n:n+(n-1)] if n > 1 else np.array([])
        S = np.zeros((n,n), dtype=float)
        for i in range(n):
            S[i,i] = diag[i]
            if i < n-1:
                S[i,i+1] = off[i]
                S[i+1,i] = off[i]
        return S

    def y_model_from_params(params):
        diag_off_len = n + max(0, n-1)
        S_raw = build_S_from_params(params[:diag_off_len])
        alpha = params[diag_off_len]
        y0    = params[diag_off_len+1]
        I = np.eye(n, dtype=complex)
        out = np.zeros_like(y_meas_norm, dtype=complex)
        for idx, ws in enumerate(omega_s):
            k = sign_omega * 1j * ws
            try:
                invM = np.linalg.inv(I + k * S_raw)
                out[idx] = y0 + alpha * (e1.conj().T @ invM @ e1)
            except np.linalg.LinAlgError:
                out[idx] = 1e12 + 0j
        return out

    def residuals(params):
        y_model = y_model_from_params(params)
        diff = y_model - y_meas_norm
        r = np.hstack([np.real(diff), np.imag(diff)])
        # Tikhonov regularization on S_raw entries
        if reg_lambda > 0:
            reg = np.sqrt(reg_lambda) * params[:(n + max(0, n-1))]
            r = np.concatenate([r, reg])
        return r

    # -----------------------------
    # Initial guess
    # -----------------------------
    diag0 = np.linspace(1.0, 0.2, n)
    off0  = 0.1 * np.ones(max(0, n-1))
    alpha0 = (np.max(np.real(y_meas_norm)) - np.min(np.real(y_meas_norm))) or 1.0
    y0_init = np.real(y_meas_norm[0])
    params0 = np.concatenate([diag0, off0, [alpha0, y0_init]])

    # bounds: diag >=0, off free, alpha free, y0 free
    lb_diag = np.zeros(n)
    ub_diag = np.full(n, 1e6)
    lb_off  = np.full(n-1, -1e6) if n>1 else np.array([])
    ub_off  = np.full(n-1,  1e6) if n>1 else np.array([])
    lb = np.concatenate([lb_diag, lb_off, [-np.inf, -np.inf]])
    ub = np.concatenate([ub_diag, ub_off, [ np.inf,  np.inf]])

    # -----------------------------
    # Run solver
    # -----------------------------
    res = least_squares(residuals, params0, bounds=(lb,ub),
                        verbose=verbose, xtol=1e-12, ftol=1e-12,
                        gtol=1e-12, max_nfev=5000)

    params_fit = res.x
    S_raw_fit = build_S_from_params(params_fit[:(n + max(0, n-1))])
    alpha_fit = params_fit[-2]
    y0_fit    = params_fit[-1]
    Sigma_fit = S_raw_fit / omega_scale

    y_fit_norm = y_model_from_params(params_fit)
    y_fit = y_fit_norm * y_scale

    residual_norm = np.linalg.norm(np.hstack([np.real(y_fit_norm - y_meas_norm),
                                              np.imag(y_fit_norm - y_meas_norm)]))

    results = {
        "params": params_fit,
        "S_raw": S_raw_fit,
        "Sigma": Sigma_fit,
        "alpha": alpha_fit,
        "y0": y0_fit,
        "y_fit": y_fit,
        "residual_norm": residual_norm,
        "omega_scale": omega_scale,
        "y_scale": y_scale,
    }

    return results


if __name__ == "__main__":
    df = read_constitutive_values()
    freq = df['res_angular_frequency'].values/2/np.pi

    omega = 2*np.pi*freq
    for param in ['nu11','nu22','sigma33']:
        y_meas = df[param].values

        results = fit_response(omega, y_meas, n=3, sign_omega=+1)

        print("\nFit results:")
        print("alpha =", results["alpha"])
        print("y0    =", results["y0"])
        print("S_raw =\n", results["S_raw"])
        print("Sigma =\n", results["Sigma"])
        print("Residual norm =", results["residual_norm"])

        # Plot
        plt.figure(figsize=(11,4))
        plt.subplot(1,2,1)
        plt.semilogx(freq, np.real(y_meas), 'o', label='Re meas')
        plt.semilogx(freq, np.real(results["y_fit"]), '-', label='Re fit')
        plt.xlabel('frequency (Hz)'); plt.ylabel(f'{param} Re(y)'); plt.legend()
        plt.subplot(1,2,2)
        plt.semilogx(freq, np.imag(y_meas), 'o', label='Im meas')
        plt.semilogx(freq, np.imag(results["y_fit"]), '-', label='Im fit')
        plt.xlabel('frequency (Hz)'); plt.ylabel(f'{param} Im(y)'); plt.legend()
        plt.tight_layout()
        plt.show()



