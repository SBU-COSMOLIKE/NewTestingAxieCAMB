"""FFTLog prototype v3: chase the high-k residual.
Fixed: edge = constant extension + analytic sici tail.
Scan: interpolation order (linear vs cubic), N_FFT (4096/8192), window,
combined extend+half-kink handling. Also locate where dev(k) peaks.
"""
import numpy as np
from scipy.special import loggamma, sici
from scipy.interpolate import CubicSpline
from scipy import integrate

N_DATA = 2000


def mellin_j0(s):
  return np.exp((s - 2.0) * np.log(2.0) + 0.5 * np.log(np.pi)
                + loggamma(s / 2.0) - loggamma((3.0 - s) / 2.0))


def fftlog_j0(r, F, q, window, n_fft):
  n_data = len(r)
  dlnr = np.log(r[-1] / r[0]) / (n_data - 1)
  a = np.zeros(n_fft)
  a[:n_data] = F
  a[n_data:] = F[-1]                      # continuous constant extension
  r_full = r[0] * np.exp(np.arange(n_fft) * dlnr)
  a *= r_full**(-q)
  c = np.fft.rfft(a)
  m = np.arange(len(c))
  eta = 2.0 * np.pi * m / (n_fft * dlnr)
  k = (1.0 / r_full[-1]) * np.exp(np.arange(n_fft) * dlnr)
  b = c * mellin_j0(q + 1j * eta) * np.exp(-1j * eta * np.log(k[0] * r[0]))
  m_cut = int((1.0 - window) * (len(c) - 1))
  idx = m > m_cut
  taper = np.ones(len(c))
  taper[idx] = 0.5 * (1.0 + np.cos(np.pi * (m[idx] - m_cut)
                                   / (len(c) - 1 - m_cut)))
  b *= taper
  I = np.fft.irfft(np.conj(b), n=n_fft) * k**(-q)
  aa = k * r[-1]
  si, ci = sici(aa)
  I -= F[-1] * (np.sin(aa) / aa - ci)
  return k, I


r_vir, r_s, r_c = 1.5, 0.2, 0.02
rho_s, rho_core = 1.0e15, 4.0e16


def rho_profile(r):
  nfw = rho_s / ((r / r_s) * (1 + r / r_s)**2)
  sol = rho_core / (1 + 0.091 * (r / r_c)**2)**8
  return np.where(sol > nfw, sol, nfw)


r = np.geomspace(1e-15, r_vir, N_DATA)
F = rho_profile(r) * r**3

k_targets = np.geomspace(1e-4, 15.0, 250)
kernel = np.sin(np.outer(k_targets, r)) / np.outer(k_targets, r)
ref = integrate.simpson(y=rho_profile(r) * r**2 * kernel, x=r, axis=-1)
norm = ref[0]

print(f"{'interp':>7} {'n_fft':>6} {'q':>5} {'win':>5} {'max dev':>10} "
      f"{'at k':>7} {'dev@k>1':>10}")
for interp in ("lin", "cub"):
  for n_fft in (4096, 8192):
    for q in (0.9, 1.1):
      for window in (0.15, 0.25):
        k_out, I_out = fftlog_j0(r, F, q, window, n_fft)
        if interp == "lin":
          got = np.interp(np.log(k_targets), np.log(k_out), I_out)
        else:
          got = CubicSpline(np.log(k_out), I_out)(np.log(k_targets))
        dev = np.abs(got - ref) / norm
        i = int(np.argmax(dev))
        hi = k_targets * r_vir > 1.0
        print(f"{interp:>7} {n_fft:6d} {q:5.2f} {window:5.2f} "
              f"{dev.max():10.2e} {k_targets[i]:7.3f} {dev[hi].max():10.2e}")
