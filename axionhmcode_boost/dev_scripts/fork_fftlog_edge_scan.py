"""FFTLog prototype v2: fix the truncation-edge error.

Variant B: half-weight the last data sample (trapezoid-consistent edge).
Variant C: fill the back padding with the constant F(r_vir) (continuous
   extension, no jump anywhere; the r^-q bias kills it long before the
   periodic wrap), then subtract the extension's contribution analytically:
     int_{r_vir}^inf F(r_vir) j0(kr) dr/r
         = F(r_vir) * [ sin(a)/a - Ci(a) ],  a = k r_vir,
   via scipy.special.sici (the same special function upstream already uses
   for the cold NFW k-space profile).
Deviations reported against upstream's rule (dense Simpson on the same
samples), split at k r_vir = 1.
"""
import numpy as np
from scipy.special import loggamma, sici
from scipy import integrate

N_DATA = 2000
N_FFT = 4096


def mellin_j0(s):
  return np.exp((s - 2.0) * np.log(2.0) + 0.5 * np.log(np.pi)
                + loggamma(s / 2.0) - loggamma((3.0 - s) / 2.0))


def fftlog_j0(r, F, q, window, n_fft=N_FFT, edge="none"):
  n_data = len(r)
  dlnr = np.log(r[-1] / r[0]) / (n_data - 1)
  a = np.zeros(n_fft)
  a[:n_data] = F
  if edge == "half":
    a[n_data - 1] *= 0.5
  elif edge == "extend":
    # constant continuation of F beyond r_vir through the padding
    a[n_data:] = F[-1]
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
  if edge == "extend":
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
lo = k_targets * r_vir < 1.0

print(f"{'edge':>8} {'q':>5} {'win':>5} {'max dev, krvir<1':>17} "
      f"{'max dev, krvir>1':>17}")
for edge in ("none", "half", "extend"):
  for q in (0.8, 1.1, 1.5):
    for window in (0.25,):
      k_out, I_out = fftlog_j0(r, F, q, window, edge=edge)
      got = np.interp(np.log(k_targets), np.log(k_out), I_out)
      dev = np.abs(got - ref) / norm
      print(f"{edge:>8} {q:5.2f} {window:5.2f} {dev[lo].max():17.3e} "
            f"{dev[~lo].max():17.3e}")
