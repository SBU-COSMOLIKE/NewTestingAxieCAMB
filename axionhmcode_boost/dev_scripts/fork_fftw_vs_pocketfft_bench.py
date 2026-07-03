"""pyfftw (FFTW) vs numpy/scipy (pocketfft) on the round-3 workload:
real 4096-point rfft + irfft pairs, one pair per halo mass.

pyfftw is exercised exactly in the plan-reuse pattern of cosmo2D.c: build
the FFTW plan once on aligned arrays, then call it on new arrays of the
same shape/dtype (the new-array execute interface behind pyfftw.FFTW.__call__).
"""
import time
import numpy as np
import scipy.fft as sfft
import pyfftw

N = 4096
REPS = 4000
rng = np.random.default_rng(7)
x = rng.standard_normal(N)
X = rng.standard_normal((100, N))   # batched shape (one z-node's masses)


def bench(fn, reps=REPS, warm=50):
  for _ in range(warm):
    fn()
  t0 = time.perf_counter()
  for _ in range(reps):
    fn()
  return (time.perf_counter() - t0) / reps * 1e6   # us per call


# ---- pocketfft -------------------------------------------------------------
c_np = np.fft.rfft(x)
t_np_f = bench(lambda: np.fft.rfft(x))
t_np_b = bench(lambda: np.fft.irfft(c_np, n=N))
t_sp_f = bench(lambda: sfft.rfft(x))
t_sp_b = bench(lambda: sfft.irfft(c_np, n=N))

# ---- pyfftw, plan reuse (cosmo2D.c pattern) --------------------------------
for effort, label in (("FFTW_ESTIMATE", "estimate"), ("FFTW_MEASURE", "measure")):
  a_in = pyfftw.empty_aligned(N, dtype="float64")
  a_out = pyfftw.empty_aligned(N // 2 + 1, dtype="complex128")
  b_in = pyfftw.empty_aligned(N // 2 + 1, dtype="complex128")
  b_out = pyfftw.empty_aligned(N, dtype="float64")
  t0 = time.perf_counter()
  plan_f = pyfftw.FFTW(a_in, a_out, flags=(effort,))
  plan_b = pyfftw.FFTW(b_in, b_out, direction="FFTW_BACKWARD", flags=(effort,))
  t_plan = (time.perf_counter() - t0) * 1e3
  # new-array execution: same size/dtype/alignment as the planned arrays
  x_al = pyfftw.empty_aligned(N, dtype="float64"); x_al[:] = x
  c_al = pyfftw.empty_aligned(N // 2 + 1, dtype="complex128"); c_al[:] = c_np
  t_f = bench(lambda: plan_f(x_al))
  t_b = bench(lambda: plan_b(c_al))
  print(f"pyfftw {label:9s}: plan {t_plan:7.1f} ms | "
        f"fwd {t_f:6.2f} us  inv {t_b:6.2f} us  pair {t_f+t_b:6.2f} us")

print(f"numpy  pocketfft :               | "
      f"fwd {t_np_f:6.2f} us  inv {t_np_b:6.2f} us  pair {t_np_f+t_np_b:6.2f} us")
print(f"scipy  pocketfft :               | "
      f"fwd {t_sp_f:6.2f} us  inv {t_sp_b:6.2f} us  pair {t_sp_f+t_sp_b:6.2f} us")

# ---- batched (100, 4096), the per-redshift shape ---------------------------
C = np.fft.rfft(X, axis=-1)
t_np_batch = bench(lambda: np.fft.rfft(X, axis=-1), reps=400)
t_sp_batch = bench(lambda: sfft.rfft(X, axis=-1), reps=400)
Xa = pyfftw.empty_aligned((100, N), dtype="float64"); Xa[:] = X
Ca = pyfftw.empty_aligned((100, N // 2 + 1), dtype="complex128")
plan_batch = pyfftw.FFTW(Xa, Ca, axes=(-1,), flags=("FFTW_MEASURE",))
t_fw_batch = bench(lambda: plan_batch(Xa), reps=400)
print(f"\nbatched (100, {N}) forward: numpy {t_np_batch:7.1f} us | "
      f"scipy {t_sp_batch:7.1f} us | pyfftw(measure) {t_fw_batch:7.1f} us")

# ---- projection onto the pipeline ------------------------------------------
pair_sp = (t_sp_f + t_sp_b) * 1e-6
pair_fw = None  # recompute from the measure plan above (last loop iteration)
