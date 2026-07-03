"""PI-requested comparison: re-evaluate every accepted point of the V5
micro-MCMC (run with aggressive_optimization: false) using
aggressive_optimization: true, and report the delta-chi2 distribution.
|delta chi2| << 1 across the posterior = chains statistically
indistinguishable between the modes."""
import sys
import numpy as np
import yaml

SCRATCH = ("/private/tmp/claude-501/-Users-vivianmiranda-data-research-WayneHu"
           "-rayne/4e32f7cb-e470-4d65-8ad6-bf48eb7553b7/scratchpad")

with open(f"{SCRATCH}/v5_chains/v5.updated.yaml") as f:
    info = yaml.safe_load(f)
info.pop("sampler", None)
info.pop("output", None)
blk = info["theory"]["axionhmcode_boost.AxionHMcodeBoost"]
blk.pop("processes", None)   # written before the option was removed
blk["legacy_root_finder"] = False
info["debug"] = False

from cobaya.model import get_model
model = get_model(info)
sampled = list(model.parameterization.sampled_params())

with open(f"{SCRATCH}/v5_chains/v5.1.txt") as f:
    header = f.readline().lstrip("#").split()
data = np.loadtxt(f"{SCRATCH}/v5_chains/v5.1.txt")
cols = {name: i for i, name in enumerate(header)}
print(f"chain rows: {len(data)}; sampled params: {sampled}", flush=True)

dchi2 = []
for irow, row in enumerate(data):
    point = {p: row[cols[p]] for p in sampled}
    stored = -row[cols["minuslogpost"]]
    lp = model.logposterior(point, cached=False)
    new = lp.logpost
    dchi2.append(2.0 * (stored - new))
    if irow % 20 == 0:
        print(f"  row {irow}: stored logpost {stored:.4f} "
              f"aggressive {new:.4f} dchi2 {dchi2[-1]:+.2e}", flush=True)

dchi2 = np.array(dchi2)
print("\n==== V5 aggressive-mode postprocess ====")
print(f"points: {len(dchi2)}")
print(f"max |dchi2| = {np.abs(dchi2).max():.3e}")
print(f"mean dchi2  = {dchi2.mean():+.3e},  rms = {dchi2.std():.3e}")
np.save(f"{SCRATCH}/v5_dchi2.npy", dchi2)
