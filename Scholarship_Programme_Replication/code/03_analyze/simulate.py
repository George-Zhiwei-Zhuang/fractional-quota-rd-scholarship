"""Fixed-design Monte Carlo experiments; all data are generated locally.

Three cluster experiments use the manuscript's existing estimators.
The fractional quota is 10 percent throughout.
No numerical selection of bandwidths, seeds, or favorable repetitions occurs.
"""
from pathlib import Path
import hashlib
import json
import math
import platform
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
C = json.loads((ROOT / "config/design.json").read_text())
TAB = ROOT / "output/tables"
RAW = ROOT / "data/raw"
LOG = ROOT / "output/logs"
for p in [TAB, RAW, LOG, ROOT / "data/processed"]:
    p.mkdir(parents=True, exist_ok=True)
ALPHA = C["quota_numerator"] / C["quota_denominator"]
REF = 1 - ALPHA
TAU = C["treatment_effect"]
SIGMA = C["outcome_noise_sd"]
Z = 1.959963984540054


def rng_for(experiment, key=0):
    return np.random.Generator(np.random.PCG64(
        np.random.SeedSequence([C["seed"], experiment, key])))


def write_csv(frame, name):
    frame.to_csv(TAB / name, index=False, float_format="%.12g", lineterminator="\n")


def quota(m):
    assert m % C["quota_denominator"] == 0
    return m * C["quota_numerator"] // C["quota_denominator"]


def side_fit(u, y, side):
    """Batch local linear/linear intercept and HC0 sandwich variance.

    u is normalized by the bandwidth for local linear fits. Each row is an
    independent repetition. The sandwich uses squared WLS weights exactly.
    """
    mask = u > 0 if side == 1 else u < 0
    w = np.maximum(1 - np.abs(u), 0) * mask
    n = (w > 0).sum(axis=1)
    s0 = w.sum(axis=1)
    s1 = (w * u).sum(axis=1)
    s2 = (w * u * u).sum(axis=1)
    t0 = (w * y).sum(axis=1)
    t1 = (w * u * y).sum(axis=1)
    det = s0 * s2 - s1 * s1
    valid = (n >= 3) & (det > 1e-12)
    a = np.divide(s2, det, out=np.full_like(s0, np.nan), where=valid)
    b = np.divide(-s1, det, out=np.full_like(s0, np.nan), where=valid)
    intercept = a * t0 + b * t1
    slope = np.divide(s0 * t1 - s1 * t0, det,
                      out=np.full_like(s0, np.nan), where=valid)
    residual = y - intercept[:, None] - slope[:, None] * u
    influence = w * (a[:, None] + b[:, None] * u)
    variance = ((influence * residual) ** 2).sum(axis=1)
    return intercept, variance, n


def rd_fit(u, y):
    plus, vp, np_ = side_fit(u, y, 1)
    minus, vm, nm = side_fit(u, y, -1)
    return plus - minus, np.sqrt(vp + vm), np_, nm


def summary(est, se, target):
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    x, s = est[ok], se[ok]
    r = len(x)
    assert r > 0
    cover = np.mean(np.abs(x - target) <= Z * s)
    return dict(valid_replications=r, failed_replications=len(est)-r,
                target=target, mc_mean=x.mean(), mc_sd=x.std(ddof=1),
                mean_mcse=x.std(ddof=1)/np.sqrt(r), bias=x.mean()-target,
                rmse=np.sqrt(np.mean((x-target)**2)),
                mean_se=s.mean(), q025=np.quantile(x, .025),
                q975=np.quantile(x, .975), coverage=cover,
                coverage_mcse=np.sqrt(cover*(1-cover)/r))


def exact_cutoff(m):
    """Exact first two moments of the midpoint of adjacent uniform order stats."""
    r = m - quota(m)
    den = (m + 1)**2 * (m + 2)
    var_r = r*(m+1-r)/den
    var_r1 = (r+1)*(m-r)/den
    cov = r*(m-r)/den
    mean = (r+.5)/(m+1)
    variance = (var_r+var_r1+2*cov)/4
    return mean, variance, math.sqrt(variance+(mean-REF)**2)


def simulate_clusters():
    nr, batch = C["replications"], C["batch_size"]
    all_rows, cut_rows, rd_rows = [], [], []
    for m in C["cluster_sizes"]:
        rng = rng_for(10, m)
        k, r = quota(m), m-quota(m)
        active_rd = m >= C["rd_minimum_cluster_size"]
        bandwidth = C["bandwidth_constant"] * m**C["bandwidth_exponent"]
        pieces = []
        for start in range(0, nr, batch):
            count = min(batch, nr-start)
            x = rng.uniform(size=(count, m))
            # Ascending positions r-1 and r are first unselected / last selected.
            partitioned = np.partition(x, (r-1, r), axis=1)
            low, high = partitioned[:, r-1], partitioned[:, r]
            cutoff = (low+high)/2
            selected = x > cutoff[:, None]
            assert np.all(selected.sum(axis=1) == k)
            records = dict(m=np.repeat(m, count), k=np.repeat(k, count),
                           repetition=np.arange(start+1, start+count+1),
                           cutoff=cutoff, gap=high-low)
            if active_rd:
                eps = rng.normal(0, SIGMA, size=x.shape)
                y0 = C["baseline_intercept"] + C["baseline_slope_normalized_score"]*x + eps
                y = y0 + TAU*selected
                u = (x-cutoff[:, None])/bandwidth
                est, se, nplus, nminus = rd_fit(u, y)
                dx = x-REF
                # Paired score/noise draws isolate curvature and effect heterogeneity.
                ynl = y0 + C["nonlinear_baseline_curvature"]*dx**2
                ynl += selected*(TAU + C["nonlinear_effect_slope"]*dx
                                 + C["nonlinear_effect_curvature"]*dx**2)
                est_nl, se_nl, _, _ = rd_fit(u, ynl)
                records.update(estimate=est, se=se, estimate_nonlinear=est_nl,
                               se_nonlinear=se_nl, n_plus=nplus, n_minus=nminus)
            pieces.append(pd.DataFrame(records))
        frame = pd.concat(pieces, ignore_index=True)
        mean, var, rmse = exact_cutoff(m)
        mse = (frame.cutoff.to_numpy()-REF)**2
        cut_rows.append(dict(m=m, k=k, quota_share=ALPHA, replications=nr,
                             cutoff_mean=frame.cutoff.mean(),
                             cutoff_mean_mcse=frame.cutoff.std(ddof=1)/np.sqrt(nr),
                             cutoff_rmse=np.sqrt(mse.mean()),
                             cutoff_rmse_mcse=mse.std(ddof=1)/np.sqrt(nr)/(2*np.sqrt(mse.mean())),
                             cutoff_mean_exact=mean, cutoff_sd_exact=np.sqrt(var),
                             cutoff_rmse_exact=rmse,
                             cutoff_rmse_first_order=np.sqrt(REF*(1-REF)/m),
                             gap_mean=frame.gap.mean(),
                             gap_mean_mcse=frame.gap.std(ddof=1)/np.sqrt(nr),
                             gap_mean_exact=1/(m+1)))
        if active_rd:
            for dgp, est_col, se_col in [("linear", "estimate", "se"),
                                        ("nonlinear", "estimate_nonlinear", "se_nonlinear")]:
                rd_rows.append(dict(m=m, k=k, bandwidth_normalized=bandwidth,
                                    bandwidth_score_points=C["score_scale"]*bandwidth,
                                    dgp=dgp, median_n_plus=frame.n_plus.median(),
                                    median_n_minus=frame.n_minus.median(),
                                    min_n_plus=int(frame.n_plus.min()),
                                    min_n_minus=int(frame.n_minus.min()),
                                    **summary(frame[est_col].to_numpy(), frame[se_col].to_numpy(), TAU)))
        all_rows.append(frame)
        print(f"Clusters m={m}: {nr} repetitions", flush=True)
    draws = pd.concat(all_rows, ignore_index=True)
    draws.to_csv(RAW / "cluster_mc.csv.gz", index=False, float_format="%.12g",
                 compression={"method":"gzip", "mtime":0})
    write_csv(pd.DataFrame(cut_rows), "cutoff_and_spacing.csv")
    write_csv(pd.DataFrame(rd_rows), "large_cluster_rd.csv")


def simulate_small_clusters():
    nr, batch = C["replications"], C["batch_size"]
    m = C["small_cluster_size"]
    counts = C["small_cluster_counts"]
    rng = rng_for(20, m)
    results = {g: [] for g in counts}
    slope = C["baseline_slope_normalized_score"]
    target = TAU + slope/(m+1)
    for start in range(0, nr, batch):
        n = min(batch, nr-start)
        # Exact adjacent spacing law, equivalent to sorting m iid Uniform[0,1].
        gaps = rng.beta(1, m, size=(n, max(counts)))
        noise = rng.normal(0, np.sqrt(2)*SIGMA, size=gaps.shape)
        differences = TAU + slope*gaps + noise
        sums = np.cumsum(differences, axis=1)
        sums_sq = np.cumsum(differences**2, axis=1)
        for g in counts:
            est = sums[:,g-1]/g
            sample_var = (sums_sq[:,g-1]-g*est**2)/(g-1)
            se = np.sqrt(np.maximum(sample_var, 0)/g)
            results[g].append(pd.DataFrame(dict(
                G=g, m=m, k=quota(m), repetition=np.arange(start+1,start+n+1),
                estimate=est, se=se)))
    frames, rows = [], []
    for g in counts:
        f = pd.concat(results[g], ignore_index=True)
        rec = summary(f.estimate.to_numpy(), f.se.to_numpy(), target)
        causal_coverage = np.mean(np.abs(f.estimate-TAU) <= Z*f.se)
        rows.append(dict(G=g, m=m, k=quota(m), quota_share=ALPHA, total_applicants=m*g, causal_effect=TAU,
                         baseline_gap=slope/(m+1),
                         sd_exact=np.sqrt((slope**2*m/((m+1)**2*(m+2))+2*SIGMA**2)/g),
                         causal_coverage=causal_coverage,
                         causal_coverage_mcse=np.sqrt(causal_coverage*(1-causal_coverage)/nr),
                         **rec))
        frames.append(f)
    pd.concat(frames, ignore_index=True).to_csv(RAW/"small_cluster_mc.csv.gz",
        index=False, float_format="%.12g", compression={"method":"gzip", "mtime":0})
    write_csv(pd.DataFrame(rows), "small_cluster_contrast.csv")
    # This table is an analytic sensitivity calculation, not a bias correction proposal.
    sensitivity = []
    for size in C["analytic_sensitivity_sizes"]:
        for b in C["analytic_sensitivity_slopes"]:
            sensitivity.append(dict(m=size, k=quota(size), baseline_slope=b,
                                    causal_effect=TAU, expected_baseline_gap=b/(size+1),
                                    observable_contrast=TAU+b/(size+1)))
    write_csv(pd.DataFrame(sensitivity), "analytic_sensitivity.csv")
    print("Small clusters: all cluster-count designs complete", flush=True)


def illustrative_clusters():
    records=[]
    for m in C["illustrative_cluster_sizes"]:
        rng=rng_for(40,m)
        x=rng.uniform(size=m)
        order=np.argsort(-x)
        rank=np.empty(m,dtype=int); rank[order]=np.arange(1,m+1)
        k=quota(m); cutoff=(x[order[k-1]]+x[order[k]])/2
        d=rank<=k
        eps=rng.normal(0,SIGMA,size=m)
        y0=C["baseline_intercept"]+C["baseline_slope_normalized_score"]*x+eps
        records.append(pd.DataFrame(dict(cluster_size=m,applicant=np.arange(1,m+1),
            score=C["score_scale"]*x,rank=rank,selected=d.astype(int),cutoff=C["score_scale"]*cutoff,
            potential_y0=y0,potential_y1=y0+TAU,observed_y=y0+TAU*d)))
    pd.concat(records,ignore_index=True).to_csv(RAW/"illustrative_applicants.csv",index=False,float_format="%.12g")


def validate():
    checks=[]
    def check(name,passed,detail):
        checks.append(dict(check=name,passed=bool(passed),detail=detail))
        if not passed:
            raise AssertionError(f"{name}: {detail}")
    cut=pd.read_csv(TAB/"cutoff_and_spacing.csv")
    gapz=np.abs((cut.gap_mean-cut.gap_mean_exact)/cut.gap_mean_mcse)
    cz=np.abs((cut.cutoff_mean-cut.cutoff_mean_exact)/cut.cutoff_mean_mcse)
    check("quota proportions",np.allclose(cut.k/cut.m,ALPHA),"Every simulated cluster selects exactly the configured fraction.")
    check("uniform spacing mean",gapz.max()<5,f"Maximum absolute Monte Carlo z-score: {gapz.max():.3f}")
    check("midpoint mean",cz.max()<5,f"Maximum absolute Monte Carlo z-score: {cz.max():.3f}")
    rd=pd.read_csv(TAB/"large_cluster_rd.csv")
    check("all RD fits usable",rd.failed_replications.sum()==0,f"Failed fits: {int(rd.failed_replications.sum())}; none silently dropped.")
    small=pd.read_csv(TAB/"small_cluster_contrast.csv")
    z=np.abs(small.bias/small.mean_mcse)
    check("small contrast expectation",z.max()<5,f"Maximum absolute Monte Carlo z-score: {z.max():.3f}")
    ratios=small.mc_sd/small.sd_exact
    check("small contrast variance",np.all((ratios>.94)&(ratios<1.06)),f"MC/exact SD range: {ratios.min():.4f} to {ratios.max():.4f}")
    # Independent numerical check of the batched intercept / HC0 implementation.
    rng=rng_for(90)
    u=rng.uniform(-1,1,200); y=2+3*u+4*(u>0)+rng.normal(size=200)
    ai,vi,_=side_fit(u[None,:],y[None,:],1)
    mask=(u>0); xx=np.column_stack([np.ones(mask.sum()),u[mask]])
    w=1-u[mask]; yy=y[mask]
    inv=np.linalg.inv(xx.T@(w[:,None]*xx)); beta=inv@(xx.T@(w*yy))
    resid=yy-xx@beta
    vv=inv@(xx.T@(((w*resid)**2)[:,None]*xx))@inv
    check("independent WLS check",np.isclose(ai[0],beta[0]) and np.isclose(vi[0],vv[0,0]),
          "Direct matrix solution and HC0 sandwich agree with the batched calculation.")
    (LOG/"validation.json").write_text(json.dumps(checks,indent=2)+"\n")
    print(f"Validation: {len(checks)} checks passed",flush=True)


def main():
    assert 0 < C["quota_numerator"] < C["quota_denominator"]
    assert C["small_cluster_size"] in C["cluster_sizes"]
    assert all(m % C["quota_denominator"] == 0 for m in C["cluster_sizes"])
    start=time.perf_counter()
    illustrative_clusters()
    simulate_clusters()
    simulate_small_clusters()
    validate()
    versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__)
    metadata=dict(design_version=C["version"],seed=C["seed"],rng="numpy PCG64 / SeedSequence",
                  config_sha256=hashlib.sha256((ROOT/"config/design.json").read_bytes()).hexdigest(),
                  versions=versions,simulation_seconds=round(time.perf_counter()-start,3),
                  data_origin="Entirely synthetic; no downloads or private inputs",
                  selection="All prespecified sizes and repetitions retained; no seed selection")
    (LOG/"run_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")


if __name__=="__main__":
    main()
