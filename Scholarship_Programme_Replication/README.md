# A Stylized Scholarship Programme

Simulation design 2.0.0 | 11 November 2026 | Entirely simulated data

This package reproduces the three figures in Section 6. Each programme awards scholarships to the top 10 percent of student applicants. A programme with 1,000 applicants awards 100 scholarships, and a programme with ten applicants awards one. The experiments compare increasing the number of applicants within one programme with observing more independent ten-applicant programmes.

## Run

Tested with Python 3.12.14.

```bash
python -m pip install -r environment/requirements.txt
python run.py
```

Windows users may use `py` instead of `python`. Linux/macOS users may also use `bash run_all.sh`. Paths are resolved from the scripts, so execution from another directory is supported. No real data, accounts, external downloads, R, Stata or LaTeX installation are required after installing the Python dependencies. Simulation runtime is recorded in `output/logs/run_metadata.json`; plotting and reporting add several seconds. Runtime varies by hardware.

`python run.py --simulation-only` regenerates simulation outputs. `python run.py --figures-only` regenerates figures and the guide from matching included tables. A configuration-hash check prevents using stale tables with a changed design. Use the full run after changing configuration.

## Model and interpretation

Each scholarship programme is one cluster. The running variable is a pre-award academic score drawn independently for each applicant from Uniform[0,100]. The realized cutoff is the midpoint between the scores of the lowest-scoring scholarship recipient and the highest-scoring nonrecipient. The population reference cutoff is 90. The outcome is a subsequent academic performance index, observed for every applicant.

With `x = score / 100`, the main model is:

```text
Y(0) = 40 + 30*x + epsilon,   epsilon ~ Normal(0, 9)
Y(1) = Y(0) + 5
```

Noise is independent across students and independent of pre-award scores. Each student's two potential outcomes share the same noise. In the main model, receiving the scholarship raises that student's index by exactly 5 points. One student's award has no effect on other students' outcomes.

## Figures

| File stem | Purpose |
| --- | --- |
| [`fig01_cutoff_and_spacing`](output/figures/fig01_cutoff_and_spacing.pdf) | Show cutoff uncertainty and the score gap shrinking as cluster size grows; scores are prominently labelled 0-100. |
| [`fig02_large_programme_rd`](output/figures/fig02_large_programme_rd.pdf) | Show local linear RD estimates and confidence-interval coverage as the number of applicants per programme increases. |
| [`fig03_small_programmes`](output/figures/fig03_small_programmes.pdf) | Show increasing precision around 7.7273 as more ten-applicant programmes are observed. The scholarship effect is 5 index points. |

Each figure is supplied as a vector PDF, an editable SVG and a 300-dpi PNG. See the [figure guide](output/Scholarship_Programme_Figure_Guide.pdf), [captions](docs/figure_captions.txt) and [LaTeX figure blocks](docs/figure_blocks.tex).

## Experiments

1. Draw full score clusters of sizes 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, with 4,000 repetitions each. Every size is divisible by 10. Exact adjacent-order-statistic moments benchmark the simulated midpoint and spacing.
2. For sizes at least 100, fit separate local linear regressions at the realized cutoff, with a triangular kernel and the unchanged deterministic bandwidth `h = 0.5*m**(-1/4)` in normalized-score units. Normal intervals use the HC0 sandwich. Record window counts and every failed fit; the delivered run has zero failed fits. With a 10 percent selection share, the treatment-side window can reach the upper score boundary at smaller sizes. The observed finite-sample coverage is retained.
3. Fix `m=10`, select one applicant, and average the outcome difference between the first- and second-ranked applicants over G = 20, 50, 100, 200, 500, 1000, 2000 independent programmes. For efficiency, generate the exact adjacent spacing as `Beta(1,10)` and add independent `Normal(0,18)` noise to `5 + 30*gap`. This gives the same distribution for the statistic as sorting full uniform clusters in this linear, independent-noise model. The standard error is the sample SD of differences divided by `sqrt(G)`. The same interval is evaluated against the outcome-difference target and the 5-point causal effect. The small-cluster experiment concerns this unadjusted adjacent-rank comparison, not all possible estimators.

The supplementary nonlinear variant adds `15*(x-0.9)**2` to the untreated mean and has treatment effect `5 + 4*(x-0.9) + 10*(x-0.9)**2`. It uses paired scores/noise and keeps the target at the population reference cutoff equal to 5. It is included in the numerical table, not as an additional figure. The analytic sensitivity table varies the known simulated untreated slope; it is not a proposed empirical bias correction.

## Results from this run

- At 10 applicants, the exact mean adjacent score gap is 9.0909; at 1,000, it is 0.0999. Monte Carlo cutoff RMSE is 9.802 and 0.942 score points, respectively.
- At 1,000 applicants, the RD mean is 4.9908, RMSE is 0.9740, and normal-interval coverage is 94.225%.
- Coverage is 84.250% at 100 applicants and 95.300% at 10,000. Finite-sample undercoverage is reported without screening.
- The small-cluster outcome-difference target is `5 + 30/11 = 7.727272727`. The baseline component is 2.7273 outcome points, or 54.5% of the true effect.
- With 2,000 small programmes, the mean is 7.7268, the SD is 0.1112, and coverage of the outcome-difference target is 94.375%. Coverage of the causal effect is 0.000% in these replications. A simulated zero is not an assertion of exactly zero finite-sample probability.

## Reproducibility

The master seed is 20261111. Independent NumPy PCG64 streams use `SeedSequence([master_seed, experiment_id, design_key])`. The batch size of 64 is part of the configuration because it affects random-number consumption. In the small-programme experiment, designs with more programmes extend the same sequence of independent draws. All declared designs and repetitions are retained. Figure bands describe variation across simulations.

| Input or script | Generated outputs |
| --- | --- |
| `config/design.json` + `code/03_analyze/simulate.py` | Per-replication statistics, illustrative applicant records, CSV summaries, validation and run metadata |
| `code/05_figures/make_figures.py` + CSV summaries | Three figures, each in PDF/PNG/SVG |
| `code/04_tables/make_report.py` + config/results | This README, figure guide, captions, Section 6 example text, documentation and LaTeX tables |

Full score/noise arrays are regenerated from seeds; they are not stored for every repetition. Files in `data/raw/` preserve simulation-level statistics and one full applicant example for each highlighted cluster size. Never edit these generated records by hand. `docs/variable-dictionary.md` defines their fields.

Run `python code/06_diagnostics/check_reproduction.py` to remove generated inputs in a temporary source copy, rerun from an external working directory and compare numerical and visual outputs byte for byte in the same installed environment. The delivered check is recorded in `output/logs/clean_rerun_check.json`. `SHA256SUMS.txt` is a delivery snapshot. Exact file hashes on another machine additionally depend on the pinned packages and installed fonts.

