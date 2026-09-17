# Semi-competing risks simulation

R code for a shared Gamma-frailty illness-death model with an interval-censored
non-terminal event and an exactly observed or right-censored terminal event.
Each simulation replication fits both Weibull and B-spline baseline hazards to
the same generated dataset.

This repository contains **simulation code only**. It generates synthetic data;
no clinical dataset, patient records, clinical analysis scripts, manuscript,
or saved study results are included. This is research code, not an R package.

## Requirements

- R 4.4.2 was used for the study.
- R packages `survival` (study version 3.7-0) and `splines` (4.4.2, included in R).
- No ARC account or RStudio installation is required.

Install `survival` if it is missing: `install.packages("survival")`.
For close numerical reproduction, use the study's R/package versions. Different
BLAS/LAPACK libraries and operating systems can produce small differences in
estimates and borderline numerical eligibility; exact binary equality is not promised.

Download or clone this repository. Run the following commands from its root,
the folder containing `reproduce.R` and `analysis0914/`.

## Quick start

Check the numerical core and synthetic observation mechanism (no model fitting):

```sh
Rscript --vanilla analysis0914/tests/test_core.R
Rscript --vanilla analysis0914/tests/test_generation.R
```

Run one full-size replication for the reference scenario and summarize it:

```sh
Rscript --vanilla reproduce.R --smoke
```

This fits both models to 5,000 synthetic subjects. It can take several minutes
or longer, depending on the computer. A one-replication check demonstrates
execution, not estimator accuracy or coverage.

## Complete study

```sh
Rscript --vanilla reproduce.R --full
```

This explicitly starts all seven scenarios, with 5,000 subjects and 1,000
replications per scenario. The convenience launcher runs **sequentially**;
on a personal computer, a full study can take days or weeks. The published-study
computation used parallel cluster jobs; its elapsed time is not a laptop estimate.

| Scenario | Frailty variance | Basis counts (transitions 1, 2, 3) | Transition-1 age coefficient |
|---|---:|---|---:|
| knot02 | 0.25 | (4, 3, 6) | 0.05 |
| knot12 | 0.25 | (4, 4, 6) | 0.05 |
| knot22 | 0.25 | (4, 5, 6) | 0.05 |
| knot23 | 0.25 | (4, 6, 6) | 0.05 |
| theta0001 | 0.0001 | (4, 5, 6) | 0.05 |
| theta05 | 0.5 | (4, 5, 6) | 0.05 |
| beta1age008 | 0.25 | (4, 5, 6) | 0.08 |

The four knot scenarios share generated datasets and event/visit seeds. Their
Weibull fits form one paired reference experiment, not four independent studies.
Additional historical scenario definitions in the core are not part of this
seven-scenario study and are not invoked by the launcher.

## Individual scenarios and parallel chunks

For example, run replications 1-10 of the 1,000-replication reference scenario:

```sh
Rscript --vanilla analysis0914/run_simulation.R run knot22 1 10 --n=5000 --n.rep=1000 --out=results/production
```

For chunk j in 1,...,100, use START=(j-1)*10+1 and END=j*10. Keep code, settings
and output root identical across chunks. Never launch overlapping ranges into
the same output directory. After all expected replications finish, combine:

```sh
Rscript --vanilla analysis0914/run_simulation.R combine knot22 --n=5000 --n.rep=1000 --out=results/production
```

For parallel execution, set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and
MKL_NUM_THREADS to 1 before starting R. The convenience launcher already does so.
Choose cluster resources from measured runtimes rather than assuming one chunk
will finish within a fixed time.

## Observation mechanism and seeds

The event-generation seed is 1000+r and the independent visit seed is 914000+r
for replication r. Visit gaps are exponential with mean seven days. The complete
visit schedule is generated independently of the event times, through the first
visit beyond administrative censoring. Only visits completed by Y2=min(T2,C)
contribute to the observed interval.

For detected illness, the interval uses the preceding negative visit (or zero)
and first positive visit, with R<=Y2. Without a positive visit, the working
observation is L=Y2, R=Inf and delta1=0. No examination is added at Y2.
The working Case 3 interpretation T1=Inf is distinct from latent event truth,
which is retained for missed-illness diagnostics. The working likelihood does
not integrate every undetected latent illness path under discrete inspection.

## Outputs and interpretation

Each scenario contains saved replication records, a combined RDS object,
CSV parameter summaries and numerical/observation diagnostics, plus hazard PDFs.
The launcher writes smoke and full-study results into separate directories.
Numerical eligibility, interval availability and computational completion are
different quantities. Conditional coverage uses eligible intervals; it should
be interpreted with the reported denominators. Numerical success does not
establish nominal coverage, especially for frailty or stronger covariate effects.

Existing replicate files are immutable and manifest-checked. For an interrupted
run, reissue the original range with identical code/settings to verify and skip
saved records. Inspect stale locks before resuming. Modified code or settings
require a new output directory. Do not overwrite saved results to hide failures.
The convenience launcher refuses an existing output directory; use the underlying
entry point for an inspected resume. Combine only after checking completeness.

## Source map

- `reproduce.R`: optional smoke/full-study launcher; no fitting unless requested.
- `analysis0914/run_simulation.R`: run and combine entry point.
- `analysis0914/generate_dataset.R`: generation-only diagnostics.
- `analysis0914/R/core.R`: scenarios, event generator, likelihood and fitting.
- `analysis0914/R/visits.R`: visit schedules and interval construction.
- `analysis0914/R/gradient.R`: derivative support.
- `analysis0914/R/io.R`: manifests, saved records, summaries and figures.
- `analysis0914/tests/`: generation and numerical checks.

The numerical source is preserved from the completed study release; the folder
name `analysis0914` is retained for code compatibility. SHA256SUMS.txt identifies
the supplied scientific source files. The optional comparison with an older
release in the generation test is skipped in this standalone repository.

For manuscript reproducibility, record the Git commit used alongside the run's
R session information. Do not treat a later modified checkout as the same release.
