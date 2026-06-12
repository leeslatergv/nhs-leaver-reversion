# NHS trust leaver-rate mean reversion

A reproducible pipeline asking one question: **does a trust's leaver rate
mean-revert?** When a trust has an unusually high-turnover period, does the next
period tend to fall back towards the trust's own usual level, how far, at what
horizon, and along what path?

One of the checks also separates genuine mean reversion from ordinary regression
to the mean (a high reading that was mostly sampling noise, which snaps back for a
different reason), since the two look alike in the raw series. But that is a
robustness step; the headline is the persistence of the rate itself.

It runs at two horizons: a monthly first pass on a short window, then the full
Oct 2009 to Mar 2026 back-series collapsed to financial years.

This repository is the **method and the inputs only**. It deliberately does not
report results, conclusions, or any commentary on individual organisations. Run
it to produce your own figures and tables locally.

## Run it

```bash
./run.sh            # the full pipeline, 01 through 06
```

Or a stage at a time, e.g. `Rscript R/03_annual.R`. R 4.4+; every package is
auto-loaded via `pacman` (tidyverse, fixest, broom, scales, readxl, lubridate,
here), so the first run installs anything missing. Paths resolve from the project
root via the `here` package, so it does not matter where you launch R from.

Result tables and any figures are written to `outputs/` (gitignored, regenerable).

## Layout

```
R/          the pipeline, numbered in run order
data/raw/   the source extracts (public, see Data below)
data/clean/ derived panels (gitignored, rebuilt by 01 and 05)
outputs/    result tables and figures (gitignored, regenerable)
```

## The pipeline

Each script is a stage; comments at the top of each explain its method in detail.

- `R/01_prep.R` - assembles the clean trust x month panel from the monthly
  benchmarking extract, rebuilds the start-of-period denominator, and writes the
  leaver rate three ways (end / start / average stock) for a denominator
  robustness test. -> `data/clean/panel.rds`
- `R/02_analysis.R` - the monthly tests: a within-trust AR(1) with trust + month
  fixed effects, a half-panel jackknife to undo the Nickell small-T bias, a
  reliability split (sampling-noise floor) to separate real reversion from
  regression to the mean, plus denominator robustness and an above/below-own-mean
  asymmetry check.
- `R/03_annual.R` - the same logic at annual horizon on the 2018+ all-groups
  back-series collapsed to financial years: annual persistence, the Nickell
  jackknife, and the reliability split.
- `R/04_frequencies.R` - the term structure: the same within-trust AR(1) run at
  monthly, quarterly, half-yearly and yearly horizons, lined up.
- `R/05_combined.R` - splices the 2009-2018 staff-group file onto the 2018+ series
  (they meet at Jul 2018) for a continuous ~16-year annual series, then re-runs
  the term structure. -> `data/clean/combined_monthly.rds`
- `R/06_dynamics.R` - the shape of the path: an AR(2) within trust and its impulse
  response, to distinguish a monotonic glide back from an overshoot.

## Data

All three raw files are from **NHS Workforce Statistics** (NHS England, formerly
NHS Digital), the "turnover from organisation benchmarking" series. They are
aggregate counts at organisation x staff-group x month level (leavers, joiners,
headcount). No individual records. Published as official statistics under the
[Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

- `data/raw/turnover_org_month.csv` - the monthly extract, headcount and FTE
  (used by 01/02).
- `data/raw/turnover_org_month_full.csv` - the Jul 2018 to Mar 2026 all-staff-group
  series (used by 03/04/05).
- `data/raw/turnover_2009_2018_staffgroup.xlsx` - the supplementary Oct 2009 to
  Jul 2018 file, FTE, that backfills the seven earlier years (used by 05).

> Contains public sector information licensed under the Open Government Licence
> v3.0. This repository redistributes the extracts for reproducibility only;
> NHS England remains the source and authority for the figures.

## Notes on method

- The monthly reliability is computed on the trust-demeaned rate, which still
  holds the seasonal swing; the annual scripts use the two-way-within residual.
- Dynamic-panel persistence is biased downward at short T (Nickell); the jackknife
  and one-offs-removed variants bracket it rather than pin it.
- The NHS leaver definition includes trust-to-trust moves and maternity/career
  breaks: "exits from this trust", not "left the profession".
- The 2009-2018 file is in FTE while the 2018+ file carries headcount and FTE; the
  splice is built on the FTE rate, which exists in both, so units match across the
  Jul 2018 seam.
