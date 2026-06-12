# ---------------------------------------------------------------------------
# 01_prep.R  --  build the clean trust x month leaver panel
#
# Takes the NHS Workforce Statistics "turnover from organisation benchmarking"
# extract (one row per org x staff group x month) and turns it into something
# we can run an autoregression on. Nothing clever here, it is mostly
# bookkeeping: parse the period, rebuild a start-of-period denominator, work out
# the leaver rate three different ways (this matters in test 3), and drop the
# cells we can't honestly use.
#
#   in  : data/raw/turnover_org_month.csv   (standard NHS benchmarking schema)
#   out : data/clean/panel.rds              (one tidy row per org x group x month)
# ---------------------------------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, here)

raw <- read_csv(
  here("data", "raw", "turnover_org_month.csv"),
  col_types = cols(
    period_end       = col_integer(),
    org_code         = col_character(),
    org_name         = col_character(),
    nhse_region_name = col_character(),
    ics_code         = col_character(),
    staff_group      = col_character(),
    leavers_hc       = col_double(),
    joiners_hc       = col_double(),
    denom_start_hc   = col_double(),
    denom_end_hc     = col_double(),
    .default         = col_double()
  )
)

# The period arrives as YYYYMM (e.g. 202308). We want two things from it: a real
# date for plotting, and a plain running counter t = 1, 2, 3, ... so we can take
# clean lags. t has to sit on a common calendar across every trust, otherwise
# the lags won't line up, so build it once from the full set of periods rather
# than per series.
periods <- raw |>
  distinct(period_end) |>
  arrange(period_end) |>
  mutate(
    date = ym(as.character(period_end)),
    t    = row_number()
  )

panel <- raw |>
  left_join(periods, by = "period_end") |>
  arrange(org_code, staff_group, t) |>
  group_by(org_code, staff_group) |>
  mutate(
    # NHS only reports an end-of-month headcount reliably (the denom_start column
    # is patchy), so rebuild the start-of-month stock as last month's end stock.
    # The if_else guards against gaps: if the previous row isn't the month
    # immediately before, we don't pretend it is.
    denom_start = if_else(t - lag(t) == 1, lag(denom_end_hc), NA_real_),

    # Three ways of writing the same rate, differing only in the denominator.
    # This is the whole of test 3: if the "reversion" we find shifts depending on
    # which denominator we choose, it is the accounting moving when people leave,
    # not behaviour.
    rate_end   = leavers_hc / denom_end_hc,
    rate_start = leavers_hc / denom_start,
    rate_avg   = leavers_hc / ((denom_end_hc + denom_start) / 2)
  ) |>
  ungroup()

# Honest exclusions, counted so we can report what fell out.
#  - suppressed leaver counts (<5) come through blank, can't impute them
#  - a couple of tiny ICB cells report a rate above 100%, drop them
#  - very short series (mergers, late ESR joiners) can't support a within-trust
#    autoregression, so require a decent run of months
before <- nrow(panel)

clean <- panel |>
  filter(
    !is.na(leavers_hc),
    !is.na(denom_end_hc), denom_end_hc > 0,
    !is.na(rate_end), rate_end <= 1
  ) |>
  group_by(org_code, staff_group) |>
  filter(n() >= 12) |>
  ungroup()

dir.create(here("data", "clean"), showWarnings = FALSE, recursive = TRUE)
write_rds(clean, here("data", "clean", "panel.rds"))

# Attrition note for whoever is QAing this.
attrition <- tibble(
  rows_raw   = before,
  rows_clean = nrow(clean),
  dropped    = before - nrow(clean),
  orgs       = n_distinct(clean$org_code),
  groups     = n_distinct(clean$staff_group),
  months     = n_distinct(clean$t)
)
print(attrition)
