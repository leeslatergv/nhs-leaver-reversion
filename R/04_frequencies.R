# ---------------------------------------------------------------------------
# 04_frequencies.R  --  the term structure of reversion
#
# Lee's question from the start: does the leaver rate revert monthly, quarterly,
# half-yearly, yearly? The persistence is horizon-dependent and that is the
# interesting bit. So aggregate the monthly series to each frequency and run the
# same within-trust AR(1) at each, with the Nickell jackknife and the reliability
# split, and line the four up.
#
# The depletion mechanism (Lee's version: if you left this period for ANY reason
# you are not in the pool to leave next period) is a horizon story, so we expect
# the reversion to look different as the window widens.
#
# Runs on the 2018+ all-groups series. The yearly column will be thin (T=7) and
# that is the point: it shows exactly why the 2009-2023 back-series is worth
# grabbing. Splice that in and re-run this same script and the yearly row firms up.
#
#   in  : data/raw/turnover_org_month_full.csv
#   out : outputs/term_structure.csv
# ---------------------------------------------------------------------------

pacman::p_load(tidyverse, fixest, broom, here)

mono <- read_csv(
  here("data", "raw", "turnover_org_month_full.csv"),
  col_types = cols(period_end = col_integer(), org_code = col_character(),
                   staff_group = col_character(), leavers_hc = col_double(),
                   denom_end_hc = col_double(), .default = col_character())
) |>
  mutate(staff_group = if_else(
    str_detect(staff_group, regex("hchs +doctors", ignore_case = TRUE)),
    "HCHS doctors", staff_group)) |>
  filter(
    str_starts(org_code, "R"),
    !is.na(denom_end_hc), denom_end_hc > 0,
    staff_group %in% c("Nurses & health visitors",
                       "Support to doctors, nurses & midwives",
                       "All staff groups")
  ) |>
  mutate(year = period_end %/% 100, month = period_end %% 100)

# Tag each month with its bucket for a given frequency, plus a continuous ordinal
# (increments by exactly 1 each period) so lags never cross a gap, and the span
# (how many months make a complete bucket). Yearly uses the NHS financial year.
bucketise <- function(df, freq) {
  q <- (df$month - 1) %/% 3 + 1
  h <- (df$month - 1) %/% 6 + 1
  fy <- if_else(df$month >= 4, df$year, df$year - 1)
  df |>
    mutate(
      span   = c(monthly = 1, quarterly = 3, halfyear = 6, yearly = 12)[[freq]],
      bucket = switch(freq,
                      monthly   = year * 100 + month,
                      quarterly = year * 10  + q,
                      halfyear  = year * 10  + h,
                      yearly    = fy),
      ord    = switch(freq,
                      monthly   = year * 12 + (month - 1),
                      quarterly = year * 4  + (q - 1),
                      halfyear  = year * 2  + (h - 1),
                      yearly    = fy)
    )
}

# Aggregate to the frequency: a bucket only counts if every constituent month is
# present and unsuppressed (otherwise the leaver sum is wrong). Rate = leavers
# over the window divided by the average headcount across it.
to_freq <- function(df, freq) {
  bucketise(df, freq) |>
    group_by(staff_group, org_code, bucket, ord, span) |>
    summarise(n_m = n(), n_sup = sum(is.na(leavers_hc)),
              lv = sum(leavers_hc), den = mean(denom_end_hc), .groups = "drop") |>
    filter(n_m == span, n_sup == 0) |>
    mutate(rate = lv / den) |>
    filter(!is.na(rate), rate <= 1) |>
    arrange(staff_group, org_code, ord) |>
    group_by(staff_group, org_code) |>
    mutate(lag_rate = if_else(ord - lag(ord) == 1, lag(rate), NA_real_)) |>
    ungroup()
}

fit_ar  <- function(d) feols(rate ~ lag_rate | org_code + bucket, d, cluster = ~org_code)
ar_coef <- function(d) {
  out <- try(coef(fit_ar(d))[["lag_rate"]], silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else out
}
hpj <- function(d) {
  cut <- median(d$ord)
  2 * ar_coef(d) - (ar_coef(filter(d, ord <= cut)) + ar_coef(filter(d, ord > cut))) / 2
}
reliability <- function(d) {
  r <- resid(feols(rate ~ 1 | org_code + bucket, d))
  (mean(r^2) - mean(d$rate * (1 - d$rate) / d$den)) / mean(r^2)
}

# Run every group x frequency combination. Half-life in periods comes from the
# noise-corrected persistence: how many periods for a shock to fall by half.
grid <- expand_grid(
  staff_group = c("All staff groups", "Nurses & health visitors",
                  "Support to doctors, nurses & midwives"),
  freq        = c("monthly", "quarterly", "halfyear", "yearly")
)

term_structure <- grid |>
  mutate(d = map2(staff_group, freq, ~ to_freq(filter(mono, staff_group == .x), .y) |>
                    filter(!is.na(lag_rate)))) |>
  mutate(
    T_bar       = map_dbl(d, ~ mean(table(.x$org_code))),
    n_obs       = map_int(d, nrow),
    b_raw       = map_dbl(d, ar_coef),
    b_hpj       = map_dbl(d, hpj),
    reliability = map_dbl(d, reliability)
  ) |>
  mutate(
    b_true     = b_hpj / reliability,
    span_m     = c(monthly = 1, quarterly = 3, halfyear = 6, yearly = 12)[freq],
    halflife_m = if_else(b_true > 0 & b_true < 1,
                         log(0.5) / log(b_true) * span_m, NA_real_)
  ) |>
  select(staff_group, freq, T_bar, n_obs, b_raw, b_hpj, reliability, b_true, halflife_m)

dir.create(here("outputs"), showWarnings = FALSE, recursive = TRUE)
write_csv(term_structure, here("outputs", "term_structure.csv"))

print(term_structure |> mutate(across(where(is.numeric), ~ round(.x, 3))), n = 40)
