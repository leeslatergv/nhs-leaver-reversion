# ---------------------------------------------------------------------------
# 03_annual.R  --  the actual test of the depletion hypothesis
#
# The monthly pass (02) found no within-trust persistence, but the depletion
# story ("a bad year drains the pool of would-be leavers, so next year is
# lower") is an ANNUAL claim, and 32 monthly points couldn't test it. This uses
# the full back-series (Jul 2018 to Mar 2026, all staff groups) collapsed to
# financial years, so we finally have a year-on-year autoregression to look at.
#
# Same logic as before, just at annual horizon:
#   b_obs       persistence of the annual rate, within trust
#   b_hpj       the same, with the Nickell short-panel bias jackknifed out
#               (this matters a lot now: T = 7 years, so the raw bias is ~ -1/6)
#   reliability real signal share, from the annual counting-noise floor
#   b_true      b_hpj / reliability
# Depletion predicts b_true well below 1 (a high year really does pull back),
# and reliability tells us whether that pull-back is real or just noise.
#
#   in  : data/raw/turnover_org_month_full.csv
#   out : outputs/annual_*.csv, outputs/figs/annual_*.png
# ---------------------------------------------------------------------------

pacman::p_load(tidyverse, fixest, broom, scales, here)

raw <- read_csv(
  here("data", "raw", "turnover_org_month_full.csv"),
  col_types = cols(
    period_end     = col_integer(),
    org_code       = col_character(),
    staff_group    = col_character(),
    leavers_hc     = col_double(),
    denom_end_hc   = col_double(),
    .default       = col_character()
  )
) |>
  mutate(across(c(leavers_hc, denom_end_hc), as.numeric))

# The doctor category was relabelled across vintages (junior -> resident
# doctors, plus capitalisation drift), so collapse every variant to one label,
# otherwise the same group looks like a break. Then keep the substantive groups
# and the provider trusts (R codes), which also drops the tiny ICB cells that
# throw the >100% rates.
prepped <- raw |>
  mutate(staff_group = if_else(
    str_detect(staff_group, regex("hchs +doctors", ignore_case = TRUE)),
    "HCHS doctors", staff_group
  )) |>
  filter(
    staff_group %in% c(
      "All staff groups", "Nurses & health visitors",
      "Support to doctors, nurses & midwives", "Midwives", "HCHS doctors",
      "Scientific, therapeutic & technical staff", "Support to ST&T staff",
      "Ambulance staff", "Managers", "Senior managers"
    ),
    str_starts(org_code, "R"),
    !is.na(denom_end_hc), denom_end_hc > 0
  ) |>
  mutate(
    year  = period_end %/% 100,
    month = period_end %%  100,
    fy    = if_else(month >= 4, year, year - 1)   # NHS financial year, Apr to Mar
  )

# Collapse to financial years. A year only counts if all 12 months are present
# and none of its leaver counts were suppressed, otherwise the annual sum is
# wrong rather than just noisy. Three denominators again for the robustness test:
# average headcount over the year, the April stock (start), the March stock (end).
fy_panel <- prepped |>
  group_by(org_code, staff_group, fy) |>
  summarise(
    n_months   = n(),
    n_suppress = sum(is.na(leavers_hc)),
    leavers_fy = sum(leavers_hc),
    denom_avg  = mean(denom_end_hc),
    denom_apr  = denom_end_hc[month == 4][1],
    denom_mar  = denom_end_hc[month == 3][1],
    .groups = "drop"
  ) |>
  filter(n_months == 12, n_suppress == 0) |>
  mutate(
    rate_avg   = leavers_fy / denom_avg,
    rate_start = leavers_fy / denom_apr,
    rate_end   = leavers_fy / denom_mar,
    series     = paste(org_code, staff_group, sep = " | ")
  ) |>
  filter(!is.na(rate_avg), rate_avg <= 1)

# One-year lags, within each trust-and-group series, never lagging across a gap.
lagged <- fy_panel |>
  arrange(series, fy) |>
  group_by(series) |>
  mutate(
    gap       = fy - lag(fy),
    lag_avg   = if_else(gap == 1, lag(rate_avg),   NA_real_),
    lag_start = if_else(gap == 1, lag(rate_start), NA_real_),
    lag_end   = if_else(gap == 1, lag(rate_end),   NA_real_)
  ) |>
  ungroup()

# Annual AR(1), trust + year fixed effects, clustered on trust.
fit_ar <- function(df, y = "rate_avg", lg = "lag_avg", unit = "org_code") {
  f <- as.formula(paste0(y, " ~ ", lg, " | ", unit, " + fy"))
  feols(f, data = df, cluster = as.formula(paste0("~", unit)))
}
ar_coef <- function(df) {
  out <- try(coef(fit_ar(df))[["lag_avg"]], silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else out
}

# Half-panel jackknife for the Nickell bias. Split the years in two, estimate on
# each half, and b_hpj = 2*full - mean(halves). Essential at T = 7.
hpj_one <- function(df) {
  cut <- median(df$fy)
  full   <- ar_coef(df)
  first  <- ar_coef(filter(df, fy <= cut))
  second <- ar_coef(filter(df, fy >  cut))
  2 * full - (first + second) / 2
}

# Reliability the honest way this time: residualise the rate on trust AND year
# (the two-way within, matching what the AR actually uses), then compare its
# variance to the annual counting-noise floor rate*(1-rate)/headcount.
reliability_one <- function(df) {
  r <- resid(feols(rate_avg ~ 1 | org_code + fy, data = df))
  var_total <- mean(r^2)
  noise_bar <- mean(df$rate_avg * (1 - df$rate_avg) / df$denom_avg)
  (var_total - noise_bar) / var_total
}

# Run the whole lot per staff group, keeping only groups with enough trusts for
# the within estimator to mean anything.
annual <- lagged |>
  filter(!is.na(lag_avg)) |>
  group_by(staff_group) |>
  nest() |>
  mutate(n_trusts = map_int(data, ~ n_distinct(.x$org_code)),
         n_obs    = map_int(data, nrow)) |>
  filter(n_trusts >= 50) |>
  mutate(
    b_raw       = map_dbl(data, ar_coef),
    b_hpj       = map_dbl(data, hpj_one),
    reliability = map_dbl(data, reliability_one),
    se          = map_dbl(data, ~ tidy(fit_ar(.x)) |>
                            filter(term == "lag_avg") |> pull(std.error))
  ) |>
  mutate(b_true = b_hpj / reliability) |>
  select(staff_group, n_trusts, n_obs, b_raw, se, b_hpj, reliability, b_true) |>
  arrange(desc(n_obs))

# Denominator robustness, on nurses (the cleanest, densest group).
nurse <- lagged |> filter(staff_group == "Nurses & health visitors", !is.na(lag_avg))
denom_robust <- tibble(
  construction = c("average stock", "start (Apr) stock", "end (Mar) stock"),
  y  = c("rate_avg",  "rate_start", "rate_end"),
  lg = c("lag_avg",   "lag_start",  "lag_end")
) |>
  mutate(b = map2_dbl(y, lg, ~ coef(fit_ar(nurse, y = .x, lg = .y))[[.y]]))

# Asymmetry, nurses: does a year above the trust's own mean snap back harder
# than a year below it?
asym <- nurse |>
  group_by(series) |>
  mutate(above_own_mean = lag_avg > mean(lag_avg)) |>
  ungroup() |>
  feols(rate_avg ~ i(above_own_mean, lag_avg, ref = FALSE) | series + fy,
        data = _, cluster = ~series)

# COVID sensitivity: FY2020/21 and FY2021/22 were structurally odd. The year FE
# already soak up the common shock, but drop them entirely and re-estimate the
# nurse AR to be sure the result isn't a COVID artefact.
nurse_excovid <- nurse |> filter(!fy %in% c(2020, 2021))
covid_check <- tibble(
  sample = c("all years", "ex-COVID (drop FY20/21 + FY21/22)"),
  b_raw  = c(ar_coef(nurse), ar_coef(nurse_excovid)),
  b_hpj  = c(hpj_one(nurse), hpj_one(nurse_excovid))
)

# write outputs
dir.create(here("outputs"), showWarnings = FALSE, recursive = TRUE)
write_csv(annual,       here("outputs", "annual_by_group.csv"))
write_csv(denom_robust, here("outputs", "annual_denominator_robustness.csv"))
write_csv(tidy(asym),   here("outputs", "annual_asymmetry.csv"))
write_csv(covid_check,  here("outputs", "annual_covid_check.csv"))

# picture: this year's nurse leaver rate against last year's, with the no-change
# line. Sitting below the line on the right = high years followed by lower ones.
p <- nurse |>
  ggplot(aes(lag_avg, rate_avg)) +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_abline(linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_continuous(labels = label_percent()) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title    = "Nurse leaver rate: this financial year against last",
    subtitle = "Provider trusts, FY2019/20 to FY2025/26. Dashed = no change.",
    x = "Leaver rate, last year", y = "Leaver rate, this year"
  ) +
  theme_minimal(base_size = 12)
ggsave(here("outputs", "figs", "annual_nurse_scatter.png"), p,
       width = 7, height = 5, dpi = 150)

print(annual)
print(denom_robust)
print(covid_check)
print(tidy(asym))
