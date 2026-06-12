# ---------------------------------------------------------------------------
# 02_analysis.R  --  is the year-to-year wobble in trust leaver rates real?
#
# The question behind the project: when a trust has a bad month/year for leavers
# and then a better one, is that
#   (a) genuine MEAN REVERSION  -- the pool of people who wanted to leave has
#       been drained, and you can't leave twice, so next period really is lower; or
#   (b) just REGRESSION TO THE MEAN -- last period's bad number was mostly noise,
#       and noise doesn't repeat, so the trust never actually changed.
# They look identical in the raw data. The point of this script is to pull them
# apart, because the policy reading is different: (a) says a bad year is real but
# partly self-curing; (b) says a bad year was barely a signal at all. Both warn
# against ranking trusts on a single year, which is the thesis.
#
# The spine is one comparison, run per staff group:
#   b_obs       : the persistence we actually measure within a trust
#   reliability : the share of a trust's month-to-month variation that is real
#                 signal rather than the counting noise of leavers out of a
#                 finite headcount
#   b_true      : b_obs / reliability, the persistence once the noise-driven
#                 attenuation is undone
# and then the observed bounce (1 - b_obs) splits cleanly into
#   genuine reversion (1 - b_true)  +  noise / RTM  (b_true - b_obs).
# ---------------------------------------------------------------------------

pacman::p_load(tidyverse, fixest, broom, scales, here)

panel <- read_rds(here("data", "clean", "panel.rds")) |>
  # one id per trust-and-group series, used as the panel unit in pooled fits
  mutate(series = paste(org_code, staff_group, sep = " | "))

# Add the one-month lag of each rate, within each series. Doing it by hand with
# dplyr rather than fixest's l() keeps it obvious to a reviewer which row is
# being treated as "last month", and the t - lag(t) == 1 check means we never
# lag across a missing month.
lagged <- panel |>
  arrange(series, t) |>
  group_by(series) |>
  mutate(
    gap       = t - lag(t),
    lag_end   = if_else(gap == 1, lag(rate_end),   NA_real_),
    lag_start = if_else(gap == 1, lag(rate_start), NA_real_),
    lag_avg   = if_else(gap == 1, lag(rate_avg),   NA_real_)
  ) |>
  ungroup()

# The workhorse: within-trust AR(1) with trust and period fixed effects. The
# trust FE soak up "some trusts are just leakier than others"; the period FE soak
# up national seasonality (the autumn churn, the August doctor rotation) and any
# common shock, so the coefficient is pure within-trust persistence. Cluster on
# the panel unit.
fit_ar <- function(df, y = "rate_end", lg = "lag_end", unit = "series") {
  f <- as.formula(paste0(y, " ~ ", lg, " | ", unit, " + period_end"))
  feols(f, data = df, cluster = as.formula(paste0("~", unit)))
}

# --- Test 1: the raw measured persistence, per staff group -------------------
by_group <- lagged |>
  group_by(staff_group) |>
  nest() |>
  mutate(model = map(data, ~ fit_ar(.x, unit = "org_code")))

ar_table <- by_group |>
  mutate(tidy = map(model, tidy)) |>
  select(staff_group, tidy) |>
  unnest(tidy) |>
  filter(term == "lag_end") |>
  transmute(staff_group, b_raw = estimate, se = std.error)

# --- Test 1b: undo the Nickell bias with a half-panel jackknife --------------
# A within-group AR is biased downward in a short panel (Nickell 1981), which
# would hand us "mean reversion" for free. T is ~32 months so the bias is modest,
# but correct it anyway: estimate on the full window and on each time-half, then
#   b_corrected = 2 * b_full - mean(b_first_half, b_second_half),
# which cancels the 1/T bias to first order (Dhaene & Jochmans 2015). No package.
ar_coef <- function(df) coef(fit_ar(df, unit = "org_code"))[["lag_end"]]

hpj_one <- function(df) {
  cut <- median(df$t)
  full   <- ar_coef(df)
  first  <- ar_coef(filter(df, t <= cut))
  second <- ar_coef(filter(df, t >  cut))
  tibble(b_full = full, b_hpj = 2 * full - (first + second) / 2)
}

bias_table <- by_group |>
  mutate(hpj = map(data, hpj_one)) |>
  select(staff_group, hpj) |>
  unnest(hpj)

# --- Test 2: reliability, and the decomposition ------------------------------
# How much of a trust's within-series variation is real signal rather than pure
# counting noise? Leavers are a count out of a finite stock, so even a trust
# whose true rate never moved would wobble from the arithmetic of small numbers.
# That sampling-noise floor is rate*(1-rate)/stock. This is the conservative
# version of the test on purpose: sampling noise is the *smallest* the
# measurement error can be (real life adds reclassification, maternity lumpiness,
# and the rest), so if even this floor explains most of the bounce, the case for
# genuine reversion is weak rather than strong.
reliability_tbl <- lagged |>
  filter(!is.na(lag_end)) |>                 # the same sample the AR runs on
  group_by(staff_group, org_code) |>
  mutate(dev = rate_end - mean(rate_end)) |> # within-trust deviation
  ungroup() |>
  mutate(noise = rate_end * (1 - rate_end) / denom_end_hc) |>
  group_by(staff_group) |>
  summarise(
    var_total   = mean(dev^2),
    noise_bar   = mean(noise),
    reliability = (var_total - noise_bar) / var_total,
    .groups = "drop"
  )

decomp <- ar_table |>
  left_join(bias_table,      by = "staff_group") |>
  left_join(reliability_tbl, by = "staff_group") |>
  mutate(
    b_true            = b_hpj / reliability,   # the bias-corrected, noise-corrected persistence
    reversion_obs     = 1 - b_hpj,             # the bounce we actually see
    reversion_genuine = 1 - b_true,            # real mean reversion (the world reverting)
    reversion_noise   = b_true - b_hpj         # the RTM / noise share (measurement reverting)
  )

# --- Test 3: does the answer depend on the denominator? ----------------------
# Refit the AR with the start-stock and average-stock rates and line the
# coefficients up. If they wander, part of the "reversion" is just the
# denominator shrinking when people leave, which is plumbing, not behaviour.
denom_robust <- tibble(
  construction = c("end stock", "start stock", "average stock"),
  y            = c("rate_end",  "rate_start",  "rate_avg"),
  lg           = c("lag_end",   "lag_start",   "lag_avg")
) |>
  mutate(b = map2_dbl(y, lg, ~ coef(fit_ar(lagged, y = .x, lg = .y))[[.y]]))

# --- Test 4: is the snap-back symmetric? -------------------------------------
# The policy point cuts both ways: trusts above their own average should drift
# down, but unusually good years should drift back up too. Interact the lag with
# whether the trust sat above its own mean last period.
asym <- lagged |>
  filter(!is.na(lag_end)) |>
  group_by(series) |>
  mutate(above_own_mean = lag_end > mean(lag_end)) |>
  ungroup() |>
  feols(rate_end ~ i(above_own_mean, lag_end, ref = FALSE) | series + period_end,
        data = _, cluster = ~series)

# --- write everything out ----------------------------------------------------
dir.create(here("outputs"), showWarnings = FALSE, recursive = TRUE)
write_csv(decomp,       here("outputs", "decomposition.csv"))
write_csv(denom_robust, here("outputs", "denominator_robustness.csv"))
write_csv(tidy(asym),   here("outputs", "asymmetry.csv"))

# --- two pictures for the deck ----------------------------------------------
# 1. the pattern itself: this month against last month, with the 45-degree line.
#    Points sitting below the line on the right = high months followed by lower.
p_scatter <- lagged |>
  filter(!is.na(lag_end)) |>
  slice_sample(n = min(4000, sum(!is.na(lagged$lag_end)))) |>
  ggplot(aes(lag_end, rate_end)) +
  geom_point(alpha = 0.12, size = 0.7) +
  geom_abline(linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ staff_group) +
  scale_x_continuous(labels = label_percent()) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title    = "Leaver rate this month against last month",
    subtitle = "Dashed = no change. The fit line below it on the right is the bounce-back.",
    x = "Leaver rate, last month", y = "Leaver rate, this month"
  ) +
  theme_minimal(base_size = 12)

# 2. the decomposition: of the bounce we see, how much is real vs noise.
p_decomp <- decomp |>
  select(staff_group, reversion_genuine, reversion_noise) |>
  pivot_longer(-staff_group, names_to = "source", values_to = "share") |>
  mutate(source = recode(source,
                         reversion_genuine = "Genuine mean reversion",
                         reversion_noise   = "Regression to the mean (noise)")) |>
  ggplot(aes(staff_group, share, fill = source)) +
  geom_col() +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title = "What is the bounce-back made of?",
    x = NULL, y = "Share of the observed reversion", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave(here("outputs", "figs", "scatter_rate_vs_lag.png"), p_scatter,
       width = 9, height = 4.5, dpi = 150)
ggsave(here("outputs", "figs", "decomposition.png"), p_decomp,
       width = 7, height = 4.5, dpi = 150)

# --- the headline tables, to the console ------------------------------------
print(decomp)
print(denom_robust)
print(tidy(asym))
