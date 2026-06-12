# ---------------------------------------------------------------------------
# 06_dynamics.R  --  is it a clean glide back, or does it overshoot?
#
# The AR(1) work answered "how much reverts and how fast". This answers a
# different question: the SHAPE of the path. When a trust spikes, does the rate
# glide monotonically back down to its average, or does it overshoot, dropping
# below the mean before drifting back up? The second shape is what stock
# depletion predicts: a big exodus drains the pool of would-be leavers below its
# steady state, so the next period undershoots, then the pool refills.
#
# Test: fit an AR(2) within trust (trust + period FE) and trace the impulse
# response. If the response crosses zero (goes negative after the initial
# positive shock), that is overshoot. If it just decays to zero from above, it
# is monotonic.
#
# One honest confound at lag 1: with an end-of-period denominator, a leaver spike
# shrinks the denominator and can induce a mechanical one-period dip that mimics
# overshoot. So alongside the headline rate we also fit the AR(1) on a rate built
# with the PREVIOUS period's stock as denominator (pre-shock, so artefact-free).
# If the lag-1 negativity survives that, and especially if undershoot shows at
# lag 2+, it is real depletion dynamics, not plumbing.
#
#   in  : data/clean/combined_monthly.rds   (from 05)
#   out : outputs/dynamics_ar2.csv, outputs/figs/impulse_response_annual.png
# ---------------------------------------------------------------------------

pacman::p_load(tidyverse, fixest, broom, here)

cm <- read_rds(here("data", "clean", "combined_monthly.rds"))

# Aggregate to a frequency and build the lags we need, never crossing a gap.
# rate       = leavers / mean stock over the period (the headline measure)
# rate_pre   = leavers / previous period's stock (denominator can't react to the
#              current shock, so a lag-1 dip here is not the mechanical artefact)
to_freq <- function(df, freq) {
  span <- if (freq == "yearly") 12L else 1L
  df |>
    mutate(fy     = if_else(month >= 4, year, year - 1),
           bucket = if (freq == "yearly") fy else year * 100 + month,
           ord    = if (freq == "yearly") fy else year * 12 + (month - 1)) |>
    group_by(staff_group, org_code, bucket, ord) |>
    summarise(n_m = n(), lv = sum(leavers_fte), den = mean(denom_fte), .groups = "drop") |>
    filter(n_m == span) |>
    arrange(staff_group, org_code, ord) |>
    group_by(staff_group, org_code) |>
    mutate(
      step      = ord - lag(ord),
      rate      = lv / den,
      rate_pre  = if_else(step == 1, lv / lag(den), NA_real_),
      l1        = if_else(step == 1, lag(rate), NA_real_),
      l2        = if_else(step == 1 & lag(step) == 1, lag(rate, 2), NA_real_),
      l1_pre    = if_else(step == 1, lag(rate_pre), NA_real_)
    ) |>
    ungroup() |>
    filter(rate <= 1)
}

ar2     <- function(d) feols(rate ~ l1 + l2 | org_code + bucket, d, cluster = ~org_code)
ar1_pre <- function(d) feols(rate_pre ~ l1_pre | org_code + bucket, d, cluster = ~org_code)

# Trace the impulse response of an AR(2): r0 = 1, then r_k = b1 r_{k-1} + b2 r_{k-2}.
# Written as a state recursion fed through accumulate (state = c(latest, previous)),
# so no explicit loop.
irf <- function(b1, b2, h = 8) {
  step <- function(state, .x) c(b1 * state[1] + b2 * state[2], state[1])
  accumulate(seq_len(h), step, .init = c(1, 0)) |> map_dbl(1)
}

groups <- c("Nurses & health visitors",
            "Support to doctors, nurses & midwives",
            "All staff groups")

# Fit AR(2) at annual and monthly, pull the coefficients and the artefact-free
# lag-1, and flag whether the impulse response dips below zero (overshoot).
dyn <- expand_grid(staff_group = groups, freq = c("yearly", "monthly")) |>
  mutate(d = map2(staff_group, freq, ~ to_freq(filter(cm, staff_group == .x), .y))) |>
  mutate(
    m2     = map(d, ~ ar2(filter(.x, !is.na(l1), !is.na(l2)))),
    b1     = map_dbl(m2, ~ coef(.x)[["l1"]]),
    b2     = map_dbl(m2, ~ coef(.x)[["l2"]]),
    b2_p   = map_dbl(m2, ~ tidy(.x) |> filter(term == "l2") |> pull(p.value)),
    b1_pre = map_dbl(d, ~ coef(ar1_pre(filter(.x, !is.na(l1_pre))))[["l1_pre"]]),
    irf    = map2(b1, b2, irf),
    trough = map_dbl(irf, ~ min(.x[-1])),                 # most negative future response
    overshoot = trough < -0.02                            # crosses meaningfully below zero?
  )

dyn_out <- dyn |>
  transmute(staff_group, freq, b1, b2, b2_p, b1_pre, irf_trough = trough, overshoot)
write_csv(dyn_out, here("outputs", "dynamics_ar2.csv"))

# Plot the annual impulse response: where does a 1-unit spike go over the next
# several years? Below the dashed zero line = overshoot.
irf_plot <- dyn |>
  filter(freq == "yearly") |>
  select(staff_group, irf) |>
  mutate(h = map(irf, ~ seq_along(.x) - 1)) |>
  unnest(c(irf, h)) |>
  ggplot(aes(h, irf, colour = staff_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = 0:8) +
  labs(
    title    = "What happens after a leaver-rate spike? (annual)",
    subtitle = "Impulse response of a 1-unit shock. Below the dashed line = overshoot past the mean.",
    x = "Years after the shock", y = "Remaining deviation", colour = NULL
  ) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(here("outputs", "figs", "impulse_response_annual.png"), irf_plot,
       width = 8, height = 5, dpi = 150)

print(dyn_out |> mutate(across(where(is.numeric), ~ round(.x, 3))))
