# ---------------------------------------------------------------------------
# 05_combined.R  --  splice 2009-2018 onto 2018+ and re-run the term structure
#
# The 2018+ build gave us a thin yearly column (T=7). The supplementary file
# "Turnover from organisation, by staff group, 2009 to 2018" backfills the seven
# years before it, monthly, at org level. The two meet at Jul 2018 (one month
# overlap), so spliced we get a continuous Oct 2009 to Mar 2026 series, about 16
# financial years, which is what the yearly autoregression needed.
#
# One wrinkle: the 2009-2018 file is in FTE (fractional leavers and denominator),
# while the 2018+ file carries both headcount and FTE. So build the whole thing
# on the FTE leaver rate, which exists in both, and the units match across the
# seam. FTE and headcount leaver rates are near-identical in practice, so this
# changes nothing dynamically, it just lets us join cleanly.
#
#   in  : data/raw/turnover_2009_2018_staffgroup.xlsx  (FTE, Oct 2009 to Jul 2018)
#         data/raw/turnover_org_month_full.csv          (Jul 2018 to Mar 2026)
#   out : outputs/combined_term_structure.csv
# ---------------------------------------------------------------------------

pacman::p_load(tidyverse, readxl, fixest, broom, here)

groups_we_want <- c("Nurses & health visitors",
                    "Support to doctors, nurses & midwives",
                    "All staff groups")

harmonise <- function(x) {
  x |>
    str_replace("^All staff group$", "All staff groups") |>
    (\(s) if_else(str_detect(s, regex("hchs +doctors", ignore_case = TRUE)),
                  "HCHS doctors", s))()
}

# --- 2009 to 2018, FTE ---  pull the five columns we need by position and parse
# the period end out of the "200909 to 200910" label (the second YYYYMM).
early <- read_excel(
  here("data", "raw", "turnover_2009_2018_staffgroup.xlsx"),
  sheet = "Data", skip = 1, col_names = FALSE, .name_repair = "minimal"
) |>
  select(period_lbl = 1, org_code = 7, staff_group = 10,
         leavers_fte = 13, denom_fte = 14) |>
  mutate(
    period_end  = as.integer(str_extract(period_lbl, "\\d+$")),
    staff_group = harmonise(staff_group),
    across(c(leavers_fte, denom_fte), as.numeric)
  ) |>
  select(period_end, org_code, staff_group, leavers_fte, denom_fte)

# --- 2018 to 2026, FTE ---  same five fields out of our existing build
late <- read_csv(
  here("data", "raw", "turnover_org_month_full.csv"),
  col_types = cols(period_end = col_integer(), org_code = col_character(),
                   staff_group = col_character(), leavers_fte = col_double(),
                   denom_end_fte = col_double(), .default = col_character())
) |>
  transmute(period_end, org_code,
            staff_group = harmonise(staff_group),
            leavers_fte, denom_fte = denom_end_fte)

# Stack them, keep the later vintage where the one overlap month (Jul 2018)
# appears in both, and keep provider trusts (R codes) only.
combined <- bind_rows(late |> mutate(src = "late"),
                      early |> mutate(src = "early")) |>
  filter(staff_group %in% groups_we_want, str_starts(org_code, "R"),
         !is.na(denom_fte), denom_fte > 0, !is.na(leavers_fte)) |>
  distinct(period_end, org_code, staff_group, .keep_all = TRUE) |>   # late wins (listed first)
  mutate(
    year  = period_end %/% 100,
    month = period_end %%  100,
    rate  = leavers_fte / denom_fte
  ) |>
  filter(rate <= 1)

# stash the spliced monthly panel for the dynamics script (06)
dir.create(here("data", "clean"), showWarnings = FALSE, recursive = TRUE)
write_rds(combined, here("data", "clean", "combined_monthly.rds"))

# tell the reviewer what the spliced series actually spans
combined |>
  summarise(periods = n_distinct(period_end),
            first = min(period_end), last = max(period_end),
            trusts = n_distinct(org_code), rows = n()) |>
  print()

# --- frequency machinery (same as 04, on the FTE rate) ----------------------
bucketise <- function(df, freq) {
  q <- (df$month - 1) %/% 3 + 1
  h <- (df$month - 1) %/% 6 + 1
  fy <- if_else(df$month >= 4, df$year, df$year - 1)
  df |>
    mutate(
      span   = c(monthly = 1, quarterly = 3, halfyear = 6, yearly = 12)[[freq]],
      bucket = switch(freq, monthly = year * 100 + month, quarterly = year * 10 + q,
                      halfyear = year * 10 + h, yearly = fy),
      ord    = switch(freq, monthly = year * 12 + (month - 1), quarterly = year * 4 + (q - 1),
                      halfyear = year * 2 + (h - 1), yearly = fy)
    )
}
to_freq <- function(df, freq) {
  bucketise(df, freq) |>
    group_by(staff_group, org_code, bucket, ord, span) |>
    summarise(n_m = n(), lv = sum(leavers_fte), den = mean(denom_fte), .groups = "drop") |>
    filter(n_m == span) |>
    mutate(rate = lv / den) |>
    filter(!is.na(rate), rate <= 1) |>
    arrange(staff_group, org_code, ord) |>
    group_by(staff_group, org_code) |>
    mutate(lag_rate = if_else(ord - lag(ord) == 1, lag(rate), NA_real_)) |>
    ungroup()
}
fit_ar  <- function(d) feols(rate ~ lag_rate | org_code + bucket, d, cluster = ~org_code)
ar_coef <- function(d) { o <- try(coef(fit_ar(d))[["lag_rate"]], silent = TRUE)
                         if (inherits(o, "try-error")) NA_real_ else o }
hpj <- function(d) { cut <- median(d$ord)
  2 * ar_coef(d) - (ar_coef(filter(d, ord <= cut)) + ar_coef(filter(d, ord > cut))) / 2 }
reliability <- function(d) { r <- resid(feols(rate ~ 1 | org_code + bucket, d))
  (mean(r^2) - mean(d$rate * (1 - d$rate) / d$den)) / mean(r^2) }

combined_term <- expand_grid(
  staff_group = groups_we_want,
  freq        = c("monthly", "quarterly", "halfyear", "yearly")
) |>
  mutate(d = map2(staff_group, freq, ~ to_freq(filter(combined, staff_group == .x), .y) |>
                    filter(!is.na(lag_rate)))) |>
  mutate(
    T_bar = map_dbl(d, ~ mean(table(.x$org_code))),
    n_obs = map_int(d, nrow),
    b_raw = map_dbl(d, ar_coef),
    b_hpj = map_dbl(d, hpj),
    reliability = map_dbl(d, reliability)
  ) |>
  mutate(
    b_true     = b_hpj / reliability,
    span_m     = c(monthly = 1, quarterly = 3, halfyear = 6, yearly = 12)[freq],
    halflife_m = if_else(b_true > 0 & b_true < 1, log(0.5) / log(b_true) * span_m, NA_real_)
  ) |>
  select(staff_group, freq, T_bar, n_obs, b_raw, b_hpj, reliability, b_true, halflife_m)

write_csv(combined_term, here("outputs", "combined_term_structure.csv"))
print(combined_term |> mutate(across(where(is.numeric), ~ round(.x, 3))), n = 40)
