#| message: false
library(lubridate)
library(janitor)
library(stringr)
library(forcats)
library(colorspace)
library(ggnewscale)
library(rstatix)
library(mgcv)
library(performance)
library(brm)
library(pdp)
library(tidyr)
library(caret)
library(ggrepel)
library(scales)
library(cowplot)
library(ggnewscale)
library(patchwork)
library(doParallel)
library(purrr)
library(grid)
library(recipes)
library(gam.hp)
library(gt)
library(scales) # for alpha()
library(ggplot2)


theme_set(theme_minimal(base_size = 12))

# file

df_master <- read.csv("raw_df_Alpine_2025.csv")

# Add C/N ratio column
df_master <- df_master %>%
  mutate(C_N_ratio = total_C / total_N)

head(df_master)

conv_factor <- 1e6 / 12 / 86400  # ≈ 0.9645


grp_cols <- c(
  "Erect Shrub"     = "#045b5a",
  "Prostrate Shrub" = "#1f7e44",
  "Forb"            = "#abdf26",
  "Graminoid"       = "#dfd626",
  "Barren"          = "#733421"
)

# ------------ Clean up ----

resp_ch4 <- "ch4_mgCH4_m2_d"
impute_missing <- TRUE

dfL <- df_master %>% 
  dplyr::filter(light_dark == "L")

dfmL <- dfL %>%
  mutate(
    prop_shrub_sum = rowSums(
      cbind(prop_Erect.Shrub, prop_Prostrate.Shrub),
      na.rm = TRUE
    ),
    prop_shrub_sum = ifelse(
      is.na(prop_Erect.Shrub) & is.na(prop_Prostrate.Shrub),
      NA,
      prop_shrub_sum
    ),
    prop_shrub = dplyr::coalesce(
      prop_shrub_sum,
      cover_shrub / 100
    ),
    prop_graminoid = dplyr::coalesce(prop_Graminoid, cover_graminoid/100),
    prop_forb      = dplyr::coalesce(prop_Forb,      cover_forb/100),
    
    even3          = dplyr::coalesce(evenness),
    rich3          = dplyr::coalesce(as.numeric(richness))
  ) %>%
  dplyr::select(-prop_shrub_sum)

candidates_base <- c(
  "ch4_mgCH4_m2_d",
  "soil_moisture_avg",
  "x5cm_soil_temp",
  "x12cm_soil_temp",
  "surface_temp",
  "air_temp",
  "avg_par",
  "chamber_mean_temp",
  "cover_shrub",
  "cover_graminoid",
  "cover_forb",
  "prop_shrub",
  "prop_graminoid",
  "prop_forb",
  "rich3","even3",
  "gcc_mean",
  "gcc_med","gcc_p90","exg_mean","vari_mean",
  "total_n","p","cu","p_h",
  "ER_gC_m2_d","GPP_gC_m2_d","light_dark","co2_gC_m2_d",
  "bryo_lichen_sum",
  "bulk_density_g_cm3","slope_aspect",
  "no3_n","nh4_n",
  "soil_depth_max",
  "soil_depth_min",
  "veg_height_mean",
  "veg_height_max",
  "soil_depth_avg",
  "C_N_ratio"
) %>% unique()

# PRS and biomass 
prs_cols <- grep("^(no3|nh4|total|p|k|ca|mg|s|fe|mn|cu|zn|b|al|cd|p_h)\\b",
                 names(dfmL), value = TRUE, ignore.case = TRUE)
bm_cols  <- grep("^bm_", names(dfmL), value = TRUE)

preds <- intersect(unique(c(candidates_base, prs_cols, bm_cols)), names(dfmL))
stopifnot(length(preds) > 1)

# ---------------- GCC -----------------

df_greenness_long <- df_master %>%
  select(functional_group,
         gcc_mean, gcc_med, gcc_p90,
         exg_mean, vari_mean) %>%
  pivot_longer(
    cols = c(gcc_mean, gcc_med, gcc_p90, exg_mean, vari_mean),
    names_to = "index",
    values_to = "value"
  )

# Visualize greenness across functional groups
ggplot(df_greenness_long, 
       aes(x = functional_group, y = value, fill = functional_group)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 21) +
  scale_fill_manual(values = grp_cols) +
  facet_wrap(~ index, scales = "free_y") +
  labs(
    x = "Functional Group",
    y = "Greenness",
    title = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


df_rel <- df_master %>%
  select(
    gcc_mean, gcc_med, gcc_p90, exg_mean, vari_mean,
    prop_Prostrate.Shrub,
    prop_Barren,
    prop_Graminoid,
    prop_Forb,
    prop_Erect.Shrub
  ) %>%
  pivot_longer(
    cols = starts_with("prop_"),
    names_to = "veg_type",
    values_to = "prop_cover"
  ) %>%
  pivot_longer(
    cols = c(gcc_mean, gcc_med, gcc_p90, exg_mean, vari_mean),
    names_to = "greenness_metric",
    values_to = "greenness_value"
  )


ggplot(df_rel, aes(x = prop_cover, y = greenness_value)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_grid(greenness_metric ~ veg_type, scales = "free_y") +
  labs(
    x = "Proportional Cover",
    y = "Greenness Metric Value"  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )


# ---------- Light Model: PRS n~40 ------------------

safe_k <- function(x, k_max = 6L) {
  nuniq <- length(unique(x[is.finite(x)]))
  as.integer(max(3L, min(k_max, nuniq - 1L)))
}


glimpse(df_master)
glimpse(dfmL)

df_clean <- dfmL %>%
  rename(
    prop_prostrate_shrub =  prop_Prostrate.Shrub,
    prop_erect_shrub     =  prop_Erect.Shrub,
    cover_prostrate_shrub=  cover_Prostrate.Shrub,
    cover_erect_shrub    =  cover_Erect.Shrub
  )


# Core non-PRS "drivers"
core_vars <- c(
  "ch4_mgCH4_m2_d",          # response
  "soil_moisture_avg",
  "x12cm_soil_temp",
  "avg_par",
  "x5cm_soil_temp",
  "bm_prop_dead",
  "prop_shrub","prop_forb","prop_graminoid","prop_Barren",
  "gcc_mean",
  "GPP_gC_m2_d","ER_gC_m2_d","co2_gC_m2_d", "C_N_ratio",
  "slope_aspect",
  "bulk_density_g_cm3", "volumetric_water_content_percent",
  "veg_height_mean","soil_depth_avg"
)

# PRS 
prs_vars <- grep("^(total_n|no3_n|nh4_n|p|k|ca|mg|s|fe|mn|cu|zn|b|al|cd|p_h)$",
                 names(df_clean), value = TRUE, ignore.case = FALSE)

# Build two var sets (to divide the subsets)
vars_all    <- unique(c(core_vars, prs_vars))
vars_no_prs <- setdiff(core_vars, prs_vars)  

mod_all <- df_clean %>%
  dplyr::select(any_of(vars_all)) 

mod_all <- dfmL %>%
  mutate(
    functional_group = stringr::str_squish(functional_group),
    functional_group = factor(
      functional_group,
      levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
    )
  ) %>%
  dplyr::select(
    ch4_mgCH4_m2_d,
    soil_moisture_avg,
    x12cm_soil_temp,
    avg_par,
    prop_shrub, 
    prop_forb, 
    prop_graminoid,
    prop_Barren = prop_Barren,  # <- watch this naming! (its Caps, I missed it)
    gcc_mean,
    GPP_gC_m2_d,
    C_N_ratio,
    co2_gC_m2_d,
    functional_group,
    everything()  # keep extras if you want, harmless
  )



mod_no_prs <- df_clean %>%
  dplyr::select(any_of(vars_no_prs)) 

# Sanity check (number of columns is what would change) 
cat("n rows (ALL):    ", nrow(mod_all),    " | n cols:", ncol(mod_all),    "\n")
cat("n rows (NO_PRS): ", nrow(mod_no_prs), " | n cols:", ncol(mod_no_prs), "\n")

k_moist <- safe_k(mod_all$soil_moisture_avg)
k_t12   <- safe_k(mod_all$x12cm_soil_temp)
k_t5   <- safe_k(mod_all$x5cm_soil_temp)
k_b   <- safe_k(mod_all$veg_height_mean)
k_C   <- safe_k(mod_all$C_N_ratio)
k_bmd   <- safe_k(mod_all$bm_prop_dead)
k_gcc   <- safe_k(mod_all$gcc_mean)
k_gpp   <- safe_k(mod_all$GPP_gC_m2_d)
k_er    <- safe_k(mod_all$ER_gC_m2_d)
k_slope <- safe_k(mod_all$slope_aspect)

# PRS k’s 
k_map <- function(nm) if (nm %in% names(mod_all)) safe_k(mod_all[[nm]]) else 4L

form_all <- as.formula(
  paste(
    "ch4_mgCH4_m2_d ~",
    paste(c(
      sprintf("s(soil_moisture_avg, k=%d, bs='cr')", k_moist),
      sprintf("s(x12cm_soil_temp,   k=%d, bs='cr')", k_t12),
      sprintf("s(gcc_mean,          k=%d, bs='cr')", k_gcc),
      sprintf("s(GPP_gC_m2_d,       k=%d, bs='cr')", k_gpp),
      # PRS
      if ("total_n" %in% prs_vars) sprintf("s(total_n, k=%d, bs='cr')", k_map("total_n")) else NULL,
      if ("p"       %in% prs_vars) sprintf("s(p,       k=%d, bs='cr')", k_map("p")) else NULL,
      if ("cu"      %in% prs_vars) sprintf("s(cu,      k=%d, bs='cr')", k_map("cu")) else NULL
    ), collapse = " + ")
  )
)

gam_all <- gam(
  form_all, data = mod_all,
  method = "REML", select = TRUE
)

print(summary(gam_all))
print(mgcv::gam.check(gam_all))
# draw(gam_all)



plot_gamcheck_panel <- function(model) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # 2x2 layout, margins adjusted so nothing gets cut off
  par(mfrow = c(2, 2),
      mar = c(4, 4, 2, 1),   # bottom, left, top, right
      oma = c(2, 2, 2, 2))   # outer margins: bottom, left, top, right
  
  # Panel (a): Residuals vs linear predictor
  fit_vals <- fitted(model)
  res_pearson <- residuals(model, type = "pearson")
  
  plot(fit_vals, res_pearson,
       xlab = "Linear predictor", ylab = "Pearson residuals",
       pch = 19, cex = 0.6, col = rgb(0,0,0,0.6))
  abline(h = 0, col = "red", lwd = 1.2)
  rug(fit_vals, col = "gray")
  mtext("(a)", side = 3, line = 0.5, adj = 0, cex = 1.2)
  
  # Panel (b): QQ-plot of deviance residuals
  res_dev <- residuals(model, type = "deviance")
  qq <- qqnorm(res_dev, plot.it = FALSE)
  plot(qq$x, qq$y,
       xlab = "Theoretical quantiles", ylab = "Deviance residuals",
       pch = 19, cex = 0.6, col = rgb(0,0,0,0.6))
  qqline(res_dev, col = "red", lwd = 1.2)
  mtext("(b)", side = 3, line = 0.5, adj = 0, cex = 1.2)
  
  # Panel (c): Histogram of deviance 
  hist(res_dev,
       main = "", xlab = "Deviance residuals",
       col = "lightgray", border = "white", freq = FALSE)
  # lines(density(res_dev), col = "blue", lwd = 2)
  mtext("(c)", side = 3, line = 0.5, adj = 0, cex = 1.2)
  
  # Panel (d): Response vs fitted with semi-transparent points
  plot(fit_vals, model$y,
       xlab = "Fitted values", ylab = "Response",
       pch = 19, cex = 0.6, col = rgb(0,0,0,0.6))
  abline(0, 1, col = "red", lwd = 1.2)
  mtext("(d)", side = 3, line = 0.5, adj = 0, cex = 1.2)
  
  # Optional: add faint outer box for neatness
  box()
}


plot_gamcheck_panel(gam_all)



# gam_ch4_core_scat <- update(gam_all, family = mgcv::scat())
# print(summary(gam_ch4_core_scat))
# print(mgcv::gam.check(gam_ch4_core_scat))


set.seed(123)
k <- 5 # number of times
folds <- sample(rep(1:k, length.out = nrow(mod_all)))

cv_slim <- do.call(rbind, lapply(1:k, function(i) {
  # Split data
  train <- mod_all[folds != i, ]
  test  <- mod_all[folds == i, ]
  
  # Refit model on training data using the same formula
  m <- gam(formula(gam_all), data = train, method = "REML")
  
  # Predict on test data
  preds <- predict(m, newdata = test)
  
  # Collect observed and predicted
  data.frame(obs = test$ch4_mgCH4_m2_d, pred = preds)
}))

# Cross-validated R²
cv_R2_slim <- cor(cv_slim$obs, cv_slim$pred, use = "complete.obs")^2
cv_R2_slim


ggplot(cv_slim, aes(pred, obs)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Predicted CH₄ flux (mg C m⁻² d⁻¹)",
       y = "Observed CH₄ flux (mg C m⁻² d⁻¹)",
       title = paste("Cross-validation (5-fold): R² =", round(cv_R2_slim, 2))) +
  theme_minimal()


gam_all_light_slim <- gam(
  ch4_mgCH4_m2_d ~ 
    s(soil_moisture_avg, k = 4, bs = "cr") +
    s(x12cm_soil_temp, k = 4, bs = "cr") +
    s(avg_par, k = 3, bs = "cr") +
    s(GPP_gC_m2_d, k = 4, bs = "cr") +
    s(total_n, k = 4, bs = "cr") +
    s(p, k = 4, bs = "cr") +
    s(cu, k = 4, bs = "cr"),
  data = mod_all,
  method = "REML", select = TRUE
)


print(summary(gam_all_light_slim))
print(mgcv::gam.check(gam_all_light_slim))
# draw(gam_all_light_slim)


plot_gamcheck_panel(gam_all_light_slim)

# png("Minigam_check.png", width = 2000, height = 1600, res = 200)
# plot_gamcheck_panel(gam_all_light_slim)
# dev.off()

set.seed(123)
k <- 5 # number of times
folds <- sample(rep(1:k, length.out = nrow(mod_all)))

cv_slim <- do.call(rbind, lapply(1:k, function(i) {
  # Split data
  train <- mod_all[folds != i, ]
  test  <- mod_all[folds == i, ]
  
  # Refit model on training data using the same formula
  m <- gam(formula(gam_all_light_slim), data = train, method = "REML")
  
  # Predict on test data
  preds <- predict(m, newdata = test)
  
  # Collect observed and predicted
  data.frame(obs = test$ch4_mgCH4_m2_d, pred = preds)
}))

# Cross-validated R²
cv_R2_slim <- cor(cv_slim$obs, cv_slim$pred, use = "complete.obs")^2
cv_R2_slim

# this is sweet no terrible outliers or anything
# model graphs look good as well, some skewness but not too bad


ggplot(cv_slim, aes(pred, obs)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Predicted CH₄ flux (mg C m⁻² d⁻¹)",
       y = "Observed CH₄ flux (mg C m⁻² d⁻¹)",
       title = paste("Cross-validation (5-fold): R² =", round(cv_R2_slim, 2))) +
  theme_minimal()



# draw(gam_all_light_slim)

vars <- c(
  "soil_moisture_avg",
  "x12cm_soil_temp",
  "avg_par",
  "GPP_gC_m2_d",
  "total_n",
  "p",
  "cu"
)

results <- map_dfr(vars, function(v){
  
  form <- as.formula(
    paste0("ch4_mgCH4_m2_d ~ s(", v, ", k = 4, bs = 'cr')")
  )
  
  m <- gam(form, data = mod_all, method = "REML")
  
  sm <- summary(m)
  
  tibble(
    variable = v,
    edf = sm$s.table[1, "edf"],
    F = sm$s.table[1, "F"],
    p_value = sm$s.table[1, "p-value"],
    r_sq = sm$r.sq,
    dev_expl = sm$dev.expl
  )
})

results %>%
  arrange(desc(r_sq))



results %>% mutate(
  variable = recode(
    variable,
    soil_moisture_avg = "Soil moisture",
    x12cm_soil_temp = "Soil temperature (12 cm)",
    avg_par = "PAR",
    GPP_gC_m2_d = "GPP",
    total_n = "Total N",
    p = "P",
    cu = "Cu"
  )
)



pred_table <- results %>%
  arrange(desc(r_sq)) %>%
  mutate(
    Significance = case_when(
      p_value <= 0.001 ~ "P ≤ 0.001",
      p_value <= 0.01  ~ "P ≤ 0.01",
      p_value <= 0.05  ~ "P ≤ 0.05",
      TRUE             ~ "P > 0.05"
    ),
    variable = recode(
      variable,
      soil_moisture_avg = "Soil moisture",
      x12cm_soil_temp = "Soil temperature (12 cm)",
      avg_par = "PAR",
      GPP_gC_m2_d = "GPP",
      total_n = "Total N",
      p = "P",
      cu = "Cu"
    )
  ) %>%
  select(variable, Significance)

pred_table

# ---------- Light Model: no prs n~100 ----------


k_moist <- safe_k(mod_all$soil_moisture_avg)
k_t12   <- safe_k(mod_all$x12cm_soil_temp)
k_t5   <- safe_k(mod_all$x5cm_soil_temp)
k_b   <- safe_k(mod_all$veg_height_mean)
k_sa  <- safe_k(mod_all$slope_aspect)
k_sd <- safe_k((mod_all$soil_depth_avg))
k_bmd   <- safe_k(mod_all$bm_prop_dead)
k_gcc   <- safe_k(mod_all$gcc_mean)
k_gpp   <- safe_k(mod_all$GPP_gC_m2_d)
k_nee   <- safe_k(mod_all$co2_gC_m2_d)
k_er    <- safe_k(mod_all$ER_gC_m2_d)
k_slope <- safe_k(mod_all$slope_aspect)
k_C <- safe_k(mod_all$C_N_ratio)
k_par <- safe_k(mod_all$avg_par)

form_no_prs <- as.formula(
  paste(
    "ch4_mgCH4_m2_d ~",
    paste(c(
      sprintf("s(soil_moisture_avg, k=%d, bs='cr')", k_moist),
      sprintf("s(x12cm_soil_temp,   k=%d, bs='cr')", k_t12),
      "s(prop_shrub,      k=4, bs='cr')",
      "s(prop_forb,       k=4, bs='cr')",
      "s(prop_graminoid,  k=4, bs='cr')",
      "s(prop_Barren,     k=4, bs='cr')",
      sprintf("s(gcc_mean,          k=%d, bs='cr')", k_gcc),
      sprintf("s(C_N_ratio,          k=%d, bs='cr')", k_C),
      sprintf("s(avg_par,          k=%d, bs='cr')", k_par),
      sprintf("s(GPP_gC_m2_d,       k=%d, bs='cr')", k_gpp)
    ), collapse = " + ")
  )
)

gam_no_prs <- gam(
  form_no_prs, data = mod_no_prs,
  method = "REML", select = TRUE
)



print(summary(gam_no_prs))
print(mgcv::gam.check(gam_no_prs))
# draw(gam_no_prs)
# png("CNLandgam_check.png", width = 2000, height = 1600, res = 200)
plot_gamcheck_panel(gam_no_prs)
# dev.off()


set.seed(123)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(mod_no_prs)))

cv_no_prs <- do.call(rbind, lapply(1:k, function(i) {
  train <- mod_no_prs[folds != i, ]
  test  <- mod_no_prs[folds == i, ]
    m <- gam(formula(gam_no_prs), data = train, method = "REML")
    preds <- predict(m, newdata = test)
    data.frame(obs = test$ch4_mgCH4_m2_d, pred = preds)
}))

# Cross-validated R²
cv_R2_no_prs <- cor(cv_no_prs$obs, cv_no_prs$pred, use = "complete.obs")^2
cv_R2_no_prs


# draw(gam_no_prs)


vars <- c(
  "soil_moisture_avg",
  "x12cm_soil_temp",
  "prop_shrub",
  "prop_forb",
  "prop_graminoid",
  "prop_Barren",
  "gcc_mean",
  "C_N_ratio",
  "avg_par",
  "GPP_gC_m2_d"
)

results <- map_dfr(vars, function(v){
  
  form <- as.formula(
    paste0("ch4_mgCH4_m2_d ~ s(", v, ", k = 4, bs = 'cr')")
  )
  
  m <- gam(form, data = mod_no_prs, method = "REML")
  
  sm <- summary(m)
  
  tibble(
    variable = v,
    edf = sm$s.table[1, "edf"],
    F = sm$s.table[1, "F"],
    p_value = sm$s.table[1, "p-value"],
    r_sq = sm$r.sq,
    dev_expl = sm$dev.expl
  )
})

results %>%
  arrange(desc(r_sq))



results %>% mutate(
  variable = recode(
    variable,
    soil_moisture_avg = "Soil moisture",
    x12cm_soil_temp = "Soil temperature (12 cm)",
    avg_par = "PAR",
    GPP_gC_m2_d = "GPP",
    total_n = "Total N",
    p = "P",
    cu = "Cu"
  )
)



pred_table <- results %>%
  arrange(desc(r_sq)) %>%
  mutate(
    Significance = case_when(
      p_value <= 0.001 ~ "P ≤ 0.001",
      p_value <= 0.01  ~ "P ≤ 0.01",
      p_value <= 0.05  ~ "P ≤ 0.05",
      TRUE             ~ "P > 0.05"
    ),
    variable = recode(
      variable,
      soil_moisture_avg = "Soil moisture",
      x12cm_soil_temp   = "Soil temperature (12 cm)",
      prop_shrub        = "Shrub cover",
      prop_forb         = "Forb cover",
      prop_graminoid    = "Graminoid cover",
      prop_Barren       = "Bare ground",
      gcc_mean          = "GCC",
      C_N_ratio         = "C:N ratio",
      avg_par           = "PAR",
      GPP_gC_m2_d       = "GPP"
    )
  ) %>%
  select(variable, Significance)

pred_table

# ---------- Dark Model: PRS n~40 ------------------



glimpse(df_master)

dfD <- df_master %>% 
  dplyr::filter(light_dark == "D")

dfmD <- dfD %>%
  mutate(
    prop_shrub_sum = rowSums(
      cbind(prop_Erect.Shrub, prop_Prostrate.Shrub),
      na.rm = TRUE
    ),
    prop_shrub_sum = ifelse(
      is.na(prop_Erect.Shrub) & is.na(prop_Prostrate.Shrub),
      NA,
      prop_shrub_sum
    ),
    prop_shrub = dplyr::coalesce(
      prop_shrub_sum,
      cover_shrub / 100
    ),
    prop_graminoid = dplyr::coalesce(prop_Graminoid, cover_graminoid/100),
    prop_forb      = dplyr::coalesce(prop_Forb,      cover_forb/100),
    
    even3          = dplyr::coalesce(evenness),
    rich3          = dplyr::coalesce(as.numeric(richness))
  ) %>%
  select(-prop_shrub_sum)


df_clean_D <- dfmD %>%
  rename(
    prop_prostrate_shrub =  prop_Prostrate.Shrub,
    prop_erect_shrub     =  prop_Erect.Shrub,
    cover_prostrate_shrub=  cover_Prostrate.Shrub,
    cover_erect_shrub    =  cover_Erect.Shrub
  )



mod_all_D <- df_clean_D %>%
  dplyr::select(any_of(vars_all)) 

mod_all_D <- dfmD %>%
  mutate(
    functional_group = stringr::str_squish(functional_group),
    functional_group = factor(
      functional_group,
      levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
    )
  ) %>%
  select(
    ch4_mgCH4_m2_d,
    soil_moisture_avg,
    x12cm_soil_temp,
    prop_shrub, 
    prop_forb, 
    prop_graminoid,
    prop_Barren = prop_Barren,  # <- watch this naming! (its Caps, I missed it)
    gcc_mean,
    GPP_gC_m2_d,
    slope_aspect,
    functional_group,
    everything()  # keep extras if you want, harmless
  )


mod_all_D <- mod_all_D %>%
  filter(ER_gC_m2_d >= 0)

mod_no_prs_D <- df_clean_D %>%
  dplyr::select(any_of(vars_no_prs)) 

mod_no_prs_D <- mod_no_prs_D %>%
  filter(ER_gC_m2_d >= 0)

# Sanity check (number of columns is what would change) 
cat("n rows (ALL):    ", nrow(mod_all_D),    " | n cols:", ncol(mod_all_D),    "\n")
cat("n rows (NO_PRS): ", nrow(mod_no_prs_D), " | n cols:", ncol(mod_no_prs_D), "\n")

k_moist <- safe_k(mod_all_D$soil_moisture_avg)
k_t12   <- safe_k(mod_all_D$x12cm_soil_temp)
k_t5   <- safe_k(mod_all_D$x5cm_soil_temp)
k_b   <- safe_k(mod_all_D$veg_height_mean)
k_C   <- safe_k(mod_all_D$C_N_ratio)
k_bmd   <- safe_k(mod_all_D$bm_prop_dead)
k_gcc   <- safe_k(mod_all_D$gcc_mean)
k_gpp   <- safe_k(mod_all_D$GPP_gC_m2_d)
k_er    <- safe_k(mod_all_D$ER_gC_m2_d)
k_slope <- safe_k(mod_all_D$slope_aspect)

# PRS k’s 
k_map <- function(nm) if (nm %in% names(mod_all_D)) safe_k(mod_all_D[[nm]]) else 4L

form_all_D <- as.formula(
  paste(
    "ch4_mgCH4_m2_d ~",
    paste(c(
      sprintf("s(soil_moisture_avg, k=%d, bs='cr')", k_moist),
      sprintf("s(x12cm_soil_temp,   k=%d, bs='cr')", k_t12),
      sprintf("s(ER_gC_m2_d,        k=%d, bs='cr')", k_er),
      sprintf("s(C_N_ratio,        k=%d, bs='cr')", k_C),
      # PRS
      if ("total_n" %in% prs_vars) sprintf("s(total_n, k=%d, bs='cr')", k_map("total_n")) else NULL,
      if ("p"       %in% prs_vars) sprintf("s(p,       k=%d, bs='cr')", k_map("p")) else NULL,
      if ("cu"      %in% prs_vars) sprintf("s(cu,      k=%d, bs='cr')", k_map("cu")) else NULL
    ), collapse = " + ")
  )
)

gam_all_D <- gam(
  form_all_D, data = mod_all_D,
  method = "REML", select = TRUE
)

print(summary(gam_all_D))
print(mgcv::gam.check(gam_all_D))


# png("Dark_Minigam_check.png", width = 2000, height = 1600, res = 200)
plot_gamcheck_panel(gam_all_D)
# dev.off()
# draw(gam_all_D)
set.seed(123)
k <- 5 # number of times
folds <- sample(rep(1:k, length.out = nrow(mod_all_D)))

cv_slim_D <- do.call(rbind, lapply(1:k, function(i) {
  # Split data
  train <- mod_all_D[folds != i, ]
  test  <- mod_all_D[folds == i, ]
  
  # Refit model on training data using the same formula
  m <- gam(formula(gam_all_D), data = train, method = "REML")
  
  # Predict on test data
  preds <- predict(m, newdata = test)
  
  # Collect observed and predicted
  data.frame(obs = test$ch4_mgCH4_m2_d, pred = preds)
}))

# Cross-validated R²
cv_R2_slim_D <- cor(cv_slim_D$obs, cv_slim_D$pred, use = "complete.obs")^2
cv_R2_slim_D
gam_all_light_slim_D <- gam(
  ch4_mgCH4_m2_d ~
    s(x12cm_soil_temp, k = 4, bs = "cr") +
    s(ER_gC_m2_d, k = 6, bs = "cr") +
    s(total_n, k = 4, bs = "cr") +
    s(p, k = 4, bs = "cr") +
    s(cu, k = 4, bs = "cr") +
    s(soil_moisture_avg, k = 4, bs = "cr"),
  data = mod_all_D,
  method = "REML", select = TRUE
)




print(summary(gam_all_light_slim_D))
print(mgcv::gam.check(gam_all_light_slim_D))
# draw(gam_all_light_slim_D)

# png("Dark_Minigam_check.png", width = 2000, height = 1600, res = 200)
plot_gamcheck_panel(gam_all_light_slim_D)
# dev.off()

set.seed(123)
k <- 5 # number of times
folds <- sample(rep(1:k, length.out = nrow(mod_all_D)))

cv_slim_D <- do.call(rbind, lapply(1:k, function(i) {
  # Split data
  train <- mod_all_D[folds != i, ]
  test  <- mod_all_D[folds == i, ]
  
  # Refit model on training data using the same formula
  m <- gam(formula(gam_all_light_slim_D), data = train, method = "REML")
  
  # Predict on test data
  preds <- predict(m, newdata = test)
  
  # Collect observed and predicted
  data.frame(obs = test$ch4_mgCH4_m2_d, pred = preds)
}))

# Cross-validated R²
cv_R2_slim_D <- cor(cv_slim_D$obs, cv_slim_D$pred, use = "complete.obs")^2
cv_R2_slim_D


ggplot(cv_slim_D, aes(pred, obs)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Predicted CH₄ flux (mg C m⁻² d⁻¹)",
       y = "Observed CH₄ flux (mg C m⁻² d⁻¹)",
       title = paste("Cross-validation (5-fold): R² =", round(cv_R2_slim_D, 2))) +
  theme_minimal()


# ---------- Dark Model: no prs n~100 ----------

k_moist <- safe_k(mod_all_D$soil_moisture_avg)
k_t12   <- safe_k(mod_all_D$x12cm_soil_temp)
k_t5   <- safe_k(mod_all_D$x5cm_soil_temp)
k_b   <- safe_k(mod_all_D$veg_height_mean)
k_C   <- safe_k(mod_all_D$C_N_ratio)
k_sa  <- safe_k(mod_all_D$slope_aspect)
k_sd <- safe_k((mod_all_D$soil_depth_avg))
k_bmd   <- safe_k(mod_all_D$bm_prop_dead)
k_gcc   <- safe_k(mod_all_D$gcc_mean)
k_gpp   <- safe_k(mod_all_D$GPP_gC_m2_d)
k_er    <- safe_k(mod_all_D$ER_gC_m2_d)
slope_aspect <- safe_k(mod_all_D$slope_aspect)
set.seed(123)


form_no_prs_D <- as.formula(
  paste(
    "ch4_mgCH4_m2_d ~",
    paste(c(
      sprintf("s(soil_moisture_avg, k=4, bs='cr')", k_moist),
      sprintf("s(x12cm_soil_temp,   k=4, bs='cr')", k_t12),
      "s(prop_shrub,      k=3, bs='cr')",
      "s(prop_forb,       k=3, bs='cr')",
      "s(prop_graminoid,  k=3, bs='cr')",
      "s(prop_Barren,     k=3, bs='cr')",
      sprintf("s(C_N_ratio,        k=%d, bs='cr')", k_C),
      sprintf("s(gcc_mean,          k=4, bs='cr')", k_gcc),
      "s(ER_gC_m2_d,     k=6, bs='cr')"

      # sprintf("s(ER_gC_m2_d,        k=%d, bs='cr')", k_er)
    ), collapse = " + ")
  )
)






gam_no_prs_D <- gam(
  form_no_prs_D, data = mod_no_prs_D,
  method = "REML", select = TRUE
)



print(summary(gam_no_prs_D))
print(mgcv::gam.check(gam_no_prs_D))

# png("CNDark_Landgam_check.png", width = 2000, height = 1600, res = 200)
# plot_gamcheck_panel(gam_no_prs_D)
# dev.off()

set.seed(123)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(mod_no_prs_D)))

cv_no_prs_D <- do.call(rbind, lapply(1:k, function(i) {
  train <- mod_no_prs_D[folds != i, ]
  test  <- mod_no_prs_D[folds == i, ]
  m <- gam(formula(gam_no_prs_D), data = train, method = "REML")
  preds <- predict(m, newdata = test)
  data.frame(obs = test$ch4_mgCH4_m2_d, pred = preds)
}))

# Cross-validated R²
cv_R2_no_prs_D <- cor(cv_no_prs_D$obs, cv_no_prs_D$pred, use = "complete.obs")^2
cv_R2_no_prs_D



# draw(gam_no_prs_D)


# ----------- prediction underestimation-----------

df_pred_no_prs_D <- df_master %>%
  slice(1:nrow(gam_no_prs_D$model)) %>%  # only rows used in model
  mutate(pred_no_prs_D = fitted(gam_no_prs_D))

df_pred_no_prs <- df_master %>%
  slice(1:nrow(gam_no_prs$model)) %>%  # only rows used in model
  mutate(pred_no_prs = fitted(gam_no_prs))

df_gam_all_light_slim_D <- df_master %>%
  slice(1:nrow(gam_all_light_slim_D$model)) %>%  # only rows used in model
  mutate(pred_gam_all_light_slim_D = fitted(gam_all_light_slim_D))

df_gam_all_light_slim <- df_master %>%
  slice(1:nrow(gam_all_light_slim$model)) %>%  # only rows used in model
  mutate(pred_gam_all_light_slim = fitted(gam_all_light_slim))




df_summary_dark <- df_pred_no_prs_D %>%
  group_by(functional_group) %>%
  summarise(mean_dark = mean(pred_no_prs_D, na.rm = TRUE))

df_summary_light <- df_pred_no_prs %>%
  group_by(functional_group) %>%
  summarise(mean_light = mean(pred_no_prs, na.rm = TRUE))

df_summary_dark_N <- df_gam_all_light_slim_D %>%
  group_by(functional_group) %>%
  summarise(mean_dark_N = mean(pred_gam_all_light_slim_D, na.rm = TRUE))

df_summary_light_N <- df_gam_all_light_slim %>%
  group_by(functional_group) %>%
  summarise(mean_light_N = mean(pred_gam_all_light_slim, na.rm = TRUE))

# Check which functional groups are in the dataframe
unique(df_gam_all_light_slim_D$functional_group)


df_all <- df_summary_dark %>%
  full_join(df_summary_light, by = "functional_group") %>%
  full_join(df_summary_dark_N, by = "functional_group") %>%
  full_join(df_summary_light_N, by = "functional_group")



df_all <- df_all %>%
  rowwise() %>%
  mutate(
    mean_dark_all  = mean(c(mean_dark, mean_dark_N), na.rm = TRUE),
    mean_light_all = mean(c(mean_light, mean_light_N), na.rm = TRUE),
    underestimation_pct = (mean_light_all - mean_dark_all) / abs(mean_light_all) * 100
  ) %>%
  ungroup()

print(df_all)

range_underestimation <- range(df_all$underestimation_pct, na.rm = TRUE)

range_underestimation

# ----------- partial effefcts ---------

AIC(gam_no_prs_D, gam_no_prs, gam_all_light_slim, gam_all_light_slim_D, gam_all_D)

sm_all <- bind_rows(
  smooth_estimates(gam_no_prs_D)        %>% mutate(model = "No PRS – Dark"),
  smooth_estimates(gam_no_prs)          %>% mutate(model = "No PRS – Light"),
  smooth_estimates(gam_all_D)%>% mutate(model = "All – Dark"),
  smooth_estimates(gam_all_light_slim)  %>% mutate(model = "All – Light")
)

sm_long <- sm_all %>%
  pivot_longer(
    cols = c(
      soil_moisture_avg, x12cm_soil_temp, prop_shrub,
      prop_forb, prop_graminoid, prop_Barren, avg_par,
      C_N_ratio, gcc_mean, ER_gC_m2_d,GPP_gC_m2_d,
      total_n, 
      p,
      cu
    ),
    names_to = "xvar",
    values_to = "x"
  ) %>%
  filter(
    (.smooth == "s(soil_moisture_avg)" & xvar == "soil_moisture_avg")   |
      (.smooth == "s(x12cm_soil_temp)"   & xvar == "x12cm_soil_temp")   |
      (.smooth == "s(prop_shrub)"        & xvar == "prop_shrub")        |
      (.smooth == "s(avg_par)"           & xvar == "avg_par")           |
      (.smooth == "s(prop_forb)"         & xvar == "prop_forb")         |
      (.smooth == "s(prop_graminoid)"    & xvar == "prop_graminoid")    |
      (.smooth == "s(prop_Barren)"       & xvar == "prop_Barren")       |
      (.smooth == "s(C_N_ratio)"         & xvar == "C_N_ratio")         |
      (.smooth == "s(gcc_mean)"          & xvar == "gcc_mean")          |
      (.smooth == "s(GPP_gC_m2_d)"       & xvar == "GPP_gC_m2_d")       |
      (.smooth == "s(ER_gC_m2_d)"        & xvar == "ER_gC_m2_d")        |
      (.smooth == "s(total_n)"           & xvar == "total_n")           |
      (.smooth == "s(p)"                 & xvar == "p")                 |
      (.smooth == "s(cu)"                & xvar == "cu")
  )%>%
  filter(!is.na(model))

sm_light <- sm_long %>%
  filter(model %in% c("No PRS – Light", "All – Light")) %>%
  mutate(
    model = factor(
      model,
      levels = c("All – Light", "No PRS – Light")
    )
  )


ch4_lab <- expression("Methane flux ("*mg~CH[4]~m^{-2}~day^{-1}*")")


sm_slim <- sm_long %>%
  filter(model %in% c("All – Light", "All – Dark")) %>%
  mutate(
    model = factor(
      model,
      levels = c("All – Light", "All – Dark")
    )
  )
sm_noprs <- sm_long %>%
  filter(model %in% c("No PRS – Light", "No PRS – Dark")) %>%
  mutate(
    model = factor(
      model,
      levels = c("No PRS – Light", "No PRS – Dark")
    )
  )


dfmL_long <- dfmL %>%
  pivot_longer(
    cols = c(
      soil_moisture_avg, x12cm_soil_temp, prop_shrub,
      prop_forb, prop_graminoid, prop_Barren,avg_par,
      C_N_ratio, gcc_mean, ER_gC_m2_d, GPP_gC_m2_d, #co2_gC_m2_d,
      slope_aspect, total_n, veg_height_mean,
      p, s, mn, cu
    ),
    names_to = "xvar",
    values_to = "x"
  )


dfmL_long_sub2 <- dfmL_long %>%
  select(plot_key, xvar, x, y = ch4_mgCH4_m2_d, functional_group) %>%
  rename(x_obs = x, y_obs = y)


# Join on xvar
sm_light_full <- sm_light %>%
  left_join(dfmL_long_sub2, by = "xvar")
unique(dfmL_long_sub2$xvar)



sm_light_single <- sm_light_full %>% 
  filter(model == "No PRS – Light" | is.na(model)) 

p_light_single <- ggplot() +
  # raw observations colored by functional group
  geom_point(data = sm_light_single %>% filter(!is.na(y_obs)),
             aes(x = x_obs, y = y_obs, color = functional_group),
             alpha = 0.1) +
  scale_color_manual(values = grp_cols) +
  # smooth lines
  geom_line(data = sm_light_single, aes(x = x, y = .estimate, group = model)) +
  geom_ribbon(data = sm_light_single, 
              aes(x = x, ymin = .estimate - 2*.se, ymax = .estimate + 2*.se, group = model), 
              alpha = 0.6) +
  facet_wrap(. ~ .smooth, scales = "free_x") +  # only facet by variable
  ylab(ch4_lab) + 
  theme_bw()



p_light_single

prs_plots <- dfmL %>%
  filter(!is.na(p) | !is.na(s) | !is.na(mn) | !is.na(cu)) %>%
  distinct(plot_key) %>%
  pull(plot_key)

single_model <- "All – Light"
single_model <- "No PRS – Light"



sm_light_single <- sm_light_full %>%
  filter(
    (model == "All – Light" & plot_key %in% prs_plots) |
      model == "No PRS – Light" |
      is.na(model)
  )




sm_light_single %>%
  distinct(model, plot_key) %>%
  count(model)


sm_light_single <- sm_light_full %>%
  filter(model == single_model | is.na(model))

sm_light_single %>%
  count(model, xvar)



sm_light_single2 <- sm_light_single %>%
  mutate(
    xvar = gsub("^s\\(|\\)$", "", .smooth)
  )


dfmL_long_sub2_filtered <- dfmL_long_sub2 %>%
  filter(if (single_model == "All – Light") plot_key %in% prs_plots else TRUE)

sm_light_plot <- sm_light_single2 %>%
  inner_join(dfmL_long_sub2_filtered, by = "xvar")

sm_light_clean <- sm_light_plot %>%
  transmute(
    .smooth, .type, .by, .estimate, .se, model, xvar, x,
    plot_key = plot_key.y,
    x_obs   = x_obs.y,
    y_obs   = y_obs.y,
    functional_group = functional_group.y
  ) %>%
  distinct()

anyDuplicated(sm_light_clean[c("plot_key","xvar","x_obs","y_obs")])

p_light_single <- ggplot(sm_light_clean) +
  
  geom_point(
    aes(x = x_obs, y = y_obs, color = functional_group),
    shape = 20,
    alpha = 0.05,
    size  = 1
  ) +
  
  scale_color_manual(
    values = grp_cols,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  
  geom_line(
    aes(x = x, y = .estimate, group = model)
  ) +
  
  geom_ribbon(
    aes(
      x = x,
      ymin = .estimate - 3 * .se,
      ymax = .estimate + 3 * .se,
      group = model
    ),
    alpha = 0.4
  ) +
  
  facet_wrap(~ .smooth, scales = "free") +
  ylab(ch4_lab) +
  theme_bw()

p_light_single

# ggsave("Light_Landsacape_partial_effects.png", p_light_single, width = 10, height = 6, dpi = 300)

# Repeat for dark
dfmD_long <- dfmD %>%
  pivot_longer(
    cols = c(
      soil_moisture_avg, x12cm_soil_temp, prop_shrub,
      prop_forb, prop_graminoid, prop_Barren,
      C_N_ratio, gcc_mean, ER_gC_m2_d, GPP_gC_m2_d, #co2_gC_m2_d,
      slope_aspect, total_n, veg_height_mean,
      p, s, mn, cu
    ),
    names_to = "xvar",
    values_to = "x"
  )

dfmD_long_sub2 <- dfmD_long %>%
  select(plot_key, xvar, x, y = ch4_mgCH4_m2_d, functional_group) %>%
  rename(x_obs = x, y_obs = y)

sm_dark_full <- sm_dark %>%
  left_join(dfmD_long_sub2, by = "xvar")

prs_plots <- dfmD %>%
  filter(!is.na(p) | !is.na(s) | !is.na(mn) | !is.na(cu)) %>%
  distinct(plot_key) %>%
  pull(plot_key)

single_model <- "All – Dark"
# single_model <- "No PRS – Dark"
sm_dark_single <- sm_dark_full %>%
  filter(
    (model == "All – Dark" & plot_key %in% prs_plots) |
      model == "No PRS – Dark" |
      is.na(model)
  )

sm_dark_single %>%
  distinct(model, plot_key) %>%
  count(model)

sm_dark_single <- sm_dark_full %>%
  filter(model == single_model | is.na(model))

sm_dark_single %>%
  count(model, xvar)

sm_dark_single2 <- sm_dark_single %>%
  mutate(
    xvar = gsub("^s\\(|\\)$", "", .smooth)
  )

bad_plots <- dfmD_long_sub2 %>%
  filter(xvar == "ER_gC_m2_d", x_obs < 0) %>%
  distinct(plot_key) %>%
  pull(plot_key)


dfmD_long_sub2_filtered <- dfmD_long_sub2 %>%
  filter(
    if (single_model == "All – Dark") plot_key %in% prs_plots else TRUE
  ) %>%
  group_by(plot_key) %>%
  filter(
    !any(xvar == "ER_gC_m2_d" & x_obs < 0, na.rm = TRUE)
  ) %>%
  ungroup()

sm_dark_plot <- sm_dark_single2 %>%
  inner_join(dfmD_long_sub2_filtered, by = "xvar")


sm_dark_clean <- sm_dark_plot %>%
  filter(!plot_key.x %in% bad_plots) %>%   
  transmute(
    .smooth, .type, .by, .estimate, .se, model, xvar, x,
    plot_key = plot_key.y,
    x_obs   = x_obs.y,
    y_obs   = y_obs.y,
    functional_group = functional_group.y
  ) %>%
  distinct()


sm_dark_clean %>%
  semi_join(
    dfmD_long_sub2 %>% 
      filter(xvar == "ER_gC_m2_d", x_obs < 0),
    by = "plot_key"
  )

anyDuplicated(sm_dark_clean[c("plot_key","xvar","x_obs","y_obs")])

p_dark_single <- ggplot(sm_dark_clean) +
  geom_point(
    aes(x = x_obs, y = y_obs, color = functional_group),
    shape = 20,
    alpha = 0.05,
    size  = 1
  ) +
  scale_color_manual(
    values = grp_cols,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  geom_line(
    aes(x = x, y = .estimate, group = model)
  ) +
  geom_ribbon(
    aes(
      x = x,
      ymin = .estimate - 3 * .se,
      ymax = .estimate + 3 * .se,
      group = model
    ),
    alpha = 0.4
  ) +
  facet_wrap(~ .smooth, scales = "free") +
  ylab(ch4_lab) +
  theme_bw()

p_dark_single


# ggsave("Dark_Nutrient_partial_effects.png", p_dark_single, width = 10, height = 6, dpi = 300)

# -------------------- Partial Figure (Paper) -------------------------

grp_cols_light <- desaturate(lighten(grp_cols, 0.35), 0.3)

sm_all <- bind_rows(
  smooth_estimates(gam_no_prs_D, n = 400)         %>% mutate(model = "No PRS – Dark"),
  smooth_estimates(gam_no_prs, n = 400)           %>% mutate(model = "No PRS – Light"),
  smooth_estimates(gam_all_light_slim_D, n = 400) %>% mutate(model = "All – Dark"),
  smooth_estimates(gam_all_light_slim, n = 400)   %>% mutate(model = "All – Light")
)

df_obs <- bind_rows(
  dfmL_long_sub2 %>% mutate(measurement = "Light"),
  dfmD_long_sub2 %>% mutate(measurement = "Dark")
)

sm_long2 <- sm_all %>%
  mutate(xvar = gsub("^s\\(|\\)$", "", .smooth))

unique(sm_long2$xvar)
unique((df_obs$xvar))

df_obs_clean <- df_obs %>%
  distinct(plot_key, xvar, x_obs, y_obs, functional_group)

unique(df_obs_clean$xvar)

sm_big_clean <- sm_long2 %>% 
  left_join(df_obs_clean, by = "xvar") %>% 
  distinct( model, .smooth,
            .estimate, .se, plot_key, x_obs, y_obs, functional_group )

sm_big_clean <- sm_big_clean %>%
  group_by(model, .smooth) %>%
  mutate(
    n_obs = n_distinct(plot_key),
    se_scaled = .se * sqrt(mean(n_obs) / n_obs)
  ) %>%
  ungroup()
unique(sm_big_clean$.smooth)


sm_big_clean <- sm_big_clean %>%
  mutate(
    model = recode(
      model,
      "All – Light"     = "Nutrient-level GAM (Light)",
      "All – Dark"      = "Nutrient-level GAM (Dark)",
      "No PRS – Light"  = "Landscape-level GAM (Light)",
      "No PRS – Dark"   = "Landscape-level GAM (Dark)"
    )
  )
unique(sm_big_clean$.smooth)

model_cols <- c(
  "Nutrient-level GAM (Light)"    = "#F4C26B",  # yellow-leaning orange
  "Nutrient-level GAM (Dark)"     = "#8FA3C7",  # deep, rich orange
  "Landscape-level GAM (Light)"   = "#E07A2F",  # greenish / bluish purple
  "Landscape-level GAM (Dark)"    = "#5B4B8A"   # deepish "grounded" purple
)

sm_big_clean <- sm_big_clean %>%
  group_by(model, .smooth) %>%
  mutate(
    x_min = min(x_obs, na.rm = TRUE),
    x_max = max(x_obs, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(
    is.na(x_obs) | (x_obs >= x_min & x_obs <= x_max)
  )


# (g C m⁻² d⁻¹)
smooth_labs <- c(
  "s(GPP_gC_m2_d)" = "GPP (g C m⁻² d⁻¹)",
  "s(ER_gC_m2_d)"  = "Ecosystem Respiration (g C m⁻² d⁻¹)",
  "s(avg_par)"  = "PAR",
  
  "s(veg_height_mean)" = "Vegetation Height (cm)",
  "s(x12cm_soil_temp)"  = "Soil Temperature (12 cm)",
  "s(soil_moisture_avg)" = "Soil Moisture (%)",
  "s(gcc_mean)" = "GCC",
  "s(C_N_ratio)" = "Soil C/N",
 
  "s(s)"  = "S (µg 10 cm⁻²)",
  "s(cu)" = "Cu (µg 10 cm⁻²)",
  "s(mn)" = "Mn (µg 10 cm⁻²)",
  "s(total_n)" = "Total N (µg 10 cm⁻²)",
  "s(p)" = "P (µg 10 cm⁻²)",
  
  "s(prop_graminoid)" = "Graminoid Cover",
  "s(prop_forb)"      = "Forb Cover",
  "s(prop_shrub)"     = "Shrub Cover",
  "s(prop_Barren)"    = "Barren Cover"
)

sm_big_clean <- sm_big_clean %>%
  mutate(
    smooth_label = smooth_labs[.smooth]
  )

sm_big_clean <- sm_big_clean %>%
  mutate(
    smooth_label = factor(
      smooth_label,
      levels = c(
        "GPP (g C m⁻² d⁻¹)",
        "Ecosystem Respiration (g C m⁻² d⁻¹)",
        "Soil Temperature (12 cm)",
        "Soil Moisture (%)",
        "PAR",
        "Cu (µg 10 cm⁻²)",
        "P (µg 10 cm⁻²)",
        "Total N (µg 10 cm⁻²)",
        "Soil C/N",
        "S (µg 10 cm⁻²)",
        "Mn (µg 10 cm⁻²)",
        "Vegetation Height (cm)",
        "GCC",
        "Barren Cover",
        "Graminoid Cover",
        "Forb Cover",
        "Shrub Cover"
      )
    )
  )


p_big <- ggplot(sm_big_clean) +
  
  ## ---- Raw observations (vegetation) ----
geom_point(
  aes(x = x_obs, y = y_obs, color = functional_group),
  color = "grey95",
  alpha = 0.025,
  size = 1
) +
  # scale_color_manual(
  #   values = grp_cols_light,
  #   name = "Vegetation",
  #   guide = guide_legend(override.aes = list(alpha = 1))
  # ) +
  
  new_scale_color() +
  new_scale_fill() +
  
  ## ---- Smooths (models) ----
geom_line(
  aes(x = x, y = .estimate, color = model, group = model),
  linewidth = 1.1
) +
  
  geom_ribbon(
    aes(
      x = x,
      ymin = .estimate - 2 * se_scaled,
      ymax = .estimate + 2 * se_scaled,
      fill = model,
      group = model
    ),
    alpha = 0.3,
    linewidth = 0
  ) +
  
  scale_color_manual(values = model_cols, name = "Model") +
  scale_fill_manual(values = model_cols, name = "Model") +
  
  scale_y_continuous(
    limits = c(-5, 5),
    breaks = c(-5, 0, 5)
  ) +

  
  facet_wrap(
    ~ smooth_label,
    scales = "free_x"
  ) +
 
  ylab(ch4_lab) +
  xlab(NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  )

# ggsave(
#   filename = "GPP_CH4_partial_effects_all_models.png",
#   plot     = p_big,
#   width    = 14,
#   height   = 10,
#   dpi      = 300
# )

df_big   <- sm_big_clean %>% filter(smooth_label %in% c("GPP (g C m⁻² d⁻¹)", "Ecosystem Respiration (g C m⁻² d⁻¹)"))
df_small <- sm_big_clean %>% filter(!smooth_label %in% c("GPP (g C m⁻² d⁻¹)", "Ecosystem Respiration (g C m⁻² d⁻¹)"))

# df_big <- df_big %>%
#   group_by(smooth_label, model) %>%
#   mutate(
#     x_min = min(xvar),
#     x_max = max(xvar),
#     x = seq(min(xvar), max(xvar), length.out = n())
#   ) %>%
#   ungroup()


df_big <- df_big %>%
  group_by(smooth_label, model) %>%
  mutate(
    x = seq(x_min[1], x_max[1], length.out = n())
  ) %>%
  ungroup()

# df_small <- df_small %>%
#   group_by(smooth_label, model) %>%
#   mutate(
#     x_min = min(xvar),
#     x_max = max(xvar),
#     x = seq(min(xvar), max(xvar), length.out = n())
#   ) %>%
#   ungroup()


df_small <- df_small %>%
  group_by(smooth_label, model) %>%
  mutate(
    x = seq(x_min[1], x_max[1], length.out = n())
  ) %>%
  ungroup()

label_df <- sm_big_clean %>%
  distinct(smooth_label, model) %>%
  group_by(smooth_label) %>%
  mutate(model_row = row_number()) %>%
  ungroup()

label_big   <- label_df %>% 
  filter(smooth_label %in% c("GPP (g C m⁻² d⁻¹)", "Ecosystem Respiration (g C m⁻² d⁻¹)"))
label_small <- label_df %>% 
  filter(!smooth_label %in% c("GPP (g C m⁻² d⁻¹)", "Ecosystem Respiration (g C m⁻² d⁻¹)"))

coord_cartesian(clip = "off")

label_big <- label_big %>%
  left_join(
    df_big %>% group_by(smooth_label) %>% summarise(max_x = max(x, na.rm = TRUE)),
    by = "smooth_label"
  )
# unique(label_big$smooth_label)
# label_big <- label_big %>%
#   dplyr::filter(
#     !(smooth_label == "Ecosystem Respiration (g C m⁻² d⁻¹)" & max_x > 10)
#   )


label_small <- label_small %>%
  left_join(
    df_small %>% group_by(smooth_label) %>% summarise(max_x = max(x, na.rm = TRUE)),
    by = "smooth_label"
  )

df_big <- df_big  %>% 
  dplyr::arrange(model, smooth_label, x)

  p_big <- ggplot(df_big, aes(x = x, y = .estimate, color = model, group = model)) +
    geom_point(
      aes(x = x_obs, y = y_obs, color = functional_group),
      color = "grey76",
      alpha = 0.025,
      size = 0.5
    ) +
    geom_segment(
      data = label_big,
      aes(
        x    = max_x * 0.94,
        xend = max_x * 0.98,
        y    = 4.6 - 0.4 * (model_row - 1),
        yend = 4.6 - 0.4 * (model_row - 1),
        color = model
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    
    geom_ribbon(
      aes(
        ymin = .estimate - 2 * se_scaled,
        ymax = .estimate + 2 * se_scaled,
        fill = model,
        group = model
      ),
      alpha = 0.3
    ) +
  geom_line(linewidth = 1.1) +

  
  facet_wrap(~ smooth_label, ncol = 2, scales = "free_x") +
  
  scale_y_continuous(
    limits = c(-5, 5),
    breaks = c(-5, 0, 5)
    # expand = expansion(mult = c(0.05, 0.18))
  ) +
  
  scale_color_manual(values = model_cols, name = "Model") +
  scale_fill_manual(values = model_cols, name = "Model") +
  
  coord_cartesian(clip = "off") +
  
  ylab(ch4_lab) +
  xlab(NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.title = element_text(size =15),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 14), face = "bold")


  
label_small <- label_small %>%
  group_by(smooth_label) %>%
  mutate(
    # segment y starting near top of panel
    seg_y = 3 - 0.25 * (model_row - 1),  # 3 = top of y-scale for small panels
    seg_yend = seg_y
  ) %>%
  ungroup()
  
  
  
  p_small <- ggplot(df_small, aes(x = x, y = .estimate, color = model, group = model)) +
    geom_point(
      aes(x = x_obs, y = y_obs, color = functional_group),
      color = "grey76",
      alpha = 0.025,
      size = 0.5
    ) +
    geom_segment(
      data = label_small,
      aes(
        x    = max_x * 0.95,
        xend = max_x * 0.98,
        y    = seg_y,
        yend = seg_yend,
        color = model
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
  # scale_color_manual(
  #   values = grp_cols_light,
  #   name = "Vegetation",
  #   guide = guide_legend(override.aes = list(alpha = 1))
  # ) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(
    aes(ymin = .estimate - 2 * se_scaled,
        ymax = .estimate + 2 * se_scaled,
        fill = model),
    alpha = 0.3
  ) +
  
  coord_cartesian(clip = "off") +
  facet_wrap(~ smooth_label,scales = "free_x") +
  scale_y_continuous(limits = c(-3, 3), breaks = c(-3, 0, 3)) +
  scale_color_manual(values = model_cols, name = "Model") +
  scale_fill_manual(values = model_cols, name = "Model") +
  ylab(ch4_lab) +
  xlab(NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.title = element_text(size =15),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 14), face = "bold")

p_both <- p_big / p_small +
  plot_layout(heights = c(1, 2.75))

# ggsave(
#   filename = "Gam_CH4_partial_effects_all_models.png",
#   plot     = p_both,
#   width    = 12,
#   height   = 16,
#   dpi      = 300
# )

# ---------- Model Partial Figures (paper) ----------

# using land (big) model
mod_all <- dfmL %>%
  mutate(
    functional_group = factor(functional_group,
                              levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
    )
  ) %>%
  dplyr::select(
    ch4_mgCH4_m2_d,
    soil_moisture_avg,
    x12cm_soil_temp,
    prop_shrub,
    prop_forb,
    prop_graminoid,
    prop_Barren = prop_Barren,  
    gcc_mean,
    GPP_gC_m2_d,
    functional_group
  ) %>%
  filter(complete.cases(.))  


# cv problem fix

final_formula_cv <- as.formula(
  ch4_mgCH4_m2_d ~
    s(soil_moisture_avg, k = 4, bs = "cr") +
    s(x12cm_soil_temp,   k = 4, bs = "cr") +
    s(prop_Barren,       k = 3, bs = "cr") +
    s(prop_shrub,        k = 3, bs = "cr") +
    s(prop_graminoid,    k = 4, bs = "cr") +
    s(prop_forb,         k = 4, bs = "cr") +
    s(gcc_mean,          k = 4, bs = "cr") +
    s(GPP_gC_m2_d,       k = 4, bs = "cr") 
)


# cross-val
set.seed(123)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(mod_all)))

cv_df_land <- purrr::map_dfr(1:k, function(i){
  train <- mod_all[folds != i, ]
  test  <- mod_all[folds == i, ]
  
  m <- mgcv::gam(
    final_formula_cv,
    data = train,
    method = "REML",
    select = TRUE
  )
  
  tibble(
    obs  = test$ch4_mgCH4_m2_d,
    pred = predict(m, newdata = test),
    functional_group = test$functional_group
  )
})

# performance stats
mae  <- mean(abs(cv_df_land$obs - cv_df_land$pred), na.rm = TRUE)
rmse <- sqrt(mean((cv_df_land$obs - cv_df_land$pred)^2, na.rm = TRUE))
r2_cv <- cor(cv_df_land$obs, cv_df_land$pred, use = "complete.obs")^2
r2_in  <- summary(gam_no_prs)$r.sq
perf_label <- paste0(
  "MAE = ", round(mae, 2), "\n",
  "RMSE = ", round(rmse, 2), "\n",
  "R² (model) = ", round(r2_in, 2), "\n",
  "R² (CV) = ", round(cv_R2_no_prs, 2), "\n",
  "n = ", nrow(gam_no_prs)
)


group_means <- cv_df_land %>%
  group_by(functional_group) %>%
  summarize(
    obs_mean  = mean(obs,  na.rm = TRUE),
    pred_mean = mean(pred, na.rm = TRUE),
    .groups   = "drop"
  )

cv_df_land$functional_group <- factor(
  cv_df_land$functional_group,
  levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
)
group_means$functional_group <- factor(
  group_means$functional_group,
  levels = levels(cv_df_land$functional_group)
)

r2_in <- summary(gam_no_prs)$r.sq

perf_label <- paste0(
  "MAE = ", round(mae, 2), "\n",
  "RMSE = ", round(rmse, 2), "\n",
  "R² = ", round(r2_in, 2), "\n",
  "R² (CV) = ", round(cv_R2_no_prs, 2), "\n",
  "n = ", nrow(cv_no_prs) - 1
)

p_skill_land <- ggplot(cv_df_land, aes(x = pred, y = obs, color = functional_group)) +
  # 1:1 line
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  
  #all CV points 
  geom_point(alpha = 0.8, size = 2, stroke = 0.4) +
  
  # hollow circles =  means
  geom_point(
    data = group_means,
    aes(x = pred_mean, y = obs_mean, color = functional_group),
    size = 4.5, stroke = 1, fill = NA, shape = 21,
    inherit.aes = FALSE
  ) +
  
  # color scale
  scale_color_manual(
    values = grp_cols,
  ) +
  
  # stat box
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = perf_label,
    hjust = 1.1, vjust = -0.1,
    size = 3.2, color = "black"
  ) +
  
  # legend note for hollow circles
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = "Open circles = mean per\ndominant vegetation type",
    hjust = 1.1, vjust = 2.3,    
    size = 3.2, color = "black"
  ) +
  
  # axes
  labs(
    x = expression("Predicted CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    y = expression("Observed CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    color = ""
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.title       = element_text(face = "bold"),
    legend.position  = "right",
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10)
  )


# missing mean circle in the legend add in fig. description
p_skill_land

# ggsave("predicted_observed_CH4.png", p_skill_land, width = 8, height = 6, dpi = 400)


# mini nutrient gam model

mod_all <- dfmL %>%
  mutate(
    functional_group = factor(functional_group,
                              levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
    )
  ) %>%
  select(
    ch4_mgCH4_m2_d,
    soil_moisture_avg,
    x12cm_soil_temp,
    prop_shrub,
    prop_forb,
    prop_graminoid,
    prop_Barren = prop_Barren, 
    p,
    cu,
    total_n,
    gcc_mean,
    GPP_gC_m2_d,
    functional_group
  ) %>%
  filter(complete.cases(.))  


# cv problem fix again (does not need it)


final_formula_cv <- as.formula(
  ch4_mgCH4_m2_d ~ 
    s(soil_moisture_avg, k = 4, bs = "cr") +
    s(x12cm_soil_temp, k = 4, bs = "cr") +
    s(prop_shrub , k = 3, bs = "cr") +
    s(prop_forb, k = 3, bs = "cr") +
    s(prop_Barren, k = 3, bs = "cr") +
    s(prop_graminoid, k = 3, bs = "cr") +
    s(GPP_gC_m2_d, k = 4, bs = "cr") +
    s(total_n, k = 4, bs = "cr") +
    s(p, k = 4, bs = "cr") +
    s(cu, k = 4, bs = "cr")
)

# cross-val
set.seed(123)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(mod_all)))

cv_df <- purrr::map_dfr(1:k, function(i){
  train <- mod_all[folds != i, ]
  test  <- mod_all[folds == i, ]
  
  m <- mgcv::gam(
    final_formula_cv,
    data = train,
    method = "REML",
    select = TRUE
  )
  
  tibble(
    obs  = test$ch4_mgCH4_m2_d,
    pred = predict(m, newdata = test),
    functional_group = test$functional_group
  )
})

# performance stats
mae  <- mean(abs(cv_df$obs - cv_df$pred), na.rm = TRUE)
rmse <- sqrt(mean((cv_df$obs - cv_df$pred)^2, na.rm = TRUE))
r2_cv <- cor(cv_df$obs, cv_df$pred, use = "complete.obs")^2
r2_in  <- summary(gam_all_light_slim)$r.sq
perf_label <- paste0(
  "MAE = ", round(mae, 2), "\n",
  "RMSE = ", round(rmse, 2), "\n",
  "R² (model) = ", round(r2_in, 2), "\n",
  "R² (CV) = ", round(cv_R2_slim, 2), "\n",
  "n = ", nrow(cv_slim)
)

group_means <- cv_df %>%
  group_by(functional_group) %>%
  summarize(
    obs_mean  = mean(obs,  na.rm = TRUE),
    pred_mean = mean(pred, na.rm = TRUE),
    .groups   = "drop"
  )

cv_df$functional_group <- factor(
  cv_df$functional_group,
  levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
)
group_means$functional_group <- factor(
  group_means$functional_group,
  levels = levels(cv_df$functional_group)
)

r2_in <- summary(gam_all_light_slim)$r.sq

perf_label <- paste0(
  "MAE = ", round(mae, 2), "\n",
  "RMSE = ", round(rmse, 2), "\n",
  "R² = ", round(r2_in, 2), "\n",
  "R² (CV) = ", round(cv_R2_slim, 2), "\n",
  "n = ", nrow(cv_df)
)

p_skill <- ggplot(cv_df, aes(x = pred, y = obs, color = functional_group)) +
  # 1:1 line
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  
  #all CV points 
  geom_point(alpha = 0.8, size = 2, stroke = 0.4) +
  
  # hollow circles =  means
  geom_point(
    data = group_means,
    aes(x = pred_mean, y = obs_mean, color = functional_group),
    size = 4.5, stroke = 1, fill = NA, shape = 21,
    inherit.aes = FALSE
  ) +
  
  
  # color scale
  scale_color_manual(
    values = grp_cols,
  ) +
  
  # stat box
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = perf_label,
    hjust = 1.1, vjust = -0.1,
    size = 3.2, color = "black"
  ) +
  
  # legend note for hollow circles
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = "Open circles = mean per\ndominant vegetation type",
    hjust = 1.1, vjust = 2.3,    
    size = 3.2, color = "black"
  ) +
  
  # axes
  labs(
    x = expression("Predicted CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    y = expression("Observed CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    color = ""
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.title       = element_text(face = "bold"),
    legend.position  = "right",
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10)
  )


# missing mean circle in the legend add in fig. description
p_skill

lims <- range(c(cv_df$obs, cv_df$pred), cv_df_land$obs, cv_df_land$pred)

p_skill_mech  <- p_skill   + coord_cartesian(xlim = lims, ylim = lims)
p_skill_land  <- p_skill_land + coord_cartesian(xlim = lims, ylim = lims)

p_skill_mech
p_skill_land


clean_mech <- p_skill_mech +
  coord_cartesian(xlim = lims, ylim = lims) +
  theme(
    legend.position = "none"
  )

clean_land <- p_skill_land +
  coord_cartesian(xlim = lims, ylim = lims) +
  theme(
    legend.position = "right",
    axis.title.y = element_blank()
  )


side_by_side <- clean_mech + clean_land +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      vjust = 1
    )
  )


side_by_side

# ggsave("predicted_observed_CH4_miniGAM_cn_ratio.png", p_skill, width = 8, height = 6, dpi = 400)
# 
# ggsave("GPP_fig_obs_preds_sidebyside_cn_ratio.png",
# side_by_side, width = 12, height = 5, dpi = 600)


# ---------- Dark Models: Partial Figures (paper) ----------


mod_all_D <- dfmD %>%
  mutate(
    functional_group = factor(functional_group,
                              levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
    )
  ) %>%
  select(
    ch4_mgCH4_m2_d,
    soil_moisture_avg,
    x12cm_soil_temp,
    prop_shrub,
    prop_forb,
    prop_graminoid,
    prop_Barren = prop_Barren,  
    C_N_ratio,
    gcc_mean,
    ER_gC_m2_d,
    functional_group
  ) %>%
  filter(complete.cases(.))  

# cross-val
set.seed(123)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(mod_all_D)))

cv_df_land_D <- purrr::map_dfr(1:k, function(i){
  train <- mod_all_D[folds != i, ]
  test  <- mod_all_D[folds == i, ]
  
  m <- mgcv::gam(
    form_no_prs_D,
    data = train,
    method = "REML",
    select = TRUE
  )
  
  tibble(
    obs  = test$ch4_mgCH4_m2_d,
    pred = predict(m, newdata = test),
    functional_group = test$functional_group
  )
})

# performance stats
mae_D  <- mean(abs(cv_df_land_D$obs - cv_df_land_D$pred), na.rm = TRUE)
rmse_D <- sqrt(mean((cv_df_land_D$obs - cv_df_land_D$pred)^2, na.rm = TRUE))
r2_cv_D <- cor(cv_df_land_D$obs, cv_df_land_D$pred, use = "complete.obs")^2
r2_in_D  <- summary(gam_no_prs_D)$r.sq
perf_label_D <- paste0(
  "MAE = ", round(mae_D, 2), "\n",
  "RMSE = ", round(rmse_D, 2), "\n",
  "R² (model) = ", round(r2_in_D, 2), "\n",
  "R² (CV) = ", round(cv_R2_no_prs_D, 2), "\n",
  "n = ", nrow(cv_no_prs_D)
)


group_means_D <- cv_df_land_D %>%
  group_by(functional_group) %>%
  summarize(
    obs_mean  = mean(obs,  na.rm = TRUE),
    pred_mean = mean(pred, na.rm = TRUE),
    .groups   = "drop"
  )

cv_df_land_D$functional_group <- factor(
  cv_df_land_D$functional_group,
  levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
)
group_means_D$functional_group <- factor(
  group_means_D$functional_group,
  levels = levels(cv_df_land_D$functional_group)
)

r2_in_D <- summary(gam_no_prs_D)$r.sq

perf_label_D <- paste0(
  "MAE = ", round(mae_D, 2), "\n",
  "RMSE = ", round(rmse_D, 2), "\n",
  "R² = ", round(r2_in_D, 2), "\n",
  "R² (CV) = ", round(cv_R2_no_prs_D, 2), "\n",
  "n = ", nrow(dfmD) -1
)

p_skill_land_D <- ggplot(cv_df_land_D, aes(x = pred, y = obs, color = functional_group)) +
  # 1:1 line
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  
  #all CV points 
  geom_point(alpha = 0.8, size = 2, stroke = 0.4) +
  
  # hollow circles =  means
  geom_point(
    data = group_means_D,
    aes(x = pred_mean, y = obs_mean, color = functional_group),
    size = 4.5, stroke = 1, fill = NA, shape = 21,
    inherit.aes = FALSE
  ) +
  
  # color scale
  scale_color_manual(
    values = grp_cols,
  ) +
  
  # stat box
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = perf_label_D,
    hjust = 1.1, vjust = -0.1,
    size = 3.2, color = "black"
  ) +
  
  # legend note for hollow circles
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = "Open circles = mean per\ndominant vegetation type",
    hjust = 1.1, vjust = 2.3,    
    size = 3.2, color = "black"
  ) +
  
  # axes
  labs(
    x = expression("Predicted CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    y = expression("Observed CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    color = ""
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.title       = element_text(face = "bold"),
    legend.position  = "right",
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10)
  )

p_skill_land_D

# ggsave("predicted_observed_CH4.png", p_skill_land, width = 8, height = 6, dpi = 400)

# mini nutrient gam model

mod_all_D <- dfmD %>%
  mutate(
    functional_group = factor(functional_group,
                              levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
    )
  ) %>%
  select(
    ch4_mgCH4_m2_d,
    soil_moisture_avg,
    x12cm_soil_temp,
    prop_shrub,
    prop_forb,
    prop_graminoid,
    prop_Barren = prop_Barren, 
    p,
    cu,
    total_n,
    gcc_mean,
    co2_gC_m2_d,
    ER_gC_m2_d,
    functional_group
  ) %>%
  filter(complete.cases(.))  

# cv problem fix again (does not need it) but i like redundancy

final_formula_cv_D <- as.formula(
  ch4_mgCH4_m2_d ~ 
    s(soil_moisture_avg, k = 4, bs = "cr") +
    s(x12cm_soil_temp, k = 4, bs = "cr") +
    s(ER_gC_m2_d, k = 4, bs = "cr") +
    s(total_n, k = 4, bs = "cr") +
    s(p, k = 4, bs = "cr") +
    s(cu, k = 4, bs = "cr")
)

# cross-val
set.seed(123)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(mod_all_D)))

cv_df_D <- purrr::map_dfr(1:k, function(i){
  train <- mod_all_D[folds != i, ]
  test  <- mod_all_D[folds == i, ]
  
  m <- mgcv::gam(
    final_formula_cv_D,
    data = train,
    method = "REML",
    select = TRUE
  )
  
  tibble(
    obs  = test$ch4_mgCH4_m2_d,
    pred = predict(m, newdata = test),
    functional_group = test$functional_group
  )
})

# performance stats
mae_D  <- mean(abs(cv_df_D$obs - cv_df_D$pred), na.rm = TRUE)
rmse_D <- sqrt(mean((cv_df_D$obs - cv_df_D$pred)^2, na.rm = TRUE))
r2_cv_D <- cor(cv_df_D$obs, cv_df_D$pred, use = "complete.obs")^2
r2_in_D  <- summary(gam_all_light_slim_D)$r.sq
perf_label_D <- paste0(
  "MAE = ", round(mae_D, 2), "\n",
  "RMSE = ", round(rmse_D, 2), "\n",
  "R² (model) = ", round(r2_in_D, 2), "\n",
  "R² (CV) = ", round(cv_R2_slim_D, 2), "\n",
  "n = ", nrow(dfmD)
)


group_means_D <- cv_df_D %>%
  group_by(functional_group) %>%
  summarize(
    obs_mean  = mean(obs,  na.rm = TRUE),
    pred_mean = mean(pred, na.rm = TRUE),
    .groups   = "drop"
  )

cv_df_D$functional_group <- factor(
  cv_df_D$functional_group,
  levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
)
group_means_D$functional_group <- factor(
  group_means_D$functional_group,
  levels = levels(cv_df_D$functional_group)
)



r2_in_D <- summary(gam_all_light_slim_D)$r.sq

perf_label_D <- paste0(
  "MAE = ", round(mae_D, 2), "\n",
  "RMSE = ", round(rmse_D, 2), "\n",
  "R² = ", round(r2_in_D, 2), "\n",
  "R² (CV) = ", round(cv_R2_slim_D, 2), "\n",
  "n = ", nrow(cv_df_D)
)

p_skill_D <- ggplot(cv_df_D, aes(x = pred, y = obs, color = functional_group)) +
  # 1:1 line
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  
  #all CV points 
  geom_point(alpha = 0.8, size = 2, stroke = 0.4) +
  
  # hollow circles =  means
  geom_point(
    data = group_means_D,
    aes(x = pred_mean, y = obs_mean, color = functional_group),
    size = 4.5, stroke = 1, fill = NA, shape = 21,
    inherit.aes = FALSE
  ) +
  
  
  # color scale
  scale_color_manual(
    values = grp_cols,
  ) +
  
  # stat box
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = perf_label_D,
    hjust = 1.1, vjust = -0.1,
    size = 3.2, color = "black"
  ) +
  
  # legend note for hollow circles
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = "Open circles = mean per\ndominant vegetation type",
    hjust = 1.1, vjust = 2.3,    
    size = 3.2, color = "black"
  ) +
  
  # axes
  labs(
    x = expression("Predicted CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    y = expression("Observed CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    color = ""
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.title       = element_text(face = "bold"),
    legend.position  = "right",
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10)
  )


# missing mean circle in the legend add in fig. description
p_skill_D

lims <- range(c(cv_df_D$obs, cv_df_D$pred), cv_df_land_D$obs, cv_df_land_D$pred)

p_skill_mech_D  <- p_skill_D   + coord_cartesian(xlim = lims, ylim = lims)
p_skill_land_D  <- p_skill_land_D + coord_cartesian(xlim = lims, ylim = lims)

p_skill_mech_D
p_skill_land_D


clean_mech_D <- p_skill_mech_D +
  coord_cartesian(xlim = lims, ylim = lims) +
  theme(
    legend.position = "none"
  )

clean_land_D <- p_skill_land_D +
  coord_cartesian(xlim = lims, ylim = lims) +
  theme(
    legend.position = "right",
    axis.title.y = element_blank()
  )

side_by_side_D <- clean_mech_D + clean_land_D +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      vjust = 1
    )
  )

side_by_side_D

# ggsave("predicted_observed_CH4_miniGAM.png", p_skill, width = 8, height = 6, dpi = 400)
# 
# ggsave("Dark_fig_obs_preds_sidebyside.png", side_by_side_D,
#        width = 12, height = 5, dpi = 600)


# --------------- LIght Variable Importance (paper) ------------------

# prep data
mod_all2 <- dfmL %>%
  mutate(
    # collapse shrubs (not enough obs in erect)
    veg_group = case_when(
      functional_group %in% c("Erect Shrub", "Prostrate Shrub") ~ "Shrub",
      functional_group %in% c("Forb", "Graminoid", "Barren")    ~ functional_group,
      TRUE                                                      ~ NA_character_
    ),
    veg_group = factor(veg_group, levels = c("Shrub","Forb","Graminoid","Barren"))
  ) %>%
  select(
    ch4_mgCH4_m2_d,
    soil_moisture_avg,
    x12cm_soil_temp,
    avg_par,
    prop_shrub,
    prop_forb,
    prop_graminoid,
    prop_Barren = prop_Barren,
    GPP_gC_m2_d, 
    total_n,
    p,
    cu,
    veg_group
  ) %>%
  filter(complete.cases(.))

summary(gam_all_light_slim)$r.sq   
summary(gam_all_light_slim)$dev.expl * 100

smooth_terms <- c(
  's(soil_moisture_avg, k = 4, bs = "cr")',
  's(x12cm_soil_temp,   k = 4, bs = "cr")',
  's(GPP_gC_m2_d,       k = 4, bs = "cr")',
  's(avg_par,       k = 4, bs = "cr")',
  's(total_n,           k = 4, bs = "cr")',
  's(p,                 k = 4, bs = "cr")',
  's(cu,                k = 4, bs = "cr")'
)


set.seed(42)

# drop-one-term importance 100 times on bootstrap resamples
# need to redo with more vals <- did minimum for checks

boot_imp <- map_dfr(1:100, function(i){
  samp <- mod_all2[sample(nrow(mod_all2), replace = TRUE), ]
  out  <- get_importance_safegam(
    data_in = samp,
    response = "ch4_mgCH4_m2_d",
    smooth_terms = smooth_terms
  )$imp %>%
    mutate(iter = i)
})

#  mean ± sd
imp_summary <- boot_imp %>%
  mutate(variable = clean_term(term)) %>%
  group_by(variable) %>%
  summarise(
    mean_delta = mean(delta_r2, na.rm = TRUE),
    sd_delta   = sd(delta_r2, na.rm = TRUE),
    se_delta   = sd_delta / sqrt(n())   # standard error
  ) %>%
  ungroup()

fig_imp_sd <- ggplot(imp_summary,
                     aes(y = fct_reorder(variable, mean_delta), x = mean_delta)
) +
  # SE whiskers
  geom_errorbarh(
    aes(xmin = mean_delta - se_delta, xmax = mean_delta + se_delta),
    height = 0, color = "grey50", linewidth = 0.8
  ) +
  geom_point(size = 3, color = "black", fill = "grey25", shape = 21) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  labs(
    x = expression(Delta*R^2*" (mean ± SE, bootstrap)"),
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y        = element_text(color = "black"),
    axis.text.x        = element_text(color = "black"),
    axis.title.x       = element_text(face = "bold"),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

fig_imp_sd

dev_exp <- summary(gam_all_light_slim)$dev.expl * 100
r2      <- summary(gam_all_light_slim)$r.sq

# needs colors (maybe - not convinced on the colors yet)
imp_summary <- imp_summary %>%
  mutate(type = case_when(
    variable %in% c("GPP", "Soil Temperature (12 cm)", "Mean Soil Moisture") ~ "Biophysical",
    variable %in% c("P", "Total N", "Cu") ~ "Available Nutrients",
    TRUE ~ "Vegetation"
  ))


# cool plot / need to make pretty 

imp_summary <- imp_summary %>%
  mutate(
    grp = case_when(
      variable %in% c("Erect Shrub", "Erect Shrub Cover")       ~ "Erect Shrub",
      variable %in% c("Prostrate Shrub", "Shrub Cover")         ~ "Prostrate Shrub",  # ALL shrubs here
      variable %in% c("Forb", "Forb Cover")                     ~ "Forb",
      variable %in% c("Graminoid", "Graminoid Cover")           ~ "Graminoid",
      variable %in% c("Barren", "Barren Cover")                 ~ "Barren",
      variable %in% c("NEE", "Soil Temperature (12 cm)", "Mean Soil Moisture") ~ "Biophysical",
      variable %in% c("P", "Total N", "Cu") ~ "Available Nutrients",
      
      TRUE ~ "Other"
    )
  )

full_cols <- c(
  grp_cols,
  "Biophysical"        = "grey50",
  "Available Nutrients" = "grey75",
  "Other"               = "grey50"
)

fig_imp <- ggplot(imp_summary, aes(
  y = fct_reorder(variable, mean_delta),
  x = mean_delta,
  color = grp
)) +
  
  # SE whiskers
  geom_errorbarh(
    aes(xmin = mean_delta - se_delta,
        xmax = mean_delta + se_delta),
    height = 0,
    linewidth = 0.8
  ) +
  
  # Vertical zero line
  geom_vline(
    xintercept = 0, linetype = "dashed",
    color = "grey65", linewidth = 0.6
  ) +
  
  # Mean points
  geom_point(size = 3.8, shape = 19, stroke = 1) +
  
  # Labels
  labs(
    x = "Variable Importance",
    y = NULL,
    caption = paste0(
      "Explained variance = ", round(dev_exp, 0),
      "%   |   R² = ", round(r2, 2)
      
      # "   |   n = ", nrow(cv_df)
    )
  ) +
  
  # Colors
  scale_color_manual(values = full_cols) +
  
  # Theme polish
  theme_minimal(base_size = 13.8) +
  theme(
    # Clean everything
    panel.grid = element_blank(),
    
    # Axes
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = ggplot2::margin(t = 5)),
    
    # Subtle border for publication
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.6
    ),
    
    # Remove legend entirely
    legend.position = "none",
    
    # Center the caption
    plot.caption = element_text(
      hjust = 0.5,
      size = 11,
      margin = ggplot2::margin(t = 10),
      color = "grey25"
    ),
    
    # Spacing
    plot.margin = ggplot2::margin(10, 18, 10, 10)
  )



fig_imp

# ggsave("fig_imp_miniGAM.png", fig_imp, width = 10, height = 6, dpi = 320)
# 
# 
# big gam 


get_delta_r2_once <- function(dat, response, vars) {
  full_form <- as.formula(
    paste0(
      response, " ~ ",
      paste(sprintf("s(%s, bs='cr')", vars), collapse = " + ")
    )
  )
  
  full_fit <- gam(full_form, data = dat, method = "REML")
  base_r2  <- summary(full_fit)$r.sq
  
  out <- map_dfr(vars, function(v) {
    reduced_vars <- setdiff(vars, v)
    red_form <- as.formula(
      paste0(
        response, " ~ ",
        paste(sprintf("s(%s, bs='cr')", reduced_vars), collapse = " + ")
      )
    )
    red_fit <- gam(red_form, data = dat, method = "REML")
    r2_red  <- summary(red_fit)$r.sq
    
    tibble(
      variable  = v,
      delta_r2  = base_r2 - r2_red
    )
  })
  
  out
}

vars_land <- c(
  "soil_moisture_avg",
  "x12cm_soil_temp",
  "prop_shrub",
  "prop_forb",
  "prop_graminoid",
  "prop_Barren",
  "gcc_mean",
  "avg_par",
  "C_N_ratio",
  "GPP_gC_m2_d"
)

set.seed(123)

B <- 100  # need to do 50 or 100 - did minimum

boot_list <- map(1:B, function(b) {
  idx <- sample(seq_len(nrow(mod_no_prs)), replace = TRUE)
  dat_b <- mod_no_prs[idx, ]
  get_delta_r2_once(dat_b, response = "ch4_mgCH4_m2_d", vars = vars_land)
})

boot_df <- bind_rows(boot_list, .id = "boot_id")


imp_land_summary <- boot_df %>%
  group_by(variable) %>%
  summarise(
    mean_delta = mean(delta_r2, na.rm = TRUE),
    sd_delta   = sd(delta_r2,   na.rm = TRUE),
    se_delta   = sd_delta / sqrt(n()),   # standard error
    .groups = "drop"
  )


pretty_names <- c(
  soil_moisture_avg = "Soil Moisture",
  x12cm_soil_temp   = "Soil Temperature (12 cm)",
  prop_shrub        = "Shrub Cover",
  prop_forb         = "Forb Cover",
  prop_graminoid    = "Graminoid Cover",
  prop_Barren       = "Barren Cover",
  C_N_ratio         = "C/N",
  avg_par  = "PAR",
  gcc_mean          = "Canopy Greenness (GCC)",
  GPP_gC_m2_d        = "GPP"
)

imp_land_summary <- imp_land_summary %>%
  mutate(
    variable_label = pretty_names[variable],
    
    type = case_when(
      variable %in% c("GPP_gC_m2_d",
                      "x12cm_soil_temp",
                      "soil_moisture_avg") ~ "Biophysical",
      
      variable %in% c("gcc_mean",
                      "prop_shrub",
                      "prop_forb",
                      "prop_graminoid",
                      "prop_Barren") ~ "Vegetation",
      
      TRUE ~ "Other"
    )
  )


imp_land_summary <- imp_land_summary %>%
  mutate(
    se_delta = sd_delta / sqrt(n()),
    
    grp = case_when(
      variable == "prop_shrub"      ~ "Prostrate Shrub",
      variable == "prop_forb"       ~ "Forb",
      variable == "prop_graminoid"  ~ "Graminoid",
      variable == "prop_Barren"     ~ "Barren",
      
      variable %in% c("GPP_gC_m2_d", "x12cm_soil_temp", "soil_moisture_avg") ~ "Biophysical",
      
      TRUE ~ "Other"
    )
  )


full_cols_land <- c(
  grp_cols,
  "Biophysical" = "grey50",
  "Other"       = "grey50"
)


dev_exp_land <- summary(gam_no_prs)$dev.expl * 100
r2_land      <- summary(gam_no_prs)$r.sq

fig_imp_land <- ggplot(
  imp_land_summary,
  aes(
    y = fct_reorder(variable_label, mean_delta),
    x = mean_delta,
    color = grp
  )
) +
  
  # SE whiskers
  geom_errorbarh(
    aes(
      xmin = mean_delta - se_delta,
      xmax = mean_delta + se_delta
    ),
    height = 0,
    linewidth = 0.8
  ) +
  
  # Mean points
  geom_point(
    size = 3.8,
    shape = 19,
    stroke = 1
  ) +
  
  # Zero vertical line
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.6
  ) +
  
  # Colors
  scale_color_manual(values = full_cols_land) +
  scale_y_discrete(position = "right", 
                   expand = expansion(add = 0.4)
                   ) +

  
  labs(
    x = "Variable Importance",
    y = NULL,
    caption = paste0(
      "Explained variance = ", round(dev_exp_land, 0), 
      "%   |   R² = ", round(r2_land, 2),
      "   |   n = ", nrow(mod_no_prs)
    )
  ) +
  
  theme_minimal(base_size = 13.8) +
  theme(
    panel.grid = element_blank(),
    
    # axis + ticks
    axis.line  = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text  = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = ggplot2::margin(t = 5)),
    axis.text.y.right = element_text(
      hjust = 0,
      size = 12
    ),
    
    # border
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.6
    ),
    
    # remove legend
    legend.position = "none",
    
    # center caption
    plot.caption = element_text(
      hjust = 0.5,
      size  = 11,
      margin = ggplot2::margin(t = 10),
      color = "grey25"
    ),
    
    plot.margin = ggplot2::margin(10, 18, 10, 10)
  )


fig_imp_land


# 
# ggsave("NEE_fig_imp_landGAM.png", fig_imp_land, width = 10, height = 6, dpi = 320)

ggplot(mod_all, aes(GPP_gC_m2_d, avg_par)) + geom_point() + geom_smooth()

concurvity(gam_all_light_slim, full = TRUE)


fig_cd <- fig_imp + fig_imp_land +
  plot_layout(ncol = 2) &
  theme(plot.tag = element_blank())

fig_cd <- fig_cd +
  plot_annotation(
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  )

fig_cd[[1]] <- fig_cd[[1]] + labs(tag = "c")
fig_cd[[2]] <- fig_cd[[2]] + labs(tag = "d")


fig_cd

# ---------------- Dark Variable Importance (paper) -----------

# prep data

mod_all2 <- dfmD %>%
  mutate(
    # collapse shrubs (not enough obs in erect)
    veg_group = case_when(
      functional_group %in% c("Erect Shrub", "Prostrate Shrub") ~ "Shrub",
      functional_group %in% c("Forb", "Graminoid", "Barren")    ~ functional_group,
      TRUE                                                      ~ NA_character_
    ),
    veg_group = factor(veg_group, levels = c("Shrub","Forb","Graminoid","Barren"))
  ) %>%
  select(
    ch4_mgCH4_m2_d,
    soil_moisture_avg,
    x12cm_soil_temp,
    prop_shrub,
    prop_forb,
    avg_par,
    prop_graminoid,
    prop_Barren = prop_Barren,
    ER_gC_m2_d,
    total_n,
    p,
    cu,
    mn,
    s,
    veg_height_mean,
    C_N_ratio,
    slope_aspect,
    veg_group
  ) %>%
  filter(complete.cases(.))



summary(gam_all_light_slim_D)$r.sq   
summary(gam_all_light_slim_D)$dev.expl * 100

smooth_terms <- c(
  's(soil_moisture_avg, k = 4, bs = "cr")',
  's(x12cm_soil_temp,   k = 4, bs = "cr")',
  's(ER_gC_m2_d,       k = 4, bs = "cr")',
  's(total_n,           k = 4, bs = "cr")',
  's(p,                 k = 4, bs = "cr")',
  's(cu,                k = 4, bs = "cr")'
)

set.seed(42)

# drop-one-term importance 100 times on bootstrap resamples
# need to redo with more vals <- did minimum for checks

boot_imp <- map_dfr(1:100, function(i){
  samp <- mod_all2[sample(nrow(mod_all2), replace = TRUE), ]
  out  <- get_importance_safegam(
    data_in = samp,
    response = "ch4_mgCH4_m2_d",
    smooth_terms = smooth_terms
  )$imp %>%
    mutate(iter = i)
})

#  mean ± sd
imp_summary <- boot_imp %>%
  mutate(variable = clean_term(term)) %>%
  group_by(variable) %>%
  summarise(
    mean_delta = mean(delta_r2, na.rm = TRUE),
    sd_delta   = sd(delta_r2, na.rm = TRUE)
  ) %>%
  ungroup()

fig_imp_sd_D <- ggplot(imp_summary,
                     aes(y = fct_reorder(variable, mean_delta), x = mean_delta)
) +
  # SD whiskers
  geom_errorbarh(
    aes(xmin = mean_delta - sd_delta, xmax = mean_delta + sd_delta),
    height = 0, color = "grey50", linewidth = 0.8
  ) +
  # central circle
  geom_point(size = 3, color = "black", fill = "grey25", shape = 21) +
  
  # zero line to show direction
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  
  labs(
    x = expression(Delta*R^2*" (mean ± SD, bootstrap)"),
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y        = element_text(color = "black"),
    axis.text.x        = element_text(color = "black"),
    axis.title.x       = element_text(face = "bold"),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

fig_imp_sd_D

dev_exp <- summary(gam_all_light_slim_D)$dev.expl * 100
r2      <- summary(gam_all_light_slim_D)$r.sq

# needs colors (maybe - not convinced on the colors yet)
imp_summary <- imp_summary %>%
  mutate(
    se_delta = sd_delta / sqrt(n()),
    
    grp = case_when(
      # Vegetation covers → use your custom grp_cols
      variable == "Shrub Cover"     ~ "Prostrate Shrub",
      variable == "Forb Cover"      ~ "Forb",
      variable == "Graminoid Cover" ~ "Graminoid",
      variable == "Barren Cover"    ~ "Barren",
      
      # Biophysical
      variable %in% c("ER", "Soil Temperature (12 cm)", "Soil Moisture") ~ "Biophysical",
      
      # Nutrients
      variable %in% c("P", "Total N", "Cu", "Mn", "S") ~ "Nutrients",
      
      TRUE ~ "Other"
    )
  )

full_cols <- c(
  grp_cols,
  "Biophysical"        = "grey50",
  "Nutrients" = "grey75",
  "Other"               = "grey50"
)

fig_imp_D <- ggplot(
  imp_summary,
  aes(
    y = fct_reorder(variable, mean_delta),
    x = mean_delta,
    color = grp
  )
) +
  
  # Zero reference (drawn first = behind)
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.6
  ) +
  
  # SE whiskers
  geom_errorbarh(
    aes(
      xmin = mean_delta - se_delta,
      xmax = mean_delta + se_delta
    ),
    height = 0,
    linewidth = 0.8
  ) +
  
  # Points
  geom_point(size = 3.8, shape = 19, stroke = 1) +
  
  scale_color_manual(values = full_cols) +
  
  labs(
    x = "Variable Importance",
    y = NULL,
    caption = paste0(
      "Explained variance = ", round(dev_exp, 0),
      "%   |   R² = ", round(r2, 2),
      "   |   n = ", nrow(mod_all_D)
    )
  ) +
  
  theme_minimal(base_size = 13.8) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text  = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = ggplot2::margin(t = 5)),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "none",
    plot.caption = element_text(
      hjust = 0.5,
      size = 11,
      margin = ggplot2::margin(t = 10),
      color = "grey25"
    ),
    plot.margin = ggplot2::margin(10, 18, 10, 10)
  )


fig_imp_D 

# ggsave("Dark_fig_imp_miniGAM.png", fig_imp_D, width = 10, height = 6, dpi = 320)

# big gam 


get_delta_r2_once <- function(dat, response, vars) {
  full_form <- as.formula(
    paste0(
      response, " ~ ",
      paste(sprintf("s(%s, bs='cr')", vars), collapse = " + ")
    )
  )
  
  full_fit <- gam(full_form, data = dat, method = "REML")
  base_r2  <- summary(full_fit)$r.sq
  
  out <- map_dfr(vars, function(v) {
    reduced_vars <- setdiff(vars, v)
    red_form <- as.formula(
      paste0(
        response, " ~ ",
        paste(sprintf("s(%s, bs='cr')", reduced_vars), collapse = " + ")
      )
    )
    red_fit <- gam(red_form, data = dat, method = "REML")
    r2_red  <- summary(red_fit)$r.sq
    
    tibble(
      variable  = v,
      delta_r2  = base_r2 - r2_red
    )
  })
  
  out
}

vars_land <- c(
  "soil_moisture_avg",
  "x12cm_soil_temp",
  "prop_shrub",
  "prop_forb",
  "prop_graminoid",
  "prop_Barren",
  "gcc_mean",
  "C_N_ratio",
  "ER_gC_m2_d"
)

set.seed(123)

B <- 100  # need to do 50 or 100 - did minimum

boot_list <- map(1:B, function(b) {
  idx <- sample(seq_len(nrow(mod_no_prs_D)), replace = TRUE)
  dat_b <- mod_no_prs_D[idx, ]
  get_delta_r2_once(dat_b, response = "ch4_mgCH4_m2_d", vars = vars_land)
})

boot_df <- bind_rows(boot_list, .id = "boot_id")


imp_land_summary_D <- boot_df %>%
  group_by(variable) %>%
  summarise(
    mean_delta = mean(delta_r2, na.rm = TRUE),
    sd_delta   = sd(delta_r2,   na.rm = TRUE),
    se_delta   = sd_delta / sqrt(n()),   # standard error
    .groups = "drop"
  )


pretty_names <- c(
  soil_moisture_avg = "Soil Moisture",
  x12cm_soil_temp   = "Soil Temperature (12 cm)",
  prop_shrub        = "Shrub Cover",
  prop_forb         = "Forb Cover",
  prop_graminoid    = "Graminoid Cover",
  prop_Barren       = "Barren Cover",
  gcc_mean          = "Canopy Greenness (GCC)",
  C_N_ratio          = "C/N",
  ER_gC_m2_d        = "ER"
)

imp_land_summary_D <- imp_land_summary_D %>%
  mutate(
    variable_label = pretty_names[variable],
    
    type = case_when(
      variable %in% c("co2_gC_m2_d",
                      "x12cm_soil_temp",
                      "soil_moisture_avg") ~ "Biophysical",
      
      variable %in% c("gcc_mean",
                      "prop_shrub",
                      "prop_forb",
                      "prop_graminoid",
                      "prop_Barren") ~ "Vegetation",
      
      TRUE ~ "Other"
    )
  )


imp_land_summary_D <- imp_land_summary_D %>%
  mutate(
    se_delta = sd_delta / sqrt(n()),
    
    grp = case_when(
      variable == "prop_shrub"      ~ "Prostrate Shrub",
      variable == "prop_forb"       ~ "Forb",
      variable == "prop_graminoid"  ~ "Graminoid",
      variable == "prop_Barren"     ~ "Barren",
      
      variable %in% c("ER_gC_m2_d", "x12cm_soil_temp", "soil_moisture_avg"),
      
      TRUE ~ "Other"
    )
  )


full_cols_land <- c(
  grp_cols,
  "Biophysical" = "grey50",
  "Other"       = "grey50"
)


dev_exp_land <- summary(gam_no_prs_D)$dev.expl * 100
r2_land      <- summary(gam_no_prs_D)$r.sq

fig_imp_land_D <- ggplot(
  imp_land_summary_D,
  aes(
    y = fct_reorder(variable_label, mean_delta),
    x = mean_delta,
    color = grp
  )
) +
  
  # SE whiskers
  geom_errorbarh(
    aes(
      xmin = mean_delta - se_delta,
      xmax = mean_delta + se_delta
    ),
    height = 0,
    linewidth = 0.8
  ) +
  
  # Mean points
  geom_point(
    size = 3.8,
    shape = 19,
    stroke = 1
  ) +
  
  # Zero vertical line
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.6
  ) +
  
  # Colors
  scale_color_manual(values = full_cols_land) +
  scale_y_discrete(position = "right", 
                   expand = expansion(add = 0.4)
  ) +
  
  labs(
    x = "Variable Importance",
    y = NULL,
    caption = paste0(
      "Explained variance = ", round(dev_exp_land, 0), 
      "%   |   R² = ", round(r2_land, 2),
      "   |   n = ", nrow(mod_no_prs)
    )
  ) +
  
  theme_minimal(base_size = 13.8) +
  theme(
    panel.grid = element_blank(),
    
    # axis + ticks
    axis.line  = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text  = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = ggplot2::margin(t = 5)),
    axis.text.y.right = element_text(
      hjust = 0,
      size = 12
    ),
    
    # border
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.6
    ),
    
    # remove legend
    legend.position = "none",
    
    # center caption
    plot.caption = element_text(
      hjust = 0.5,
      size  = 11,
      margin = ggplot2::margin(t = 10),
      color = "grey25"
    ),
    
    plot.margin = ggplot2::margin(10, 18, 10, 10)
  )


fig_imp_land_D

# ggsave("Dark_fig_imp_landGAM_CN.png", fig_imp_land_D, width = 10, height = 6, dpi = 320)



fig_cd_D <- fig_imp_D + fig_imp_land_D +
  plot_layout(ncol = 2) &
  theme(plot.tag = element_blank())

fig_cd_D <- fig_cd_D +
  plot_annotation(
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  )

fig_cd_D[[1]] <- fig_cd_D[[1]] + labs(tag = "c")
fig_cd_D[[2]] <- fig_cd_D[[2]] + labs(tag = "d")


fig_cd_D

# ggsave(
#   filename = "DARK_Fig_c_d_importance.png",
#   plot     = fig_cd_D,
#   width    = 16,
#   height   = 6,
#   dpi      = 600,
#   units    = "in",
#   bg       = "white"
# )


# ------------------ PRS Figures --------------------
# order 
levs <- c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")

theme_pub <- function(base_size = 10){
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(linewidth = 0.3, color = "grey85"),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = base_size + 2),
      legend.position  = "none",
      strip.text       = element_text(face = "bold")
    )
}

avail <- df_master %>%
  select(functional_group,
         total_n, p, cu, ch4_mgCH4_m2_d # need methane for later 
  ) %>%
  filter(!is.na(functional_group))

# order groups 
avail <- avail %>% mutate(functional_group = factor(functional_group, levels = levs))

make_avail_plot <- function(df, value_col, title_text, element_symbol){
  stats_df <- df %>%
    group_by(functional_group) %>%
    summarise(
      mean_val = mean(.data[[value_col]], na.rm = TRUE),
      sd_val   = sd(.data[[value_col]],  na.rm = TRUE),
      .groups  = "drop"
    )
  
  ggplot(df, aes(x = functional_group, y = .data[[value_col]], fill = functional_group)) +
    # background boxplot
    geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.4, alpha = 0.75) +
    geom_jitter(aes(color = functional_group),
                width = 0.15, height = 0,
                alpha = 0.6, size = 2, stroke = 0, show.legend = FALSE) +
    # mean ± SD 
    geom_errorbar(data = stats_df,
                  aes(x = functional_group, ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.5, color = "black") +
    geom_point(data = stats_df,
               aes(x = functional_group, y = mean_val),
               inherit.aes = FALSE, shape = 21, size = 2.8, stroke = 0.35, fill = "black") +
    # color scales
    scale_fill_manual(values = grp_cols, drop = FALSE) +
    scale_color_manual(values = grp_cols, drop = FALSE) +
    labs(
      title = title_text,
      x = NULL,
      y = bquote(.(title_text)~"availability [ "*mg~.(as.name(element_symbol))*m^{-2}*d^{-1}*"]")
    ) +
    theme_pub() +
    theme(
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      panel.grid.major.x = element_blank(),   # remove vertical grid lines
      panel.grid.minor.x = element_blank(),   # remove lil vertical lines
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),  # keep horizontal lines
      panel.grid.minor.y = element_blank()
    )
}

# panels with borders
pN  <- make_avail_plot(
  avail, "total_n",
  bquote("Total Nitrogen ("*NH[4]^"+"*" + "*NO[3]^"-"*")"),
  "N"
)

pP  <- make_avail_plot(
  avail, "p",
  bquote("Phosphorus ("*H[2]*PO[4]^"-"*")"),
  "P"
)


pCu <- make_avail_plot(
  avail, "cu",
  bquote("Copper ("*Cu^"2+"*")"),
  "Cu"
)

fig_avail <- pN | pCu | pP
fig_avail


# ggsave("fig_PRS_availability_boxplots.png", fig_avail, width = 14, height = 6, dpi = 320)


# fluxes vs prs

nutrient_vars <- c("cu", "total_n", "p")

df_nutr <- mod_all %>%
  select(ch4_mgCH4_m2_d, all_of(nutrient_vars)) %>%
  tidyr::pivot_longer(cols = all_of(nutrient_vars),
                      names_to = "nutrient",
                      values_to = "value") %>%
  mutate(
    nutrient = factor(
      nutrient,
      levels = c("total_n", "p", "cu"),
      labels = c("Total N (mg/kg)", "Soil P (mg/kg)", "Soil Cu (mg/kg)")
    )
  )

# Base plot
p_nutrients <- ggplot(df_nutr, aes(x = value, y = ch4_mgCH4_m2_d)) +
  geom_point(size = 2, alpha = 0.6, color = "grey30") +
  geom_smooth(
    method = "gam", formula = y ~ s(x, bs = "cr"),
    color = "#0072B2", fill = "#0072B2", alpha = 0.25, linewidth = 1
  ) +
  facet_wrap(~ nutrient, scales = "free_x", nrow = 1) +
  labs(
    x = NULL,
    y = expression("CH"[4] * " flux (mg C m"^{-2}*" d"^{-1}*")"),
    title = "Nutrient controls on methane flux"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

p_nutrients


# long format for faceting
avail_long <- avail %>%
  select(ch4_mgCH4_m2_d, functional_group, total_n, p, cu) %>%
  pivot_longer(cols = c(total_n, p, cu),
               names_to = "nutrient", values_to = "value") %>%
  mutate(
    nutrient = factor(
      nutrient,
      levels = c("total_n", "p", "cu"),
      labels = c(
        expression("Total Nitrogen ("*NH[4]^"+"*" + "*NO[3]^"-"*")"),
        expression("Phosphorus ("*H[2]*PO[4]^"-"*")"),
        expression("Copper ("*Cu^"+"*")")
      )
    )
  )


# plot
p_nutrients_veg <- ggplot(avail_long, aes(x = value, y = ch4_mgCH4_m2_d, color = functional_group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(
    method = "gam", formula = y ~ s(x, bs = "cr"),
    se = TRUE, color = "black", fill = "grey80", alpha = 0.3, linewidth = 0.8
  ) +
  facet_wrap(~ nutrient, scales = "free_x", nrow = 1,
             labeller = label_parsed) + 
  scale_color_manual(values = grp_cols) +
  labs(
    y = expression(CH[4]*" flux (mg C m"^{-2}*" d"^{-1}*")"),
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

p_nutrients_veg


p_nutrients_veg_global <- ggplot(
  avail_long,
  aes(x = value, y = ch4_mgCH4_m2_d, color = functional_group)
) +
  geom_point(size = 1.7, alpha = 0.35) +  #transparent-ish dots
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cr"),
    se = TRUE,
    color = "black",          
    fill  = "grey80",
    alpha = 0.4,
    linewidth = 1.0
  ) +
  facet_wrap(
    ~ nutrient,
    scales   = "free_x",
    nrow     = 1,
    labeller = label_parsed  
  ) +
  scale_color_manual(values = grp_cols, name = " ") +
  labs(
    y = expression(CH[4]*" flux (mg C m"^{-2}*" d"^{-1}*")"),
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    strip.text       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold")
  )

p_nutrients_veg_global

p_nutrients_veg_global<- p_nutrients_veg_global +
  coord_cartesian(ylim = c(-5, 2)) 
# ggsave("flux_nuts.png", p_nutrients_veg_global, width = 10, height = 4.9, units = "in", dpi = 300)


# --------------- variable Importance (FROM Lai et al., 2024) -----------------


run_gam_hp <- function(fit) {
  if (exists("gam_hp", where = asNamespace("gam.hp"), inherits = FALSE)) {
    gam.hp::gam_hp(fit)
  } else if (exists("gam.hp", where = asNamespace("gam.hp"), inherits = FALSE)) {
    gam.hp::gam.hp(fit)
  } else {
    stop("gam.hp function not found. Run ls('package:gam.hp') to see exported names.")
  }
}

hp_land <- run_gam_hp(gam_no_prs)
hp_mini <- run_gam_hp(gam_all_light_slim)


land_tbl <- as.data.frame(hp_land$individual) %||% as.data.frame(hp_land$R2_individual) %||% hp_land
mini_tbl <- as.data.frame(hp_mini$individual) %||% as.data.frame(hp_mini$R2_individual) %||% hp_mini

clean_term <- function(x_vec) {
  base <- x_vec %>%
    str_replace("^s\\(", "") %>%
    str_replace(",.*$", "")
  lookup <- c(
    soil_moisture_avg = "Soil Moisture",
    x12cm_soil_temp   = "Soil Temperature (12 cm)",
    prop_shrub        = "Shrub Cover",
    prop_forb         = "Forb Cover",
    prop_Barren       = "Barren Cover",
    prop_graminoid    = "Graminoid Cover",
    GPP_gC_m2_d       = "GPP",
    total_n           = "Total N",
    veg_height_mean   = "Mean Veg. Height",
    p                 = "P",
    cu                = "Cu",
    gcc_mean          = "Canopy Greenness (GCC)"
  )
  ifelse(base %in% names(lookup), lookup[base], base)
}

guess_cols <- function(df) {
  nm <- names(df)
  term_col <- nm[which.max(sapply(nm, function(n) mean(grepl("s\\(", df[[n]]))))] %||% nm[1]
  r2_col   <- nm[which.max(sapply(nm, function(n) is.numeric(df[[n]]) && max(df[[n]], na.rm=TRUE) <= 1.0001))]
  list(term = term_col, r2 = r2_col)
}

pick_land <- guess_cols(land_tbl)
pick_mini <- guess_cols(mini_tbl)

hp_land_clean <- land_tbl %>%
  transmute(variable = clean_term(.data[[pick_land$term]]),
            hp_R2     = as.numeric(.data[[pick_land$r2]])) %>%
  group_by(variable) %>% summarise(hp_R2 = sum(hp_R2, na.rm=TRUE), .groups="drop")

hp_mini_clean <- mini_tbl %>%
  transmute(variable = clean_term(.data[[pick_mini$term]]),
            hp_R2     = as.numeric(.data[[pick_mini$r2]])) %>%
  group_by(variable) %>% summarise(hp_R2 = sum(hp_R2, na.rm=TRUE), .groups="drop")

imp_mini_ready <- imp_summary %>%
  select(variable, mean_delta, sd_delta) %>%
  mutate(model = "Mini (Nutrients)")

imp_land_ready <- imp_land_summary %>%
  transmute(variable = pretty_names[variable], mean_delta, sd_delta) %>%
  mutate(model = "Landscape")

mini_join <- imp_mini_ready %>%
  left_join(hp_mini_clean, by = "variable") %>%
  mutate(hp_R2 = 100 * hp_R2)  # convert to percent if HP returns proportions

land_join <- imp_land_ready %>%
  left_join(hp_land_clean, by = "variable") %>%
  mutate(hp_R2 = 100 * hp_R2)

plot_imp <- function(df, title_txt, dev_exp, r2_model, n_label = NULL) {
  ggplot(df, aes(y = fct_reorder(variable, mean_delta))) +
    geom_errorbarh(aes(xmin = mean_delta - sd_delta,
                       xmax = mean_delta + sd_delta),
                   height = 0, color = "grey55") +
    geom_point(aes(x = mean_delta), size = 3.4, shape = 16) +
    geom_point(aes(x = hp_R2/100), size = 3.4, shape = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    labs(
      title   = title_txt,
      x       = "Importance (unique ΔR² • solid;   total HP R² • open)",
      y       = NULL,
      caption = paste0("Explained variance = ", round(dev_exp,0), "%;  R² = ",
                       round(r2_model,2),
                       if (!is.null(n_label)) paste0(";  n = ", n_label) else "")
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.line  = element_line(color="black"),
          panel.border = element_rect(color="black", fill=NA, linewidth=0.6))
}

p_mini <- plot_imp(mini_join, "Variable importance — Mini (nutrients)",
                   dev_exp = summary(gam_all_light_slim)$dev.expl*100,
                   r2_model = summary(gam_all_light_slim)$r.sq,
                   n_label = nrow(mod_all2))

p_land <- plot_imp(land_join, "Variable importance — Landscape",
                   dev_exp = summary(gam_no_prs)$dev.expl*100,
                   r2_model = summary(gam_no_prs)$r.sq,
                   n_label = nrow(mod_no_prs))

print(p_mini); print(p_land)

# ------------ Soil Bulk Characteristics ----------

glimpse(dfmL)

bulk_dat <- dfmL %>%
  select(functional_group, bulk_density_g_cm3) %>%
  filter(!is.na(bulk_density_g_cm3),
         !is.na(functional_group)) %>%
  mutate(
    functional_group = factor(
      functional_group,
      levels = c("Erect Shrub",
                 "Prostrate Shrub",
                 "Forb",
                 "Graminoid",
                 "Barren")
    )
  )


p_bulk_final <- ggplot(
  bulk_dat,
  aes(x = functional_group,
      y = bulk_density_g_cm3,
      fill = functional_group)
) +
  geom_boxplot(
    width = 0.6,
    alpha = 1,               
    color = "black",
    linewidth = 1,
    outlier.shape = NA
  ) +
  # raw points
  geom_jitter(
    width = 0.12,
    size = 2.2,
    alpha = 0.8,
    stroke = 0.4,
    shape = 21,
    color = "black"
  ) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  
  labs(
    x = NULL,
    y = expression("Bulk density (g cm"^-3*")")
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black", linewidth = 0.6),
    axis.text.x      = element_text(
      angle = 0,               
      hjust = 0.5,
      vjust = 0.5,
      face  = "bold"
    ),
    axis.text.y      = element_text(color = "black"),
    axis.title.y     = element_text(face = "bold")
  )

p_bulk_final


# ggsave("bulkdensity.png", p_bulk_final, width = 10, height = 4.9, units = "in", dpi = 300)


# ----------- Biomass stack --------

bio_comp <- dfmL %>%
  select(
    functional_group,
    bm_prop_shrub,
    bm_prop_graminoid,
    bm_prop_forb,
    bm_prop_moss_lichen,
    bm_prop_dead
  ) %>%
  filter(
    !is.na(functional_group)
  ) %>%
  filter(functional_group %in% c("Prostrate Shrub","Forb","Graminoid","Barren")) %>%
  mutate(
    functional_group = factor(
      functional_group,
      levels = c("Prostrate Shrub","Forb","Graminoid","Barren")
    )
  )

bio_comp_long <- bio_comp %>%
  pivot_longer(
    cols = starts_with("bm_prop_"),
    names_to = "component",
    values_to = "prop"
  ) %>%
  filter(!is.na(prop)) %>%
  mutate(
    component = recode(
      component,
      bm_prop_shrub        = "Shrub",
      bm_prop_graminoid    = "Graminoid",
      bm_prop_forb         = "Forb",
      bm_prop_moss_lichen  = "Moss/Lichen",
      bm_prop_dead         = "Standing dead"
    )
  ) %>%
  group_by(functional_group, component) %>%
  summarize(
    mean_prop = mean(prop, na.rm = TRUE),
    .groups   = "drop"
  )

# bio_comp_long %>% group_by(functional_group) %>% summarize(total = sum(mean_prop))

p_biomass_stack <- ggplot(
  bio_comp_long,
  aes(
    x = functional_group,
    y = mean_prop,
    fill = component
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.6,
    width = 0.7,
    alpha = 0.7   
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 1)
  ) +
  
  scale_fill_manual(
    values = c(
      "Shrub"         = grp_cols[["Prostrate Shrub"]],
      "Graminoid"     = grp_cols[["Graminoid"]],
      "Forb"          = grp_cols[["Forb"]],
      "Moss/Lichen"   = grp_cols[["Barren"]],
      "Standing dead" = "#4D4D4D"
    ),
    name = "Aboveground biomass\ncomposition"
  ) +
  
  labs(
    x = NULL,
    y = "Mean biomass composition (%)"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black", linewidth = 0.6),
    axis.text.x      = element_text(face = "bold", color = "black"),
    axis.text.y      = element_text(color = "black"),
    axis.title.y     = element_text(face = "bold"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold")
  )

p_biomass_stack

# ggsave("biomassComp.png", p_biomass_stack, width = 8, height = 4.9, units = "in", dpi = 300)

# ----------- PCV stack --------

pcv_comp <- dfmL %>%
  dplyr::select(
    functional_group,
    prop_Prostrate.Shrub,
    prop_graminoid,
    prop_forb,
    prop_Barren,
    prop_Erect.Shrub
  ) %>%
  filter(
    !is.na(functional_group)
  ) %>%
  filter(functional_group %in% c("Erect Shrub", "Prostrate Shrub","Forb","Graminoid","Barren")) %>%
  mutate(
    functional_group = factor(
      functional_group,
      levels = c("Erect Shrub","Prostrate Shrub","Forb","Graminoid","Barren")
    )
  )

pcv_comp_long <- pcv_comp %>%
  pivot_longer(
    cols = starts_with("prop_"),
    names_to = "component",
    values_to = "prop"
  ) %>%
  filter(!is.na(prop)) %>%
  mutate(
    component = recode(
      component,
      prop_Prostrate.Shrub        = "Prostrate Shrub",
      prop_graminoid    = "Graminoid",
      prop_forb         = "Forb",
      prop_Barren       = "Barren",
      prop_Erect.Shrub = "Erect Shrub")
  ) %>%
  group_by(functional_group, component) %>%
  summarize(
    mean_prop = mean(prop, na.rm = TRUE),
    .groups   = "drop"
  )

# bio_comp_long %>% group_by(functional_group) %>% summarize(total = sum(mean_prop))

p_stack <- ggplot(
  pcv_comp_long,
  aes(
    x = functional_group,
    y = mean_prop,
    fill = fct_reorder(component, mean_prop, .fun = min)
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.6,
    width = 0.7,
    alpha = 0.7
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 1)
  ) +
  
  scale_fill_manual(
    values = c(
      "Prostrate Shrub" = grp_cols[["Prostrate Shrub"]],
      "Erect Shrub"     = grp_cols[["Erect Shrub"]],
      "Graminoid"       = grp_cols[["Graminoid"]],
      "Forb"            = grp_cols[["Forb"]],
      "Barren"          = grp_cols[["Barren"]]
    ),
    name = "PCV"
  ) +
  
  labs(
    x = NULL,
    y = "Mean PCV composition (%)"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line       = element_line(color = "black"),
    axis.ticks      = element_line(color = "black", linewidth = 0.6),
    axis.text.x     = element_text(face = "bold", color = "black"),
    axis.text.y     = element_text(color = "black"),
    axis.title.y    = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(face = "bold")
  )


p_stack

# ggsave("PCVComp.png", p_stack, width = 8, height = 4.9, units = "in", dpi = 300)

# ----------- Light and Dark Figure (paper) ---------

# ---------- CH4: Light vs Dark paired Wilcoxon ----------
paired_ann_LD <- function(d, ycol) {
  d_collapsed <- d %>%
    dplyr::select(plot_key, functional_group, light_dark, !!rlang::sym(ycol)) %>%
    dplyr::filter(!is.na(.data[[ycol]])) %>%
    dplyr::group_by(plot_key, functional_group, light_dark) %>%
    dplyr::summarise(
      val = mean(.data[[ycol]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Pivot to wide format for paired test
  wide <- d_collapsed %>%
    dplyr::group_by(functional_group, plot_key) %>%
    dplyr::filter(dplyr::n_distinct(light_dark) == 2) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = light_dark, values_from = val)
  
  # Paired Wilcoxon test per functional group
  pv <- wide %>%
    dplyr::group_by(functional_group) %>%
    dplyr::summarise(
      n_pairs = sum(stats::complete.cases(Light, Dark)),
      p = if (n_pairs >= 2) stats::wilcox.test(Light, Dark, paired = TRUE)$p.value else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_lab = p_stars(p))
  
  # span and ymax for annotation
  span <- diff(range(d[[ycol]], na.rm = TRUE))
  ymax <- d %>%
    dplyr::group_by(functional_group) %>%
    dplyr::summarise(y = max(.data[[ycol]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(y = y + 0.10 * span)
  
  # Combine stats and annotation
  dplyr::left_join(pv, ymax, by = "functional_group") %>%
    dplyr::filter(!is.na(p)) %>%  # drop groups with no pairs
    dplyr::mutate(
      xmin = 0.9,
      xmax = 2.1,
      xmid = 1.5,
      span = span
    )
}

# ---------- CO2: GPP vs ER paired Wilcoxon ----------
paired_ann_ERGPP <- function(d_long) {
  d_collapsed <- d_long %>%
    dplyr::group_by(plot_key, functional_group, component) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  wide <- d_collapsed %>%
    tidyr::pivot_wider(names_from = component, values_from = value)
  
  # Paired Wilcoxon test per functional group
  pv <- wide %>%
    dplyr::group_by(functional_group) %>%
    dplyr::summarise(
      n_pairs = sum(stats::complete.cases(GPP, ER)),
      p = if (n_pairs >= 2) stats::wilcox.test(GPP, ER, paired = TRUE)$p.value else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_lab = p_stars(p))
  
  span <- diff(range(d_long$value, na.rm = TRUE))
  ymax <- d_long %>%
    dplyr::group_by(functional_group) %>%
    dplyr::summarise(y = max(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(y = y + 0.10 * span)
  
  dplyr::left_join(pv, ymax, by = "functional_group") %>%
    dplyr::filter(!is.na(p)) %>%  # drop groups with no pairs
    dplyr::mutate(
      xmin = 0.9,
      xmax = 2.1,
      xmid = 1.5,
      span = span
    )
}



theme_pub <- function(base_size = 10){
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(linewidth = 0.3, color = "grey85"),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = base_size + 2),
      legend.position  = "none",
      strip.text       = element_text(face = "bold"),
      panel.border     = element_rect(color = "grey85", fill = NA, linewidth = 0.8)
    )
}

# p value symbols
p_stars <- function(p) ifelse(p <= 0.001, "***",
                              ifelse(p <= 0.01,  "**",
                                     ifelse(p <= 0.05,  "*", "ns")))

dat <- df_master %>%
  mutate(light_dark = case_when(
    toupper(as.character(light_dark)) %in% c("L", "LIGHT") ~ "Light",
    toupper(as.character(light_dark)) %in% c("D", "DARK")  ~ "Dark",
    TRUE ~ as.character(light_dark)
  )) %>%
  # trim whitespace 
  mutate(functional_group = str_squish(functional_group)) %>%
  filter(!is.na(ch4_mgCH4_m2_d)) %>%
  filter(functional_group %in% names(grp_cols), light_dark %in% c("Light","Dark")) %>%
  mutate(functional_group = factor(functional_group, levels = names(grp_cols)),
         light_dark = factor(light_dark, levels = c("Light","Dark")))

# quick check
cat("CH4 dat rows:", nrow(dat), "\n")
print(table(dat$functional_group, dat$light_dark))


ann_ch4 <- paired_ann_LD(dat, "ch4_mgCH4_m2_d")

# keep only significant (p <= 0.05) 
ann_ch4_sig <- ann_ch4 %>%
  dplyr::filter(!is.na(p), p <= 0.05)

cat("paired groups (CH4):\n")
print(ann_ch4 %>% dplyr::select(functional_group, n_pairs, p, p_lab))
cat("significant groups (CH4, p <= 0.05):\n")
print(ann_ch4_sig %>% dplyr::select(functional_group, n_pairs, p, p_lab))
dat_co2 <- df_master %>%
  mutate(
    functional_group = str_squish(functional_group),
    GPP_gC_m2_d = -abs(GPP_gC_m2_d)   # changed to negative
  ) %>%
  filter(functional_group %in% names(grp_cols)) %>%
  select(plot_key, functional_group, GPP_gC_m2_d, ER_gC_m2_d) %>%
  distinct() %>%
  pivot_longer(
    cols = c(GPP_gC_m2_d, ER_gC_m2_d),
    names_to = "component",
    values_to = "value"
  ) %>%
  mutate(
    component = recode(component, GPP_gC_m2_d = "GPP", ER_gC_m2_d = "ER"),
    functional_group = factor(functional_group, levels = names(grp_cols)),
    component = factor(component, levels = c("GPP","ER"))
  ) %>%
  filter(!is.na(value)) %>% 
  filter(!(component == "ER" & value < 0))

cat("CO2 dat rows:", nrow(dat_co2), "\n")
print(table(dat_co2$functional_group, dat_co2$component))

# wilcoxon 
ann_co2 <- paired_ann_ERGPP(dat_co2)
ann_co2_sig <- ann_co2 %>% dplyr::filter(!is.na(p), p <= 0.05)

cat("paired groups (CO2):\n")
print(ann_co2 %>% dplyr::select(functional_group, n_pairs, p, p_lab))
cat("significant groups (CO2, p <= 0.05):\n")
print(ann_co2_sig %>% dplyr::select(functional_group, n_pairs, p, p_lab))

cat("\nHEAD dat (CH4):\n"); print(head(dat, 8))
cat("\nHEAD dat_co2 (CO2):\n"); print(head(dat_co2, 8))

# dat, ann_ch4, ann_ch4_sig, dat_co2, ann_co2, ann_co2_sig

ann_ch4_sig$y <- 5
ann_ch4_sig$span <- 12  
tick_height <- 0.3

ann_ch4_sig <- ann_ch4_sig  %>% 
  dplyr::mutate(
    y = 5,
    span = 12
  )


p_ch4 <- ggplot(
  dat,
  aes(x = light_dark, y = ch4_mgCH4_m2_d,
      fill = functional_group, alpha = light_dark)
) +
  # light gridlines
  geom_hline(yintercept = c(-6, -3, 0, 3,  6), colour = "grey85", linewidth = 0.5) +
  
  # boxplot
  geom_boxplot(width = 0.60, outlier.shape = NA, colour = "black", linewidth = 0.4) +
  
  # jitter
  geom_jitter(aes(color = functional_group), width = 0.12,
              shape = 16, size = 2, alpha = 0.70, show.legend = FALSE) +
  
  # significance bars
  geom_segment(data = ann_ch4_sig,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = ann_ch4_sig,
               aes(x = xmin, xend = xmin, y = y, yend = y - tick_height),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = ann_ch4_sig,
               aes(x = xmax, xend = xmax, y = y, yend = y - tick_height),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = ann_ch4_sig,
            aes(x = xmid, y = y + 0.02*span, label = p_lab),
            inherit.aes = FALSE, size = 4.2) +
  
  
  # scales
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_color_manual(values = grp_cols, guide = "none") +
  scale_alpha_manual(values = c(Light = 0.35, Dark = 1.00), guide = "none") +
  scale_x_discrete(limits = c("Light","Dark")) +
  scale_y_continuous(
    limits = c(-6, 6),
    breaks = c(-6, 0, 6),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  labs(x = NULL,
       y = expression("Methane flux ("*mg~CH[4]~m^{-2}~day^{-1}*")")) +
  
  facet_wrap(~ functional_group, nrow = 1, strip.position = "bottom") +
  theme_pub(13) +
  theme(
    panel.border    = element_blank(),  # remove the panel border
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, margin = ggplot2::margin(t = 6)),
    axis.title.y     = element_text(size = 14, margin = ggplot2::margin(r = 8))
  )


p_ch4



dat_light <- dat %>% 
  dplyr::filter(light_dark == "Light")
kruskal.test(ch4_mgCH4_m2_d ~ functional_group, data = dat_light)

library(FSA)

dunn_res <- dunnTest(
  ch4_mgCH4_m2_d ~ functional_group,
  data = dat_light,
  method = "bh"  
)

dunn_res

library(rstatix)

pairwise <- dat_light %>% 
  rstatix::dunn_test(ch4_mgCH4_m2_d ~ functional_group, p.adjust.method = "BH") %>% 
  rstatix::add_significance("p.adj")

pairwise_gpp <- dat_light %>% 
  rstatix::dunn_test(GPP_gC_m2_d ~ functional_group, p.adjust.method = "BH") %>% 
  rstatix::add_significance("p.adj")

pairwise_er <- dat_light %>% 
  rstatix::dunn_test(ER_gC_m2_d  ~ functional_group, p.adjust.method = "BH") %>% 
  rstatix::add_significance("p.adj")

pairwise_nee <- dat_light %>% 
  rstatix::dunn_test(co2_gC_m2_d  ~ functional_group, p.adjust.method = "BH") %>% 
  rstatix::add_significance("p.adj")

pairwise

pairwise <- pairwise %>% 
  mutate(
    sig_label = case_when(
      p.adj < 0.05 ~ "*",
      p.adj < 0.1  ~ ".",
      TRUE ~ ""
    )
  )

ann_fg_sig <- pairwise  %>% 
  mutate(
    functional_group = ifelse(group1 == "Forb", group2, group1),
    x = 1.5,
    y = -5.5,
    label = p.adj.signif
  )

ann_fg_sig2 <- ann_fg_sig  %>% 
  dplyr::filter(sig_label != "")

ann_fg_sig2 <- ann_fg_sig2  %>% 
  dplyr::mutate(
    functional_group = case_when(
      group1 == "Forb" ~ group2,
      group2 == "Forb" ~ group1,
      TRUE ~ functional_group
    ),
    x = 1.5,     # center between Light/Dark
    y = -5.5,    # bottom of panel
    xmin = 1.2,
    xmax = 1.8
  )

dat$functional_group <- factor(
  dat$functional_group,
  levels = c("Erect Shrub", "Prostrate Shrub", "Forb", "Graminoid", "Barren")
)

ann_fg_sig2 <- pairwise %>% 
  dplyr::filter(sig_label != "")  %>% 
  dplyr::mutate(
    functional_group = case_when(
      group1 == "Forb" ~ group2,
      group2 == "Forb" ~ group1
    )
  )

ann_fg_sig2 <- ann_fg_sig2  %>% 
  dplyr::mutate(
    x = 1.5,
    y = -5.5,
    xmin = 1.25,
    xmax = 1.75
  )

p_ch4 <- ggplot(
  dat,
  aes(x = light_dark, y = ch4_mgCH4_m2_d,
      fill = functional_group, alpha = light_dark)
) +
  # light gridlines
  geom_hline(yintercept = c(-6,-3, 0,3, 6), colour = "grey85", linewidth = 0.5) +
  
  # boxplot
  geom_boxplot(width = 0.60, outlier.shape = NA, colour = "black", linewidth = 0.4) +
  
  # jitter
  geom_jitter(aes(color = functional_group), width = 0.12,
              shape = 16, size = 2, alpha = 0.70, show.legend = FALSE) +
  
  # significance bars
  geom_segment(data = ann_ch4_sig,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = ann_ch4_sig,
               aes(x = xmin, xend = xmin, y = y, yend = y - tick_height),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = ann_ch4_sig,
               aes(x = xmax, xend = xmax, y = y, yend = y - tick_height),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = ann_ch4_sig,
            aes(x = xmid, y = y + 0.02*span, label = p_lab),
            inherit.aes = FALSE, size = 4.2) +
  
  geom_text(
    data = data.frame(
      functional_group = factor(
        "Forb",
        levels = levels(dat$functional_group)
      )
    ),
    aes(x = 1.5, y = -5.5, label = "*"),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  
  # scales
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_color_manual(values = grp_cols, guide = "none") +
  scale_alpha_manual(values = c(Light = 0.35, Dark = 1.00), guide = "none") +
  scale_x_discrete(limits = c("Light","Dark")) +
  scale_y_continuous(
    limits = c(-6, 6),
    breaks = c(-6,-3, 0,3, 6),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  labs(x = NULL,
       y = expression("Methane flux ("*mg~CH[4]~m^{-2}~day^{-1}*")")) +
  
  facet_wrap(~ functional_group, nrow = 1, strip.position = "bottom") +
  theme_pub(13) +
  theme(
    panel.border    = element_blank(),  # remove the panel border
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, margin = ggplot2::margin(t = 6)),
    axis.title.y     = element_text(size = 14, margin = ggplot2::margin(r = 8))
  )


p_ch4

# ggsave("ch4_fig.png", p_ch4, width = 16, height = 9, units = "in", dpi = 300)

# ------------ CO2 


p_co2 <- ggplot(
  dat_co2,
  aes(x = component, y = value,
      fill = functional_group, alpha = component)
) +
  # light gridlines
  geom_hline(yintercept = c( -10, 0, 5, 15, 25), colour = "grey85", linewidth = 0.5) +
  
  # boxplot
  geom_boxplot(width = 0.60, outlier.shape = NA, colour = "black", linewidth = 0.4) +
  
  # jitter
  geom_jitter(aes(color = functional_group), width = 0.12,
              shape = 16, size = 2, alpha = 0.70, show.legend = FALSE) +
  
  # significance bars (only p ≤ 0.05)
  # geom_segment(data = ann_co2_sig,
  #              aes(x = xmin, xend = xmax, y = y, yend = y),
  #              inherit.aes = FALSE, linewidth = 0.6) +
  # geom_segment(data = ann_co2_sig,
  #              aes(x = xmin, xend = xmin, y = y, yend = y - 0.02*span),
  #              inherit.aes = FALSE, linewidth = 0.6) +
  # geom_segment(data = ann_co2_sig,
  #              aes(x = xmax, xend = xmax, y = y, yend = y - 0.02*span),
  #              inherit.aes = FALSE, linewidth = 0.6) +
  # geom_text(data = ann_co2_sig,
  #           aes(x = xmid, y = y + 0.02*span, label = p_lab),
  #           inherit.aes = FALSE, size = 4.2) +
  
  # scales
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_color_manual(values = grp_cols, guide = "none") +
  scale_alpha_manual(values = c(GPP = 0.35, ER = 1.00), guide = "none") +
  scale_x_discrete(limits = c("GPP","ER")) +
  scale_y_continuous(breaks = c(-10, 0, 10, 20),
                     expand = expansion(mult = c(0.02, 0.18))) +
  
  labs(x = NULL,
       y = expression(C~flux~(mg~C~m^{-2}~day^{-1}))) +
  
  facet_wrap(~ functional_group, nrow = 1, strip.position = "bottom") +
  theme_pub(13) +
  theme(
    panel.border    = element_blank(),  # remove the panel border
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, margin = ggplot2::margin(t = 6)),
    axis.title.y     = element_text(size = 14, margin = ggplot2::margin(r = 8))
  )

p_co2

# ---- NEE


dat_NEE <- df_master %>%
  mutate(functional_group = str_squish(functional_group)) %>%
  filter(functional_group %in% names(grp_cols)) %>%
  select(plot_key, functional_group, co2_gC_m2_d, ER_gC_m2_d) %>%
  distinct() %>%
  pivot_longer(
    cols = c(co2_gC_m2_d, ER_gC_m2_d),
    names_to = "component",
    values_to = "value"
  ) %>%
  mutate(component = recode(component, co2_gC_m2_d = "NEE", ER_gC_m2_d = "ER"),
         functional_group = factor(functional_group, levels = names(grp_cols)),
         component = factor(component, levels = c("NEE","ER"))) %>%
  filter(!is.na(value)) %>% 
  filter(!(component == "ER" & value < 0))


dat_NEE <- df_master %>%
  mutate(functional_group = str_squish(functional_group)) %>%
  filter(functional_group %in% names(grp_cols)) %>%
  select(plot_key, functional_group, light_dark,
         co2_gC_m2_d, ER_gC_m2_d) %>%
  distinct() %>%
  pivot_longer(
    cols = c(co2_gC_m2_d, ER_gC_m2_d),
    names_to = "component",
    values_to = "value"
  ) %>%
  mutate(
    component = recode(component,
                       co2_gC_m2_d = "NEE",
                       ER_gC_m2_d   = "ER"),
    functional_group = factor(functional_group,
                              levels = names(grp_cols)),
    component = factor(component,
                       levels = c("NEE", "ER"))
  ) %>%
  filter(!is.na(value)) %>% 
  filter(!(component == "NEE" & light_dark != "L")) %>%
  filter(!(component == "ER" & value < 0))


# ann_NEE <- paired_ann_ERGPP(dat_NEE)
# ann_NEE_sig <- ann_NEE %>% dplyr::filter(!is.na(p), p <= 0.05)

p_NEE <- ggplot(
  dat_NEE,
  aes(x = component, y = value,
      fill = functional_group, alpha = component)
) +
  # light gridlines
  geom_hline(yintercept = c(-10, 0, 5, 10, 15), colour = "grey85", linewidth = 0.5) +
  
  # boxplot
  geom_boxplot(width = 0.60, outlier.shape = NA, colour = "black", linewidth = 0.4) +
  
  # jitter
  geom_jitter(aes(color = functional_group), width = 0.12,
              shape = 16, size = 2, alpha = 0.70, show.legend = FALSE) +
  
  # significance bars (commented out)
  # geom_segment(data = ann_NEE_sig,
  #              aes(x = xmin, xend = xmax, y = y, yend = y),
  #              inherit.aes = FALSE, linewidth = 0.6) +
  # geom_segment(data = ann_NEE_sig,
  #              aes(x = xmin, xend = xmin, y = y, yend = y - 0.02*span),
  #              inherit.aes = FALSE, linewidth = 0.6) +
  # geom_segment(data = ann_NEE_sig,
  #              aes(x = xmax, xend = xmax, y = y, yend = y - 0.02*span),
  #              inherit.aes = FALSE, linewidth = 0.6) +
  # geom_text(data = ann_NEE_sig,
  #           aes(x = xmid, y = y + 0.02*span, label = p_lab),
  #           inherit.aes = FALSE, size = 4.2) +
  
  # scales
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_color_manual(values = grp_cols, guide = "none") +
  scale_alpha_manual(values = c(NEE = 0.35, ER = 1.00), guide = "none") +
  scale_x_discrete(limits = c("NEE","ER")) +
  scale_y_continuous(breaks = c(-5, 5, 15),
                     expand = expansion(mult = c(0.02, 0.18))) +
  
  labs(x = NULL,
       y = expression(C~flux~("g C m"^{-2}~day^{-1}))) +
  
  facet_wrap(~ functional_group, nrow = 1, strip.position = "bottom") +
  theme_pub(13) +
  theme(
    panel.border    = element_blank(),  # remove the panel border
    strip.placement  = "outside",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, margin = ggplot2::margin(t = 6)),
    axis.title.y     = element_text(size = 14, margin = ggplot2::margin(r = 8))
  )


p_NEE

stack_two_panels <- function(top_plot, bottom_plot,
                             pad_left = 0.075,
                             pad_right = 0.005,
                             pad_bottom = 0.115,
                             pad_top = 0.015,
                             lwd = 0.9) {
  
  add_panel_border <- function(p) {
    ggdraw() +
      draw_plot(p + theme(panel.border = element_blank())) +
      draw_line(
        x = c(pad_left, 1 - pad_right, 1 - pad_right, pad_left, pad_left),
        y = c(pad_bottom, pad_bottom, 1 - pad_top, 1 - pad_top, pad_bottom),
        size = 1, color = "grey85"
      )
  }
  
  # apply border
  top <- add_panel_border(top_plot)
  bottom <- add_panel_border(bottom_plot)
  top / bottom +
    plot_annotation(
      tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold", size = 14),
        plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)  # <- specify t,r,b,l
      )
    )
}


final_two_CO2 <- stack_two_panels(p_ch4, p_co2)
final_two_CO2

final_two_NEE <- stack_two_panels(p_ch4, p_NEE)
final_two_NEE

# ggsave("GPP_ch4_fig.png", final_two_CO2, width = 10, height = 9.8, units = "in", dpi = 300)
# ggsave("NEE_ch4_fig.png", final_two_NEE, width = 10, height = 9.8, units = "in", dpi = 300)


## ------- Kenzies 3 panel idea (paper)


dat_GPP <- dat_co2 %>% filter(component == "GPP")
dat_ER  <- dat_co2 %>% filter(component == "ER")
dat_NEE_plot <- dat_NEE %>% filter(component == "NEE")


p_ch4_A <- ggplot(
  dat,
  aes(x = light_dark, y = ch4_mgCH4_m2_d,
      fill = functional_group, alpha = light_dark)
) +
  
  geom_hline(yintercept = c(-10, -5, 0, 5),
             colour = "grey85", linewidth = 0.5) +
  
  # Boxplots
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    colour = "black",
    linewidth = 0.40
  ) +
  
  # Points
  geom_jitter(
    aes(color = functional_group),
    width = 0.10,
    size = 2,
    alpha = 0.70,
    show.legend = FALSE
  ) +
  
  # Significance bars
  geom_segment(
    data = ann_ch4_sig,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_segment(
    data = ann_ch4_sig,
    aes(x = xmin, xend = xmin, y = y, yend = y - 0.02*span),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_segment(
    data = ann_ch4_sig,
    aes(x = xmax, xend = xmax, y = y, yend = y - 0.02*span),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = ann_ch4_sig,
    aes(x = xmid, y = y + 0.02*span, label = p_lab),
    inherit.aes = FALSE,
    size = 4.2
  ) +
  
  # Scales
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_color_manual(values = grp_cols, guide = "none") +
  scale_alpha_manual(values = c(Light = 0.35, Dark = 1.00), guide = "none") +
  scale_x_discrete(limits = c("Light","Dark")) +
  scale_y_continuous(
    breaks = c(-10,-5, 0, 5),
    expand = expansion(mult = c(0.02, 0.25))  # MATCH expansion with GPP/ER
  ) +
  
  # Labels
  labs(
    x = NULL,
    y = expression("Methane flux ("*mg~CH[4]~m^{-2}~day^{-1}*")")
  ) +
  
  # Facets
  facet_wrap(
    ~ functional_group,
    nrow = 1,
    strip.position = "bottom"
  ) +
  
  # Theme
  theme_pub(13) +
  theme(
    panel.border = element_blank(),
    
    #  panel widths
    panel.spacing.x = unit(0.35, "lines"),
    panel.spacing.y = unit(0.55, "lines"),
    
    strip.placement  = "outside",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", margin = ggplot2::margin(t = 6)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 10)),
    axis.text.y  = element_text(margin = ggplot2::margin(r = -1)),
    plot.margin = ggplot2::margin(t = 10, r = 15, b = -1, l = 10)
  )


make_co2_panel <- function(data, breaks, ylabel, alpha_val = 1) {
  ggplot(
    data,
    aes(x = functional_group, y = value,
        fill = functional_group)
  ) +
    geom_hline(yintercept = breaks, colour = "grey85", linewidth = 0.5,
               alpha = alpha_val) +
    geom_boxplot(width = 0.60, outlier.shape = NA,
                 colour = "black", linewidth = 0.4,
                 alpha = alpha_val) +
    geom_jitter(aes(color = functional_group),
                width = 0.15, size = 2, alpha = 0.7) +
    scale_fill_manual(values = grp_cols, guide = "none") +
    scale_color_manual(values = grp_cols, guide = "none") +
    # scale_alpha_manual(values = c("NEE" = 0.35, "ER" = 1.00), guide = "none") +
    # scale_y_continuous(breaks = c(-10, 0, 10, 20),
    #                    expand = expansion(mult = c(0.02, 0.18))) +
    labs(x = NULL, y = ylabel) +
    theme_pub(13) +
    theme(
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", margin = ggplot2::margin(t = 6)),
      axis.text.x = element_text(face = "bold", colour = "black"),
      axis.title.y = element_text(
        margin = ggplot2::margin(t = 0, r = 8, b = 0, l = 0)
      )
    )
}

stack_three_panels <- function(p1, p2, p3,
                               pad_left = 0.075,
                               pad_right = 0.005,
                               pad_bottom = 0.085,
                               pad_top = 0.015,
                               lwd = 0.9) {
  
  add_border <- function(p) {
    ggdraw() +
      draw_plot(p + theme(panel.border = element_blank())) +
      draw_line(
        x = c(pad_left, 1 - pad_right, 1 - pad_right, pad_left, pad_left),
        y = c(pad_bottom, pad_bottom, 1 - pad_top, 1 - pad_top, pad_bottom),
        size = lwd, color = "grey85"
      )
  }
  
  p1b <- add_border(p1)
  p2b <- add_border(p2)
  p3b <- add_border(p3)
  
  p1b / p2b / p3b +
    plot_annotation(
      tag_levels = "a",
      tag_prefix = "(",
      tag_suffix = ")",
      theme = theme(
        plot.tag = element_text(face = "bold", size = 16),
        plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
      )
    )
}

p_GPP <- make_co2_panel(dat_GPP,
                        breaks = c(-5, 0, 10),
                        ylabel = expression(GPP~(g~C~m^{-2}~day^{-1})),
                        alpha_val = 0.5)

p_ER <- make_co2_panel(dat_ER,
                       breaks = c(-10, 0, 10, 20),
                       ylabel = expression(ER~(g~C~m^{-2}~day^{-1})))

p_NEE_only <- make_co2_panel(
  dat_NEE_plot,
  breaks = c(-10, 0, 10),
  ylabel = expression(NEE~(g~C~m^{-2}~day^{-1})),
  alpha_val = 1
)

p_GPP
# ggsave("p_NEE_only_fig.png", p_NEE_only, width = 16, height = 9, units = "in", dpi = 300)

fig_CH4_NEE_ER <- stack_three_panels(
  p_ch4_A,
  p_NEE_only,
  p_ER
)

fig_CH4_NEE_ER

# ggsave("NEE_ch4_3stack_fig_big.png", fig_CH4_NEE_ER, width = 12, height = 16, units = "in", dpi = 300)

fig_CH4_GPP_ER <- stack_three_panels(
  p_ch4_A,
  p_GPP,
  p_ER
)

fig_CH4_GPP_ER

# ggsave("GPP_ch4_3stack_fig.png", fig_CH4_GPP_ER, width = 12, height = 16, units = "in", dpi = 300)

# ------------ summary stats ----------------


df_shrubtest <- dfmL %>%
  mutate(
    shrub_group = case_when(
      functional_group %in% c("Erect Shrub", "Prostrate Shrub") ~ "Shrub",
      functional_group %in% c("Graminoid", "Forb", "Barren")    ~ "No shrub",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(shrub_group),
    !is.na(ch4_mgCH4_m2_d)
  )

table(df_shrubtest$shrub_group)

flux_shrub    <- df_shrubtest %>% filter(shrub_group == "Shrub")    %>% pull(ch4_mgCH4_m2_d)
flux_noshrub  <- df_shrubtest %>% filter(shrub_group == "No shrub") %>% pull(ch4_mgCH4_m2_d)

wilcox_res <- wilcox.test(flux_shrub, flux_noshrub, exact = FALSE)

wilcox_res

summaries <- df_shrubtest %>%
  group_by(shrub_group) %>%
  summarize(
    n        = n(),
    median   = median(ch4_mgCH4_m2_d, na.rm = TRUE),
    IQR_low  = quantile(ch4_mgCH4_m2_d, 0.25, na.rm = TRUE),
    IQR_high = quantile(ch4_mgCH4_m2_d, 0.75, na.rm = TRUE)
  )

summaries

df_shrubplot <- dfmL %>%
  mutate(
    shrub_group = case_when(
      functional_group %in% c("Erect Shrub", "Prostrate Shrub") ~ "Shrub",
      functional_group %in% c("Graminoid", "Forb", "Barren")    ~ "No shrub",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(shrub_group),
    !is.na(ch4_mgCH4_m2_d)
  )

df_shrubplot$shrub_group <- factor(df_shrubplot$shrub_group,
                                   levels = c("Shrub","No shrub"))

p_shrub_effect <- ggplot(df_shrubplot,
                         aes(x = shrub_group,
                             y = ch4_mgCH4_m2_d,
                             fill = shrub_group)) +
  
  geom_boxplot(width = 0.5,
               alpha = 0.6,
               color = "black",
               outlier.shape = NA,
               linewidth = 0.8) +
    geom_jitter(width = 0.12,
              size = 2,
              alpha = 0.4,
              stroke = 0.4,
              color = "black",
              aes(color = shrub_group)) +
    geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey40",
             linewidth = 0.6) +
  
  scale_fill_manual(
    values = c(
      "Shrub"     = grp_cols[["Prostrate Shrub"]],
      "No shrub"  = grp_cols[["Graminoid"]]
    ),
    guide = "none"
  ) +
  scale_color_manual(
    values = c(
      "Shrub"     = grp_cols[["Prostrate Shrub"]],
      "No shrub"  = grp_cols[["Graminoid"]]
    ),
    guide = "none"
  ) +
  
  labs(
    x = NULL,
    y = expression(CH[4]*" flux (mg C m"^{-2}*" d"^{-1}*")"),
    title = ""
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.text.x      = element_text(face = "bold", color = "black"),
    axis.text.y      = element_text(color = "black"),
    axis.title.y     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )

p_shrub_effect


# ------- results stats------

df_master %>%
  group_by(functional_group) %>%
  summarise(
    n = n(),
    mean_ch4 = mean(ch4_mgCH4_m2_d, na.rm = TRUE),
    sd_ch4 = sd(ch4_mgCH4_m2_d, na.rm = TRUE),
    median_ch4 = median(ch4_mgCH4_m2_d, na.rm = TRUE),
    iqr_ch4 = IQR(ch4_mgCH4_m2_d, na.rm = TRUE)
  )


df_master %>%
  filter(light_dark == "L") %>%
  group_by(functional_group) %>%
  summarise(
    n = n(),
    mean_ch4 = mean(ch4_mgCH4_m2_d, na.rm = TRUE),
    sd_ch4 = sd(ch4_mgCH4_m2_d, na.rm = TRUE),
    median_ch4 = median(ch4_mgCH4_m2_d, na.rm = TRUE),
    iqr_ch4 = IQR(ch4_mgCH4_m2_d, na.rm = TRUE)
  )

df_master %>%
  filter(light_dark == "D") %>%
  group_by(functional_group) %>%
  summarise(
    n = n(),
    mean_ch4 = mean(ch4_mgCH4_m2_d, na.rm = TRUE),
    sd_ch4 = sd(ch4_mgCH4_m2_d, na.rm = TRUE),
    median_ch4 = median(ch4_mgCH4_m2_d, na.rm = TRUE),
    iqr_ch4 = IQR(ch4_mgCH4_m2_d, na.rm = TRUE)
  )

df_master %>%
  group_by(light_dark) %>%
  summarise(
    n = n(),
    mean_ch4 = mean(ch4_mgCH4_m2_d, na.rm = TRUE),
    sd_ch4 = sd(ch4_mgCH4_m2_d, na.rm = TRUE)
  )

wilcox.test(ch4_mgCH4_m2_d ~ light_dark, data = df_master)


df_light <- df_master %>% filter(light_dark == "L")

cor.test(df_light$ch4_mgCH4_m2_d, df_light$soil_moisture_avg, use = "complete.obs")
cor.test(df_light$ch4_mgCH4_m2_d, df_light$air_temp, use = "complete.obs")


wilcox_by_group <- df_master %>%
  group_by(functional_group) %>%
  wilcox_test(
    ch4_mgCH4_m2_d ~ light_dark,
    detailed = TRUE
  ) %>%
  adjust_pvalue(method = "BH") %>%   # FDR correction
  add_significance()

wilcox_by_group


wilcox_table <- wilcox_by_group %>%
  select(
    functional_group,
    p
  ) %>%
  rename(
    `Functional Group` = functional_group,
    `p-value` = p
  ) %>%
  gt() %>%
  
  tab_header(
    title = "Light vs Dark Methane Flux Differences",
    subtitle = "Wilcoxon rank-sum test by functional group"
  ) %>%
  
  fmt_number(
    columns = `p-value`,
    decimals = 3
  ) %>%
  
  cols_align(
    align = "center",
    columns = `Functional Group`
  ) %>%
  
  tab_source_note(
    source_note = "* P-values from Wilcoxon rank-sum test"
  ) 

for (fg in names(grp_cols)) {
  wilcox_table <- wilcox_table %>%
    tab_style(
      style = cell_fill(color = alpha(grp_cols[[fg]], 0.15)),
      locations = cells_body(rows = `Functional Group` == fg)
    )
}

wilcox_table


df_clean <- df_master %>%
  mutate(
    light_dark = case_when(
      toupper(light_dark) %in% c("L","LIGHT") ~ "Light",
      toupper(light_dark) %in% c("D","DARK")  ~ "Dark",
      TRUE ~ light_dark
    ),
    functional_group = str_squish(functional_group)
  ) %>%
  filter(light_dark %in% c("Light","Dark")) %>%
  group_by(plot_key, functional_group, light_dark) %>%
  summarise(
    ch4 = mean(ch4_mgCH4_m2_d, na.rm = TRUE),
    co2 = mean(co2_gC_m2_d, na.rm = TRUE),
    .groups = "drop"
  )

df_wide <- df_clean %>%
  pivot_wider(names_from = light_dark, values_from = c(ch4, co2))

df_wide <- df_wide %>%
  mutate(
    GPP = co2_Light - co2_Dark,
    ER  = co2_Dark
  )

p_ch4 <- df_wide %>%
  group_by(functional_group) %>%
  summarise(
    p = wilcox.test(ch4_Light, ch4_Dark, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(variable = "CH4")

p_nee <- df_wide %>%
  group_by(functional_group) %>%
  summarise(
    p = wilcox.test(co2_Light, co2_Dark, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(variable = "NEE")

p_gpp <- df_wide %>%
  group_by(functional_group) %>%
  summarise(
    p = wilcox.test(GPP, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  mutate(variable = "GPP")

p_er <- df_wide %>%
  group_by(functional_group) %>%
  summarise(
    p = wilcox.test(ER, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  mutate(variable = "ER")

p_values <- bind_rows(p_ch4, p_nee, p_gpp, p_er)

p_values


p_values_wide <- p_values %>%
  select(functional_group, variable, p) %>%
  pivot_wider(names_from = functional_group, values_from = p) %>%
  select(variable, all_of(fg_order))  # reorder columns
plot_area <- 0.25
n_plots   <- 100

df_master %>%
  summarise(total_mg_day =
              mean(ch4_mgCH4_m2_d, na.rm = TRUE) *
              plot_area * n_plots)



df_master %>%
  summarise(net_flux = mean(ch4_mgCH4_m2_d, na.rm = TRUE))

df_master %>%
  summarise(total_net = sum(ch4_mgCH4_m2_d, na.rm = TRUE))

df_master %>%
  group_by(plot_key) %>%
  summarise(cumulative = sum(ch4_mgCH4_m2_d, na.rm = TRUE))

mean(df_master$ch4_mgCH4_m2_d)
sd(df_master$ch4_mgCH4_m2_d)



effect_sizes <- df_master %>%
  group_by(functional_group) %>%
  wilcox_effsize(
    ch4_mgCH4_m2_d ~ light_dark
  )

effect_sizes


kruskal.test(
  ch4_mgCH4_m2_d ~ interaction(functional_group, light_dark),
  data = df_master
)


# Select nutrients + functional group, filter complete cases
prs_vars <- c("no3_n", "nh4_n", "ca", "mg", "k", "p", "fe", "mn", "cu", "zn", "b", "s", "cd", "al")
df_prs <- df_master %>%
  dplyr::select(functional_group, all_of(prs_vars)) %>%
  filter(complete.cases(.))  # only plots with data

# Summarize: mean nutrient per functional group
prs_summary <- df_prs %>%
  group_by(functional_group) %>%
  summarise(
    across(all_of(prs_vars), 
           ~ paste0(round(mean(.x, na.rm = TRUE), 2), " ± ", round(sd(.x, na.rm = TRUE), 2))),
    .groups = "drop"
  )

chem_labels <- c(
  no3_n = "NO₃⁻",
  nh4_n = "NH₄⁺",
  ca    = "Ca²⁺",
  mg    = "Mg²⁺",
  k     = "K⁺",
  p     = "H₂PO₄⁻",
  fe    = "Fe³⁺",
  mn    = "Mn²⁺",
  cu    = "Cu²⁺",
  zn    = "Zn²⁺",
  b     = "B³⁺",
  s     = "SO₄²⁻",
  cd    = "Pb²⁺", 
  al    = "Al³⁺"
)


prs_long <- prs_summary %>%
  pivot_longer(
    -functional_group, 
    names_to = "nutrient", 
    values_to = "mean_value"
  ) %>%
  mutate(nutrient = recode(nutrient, !!!chem_labels))  

# order
fg_order <- c("Erect Shrub", "Prostrate Shrub", "Forb", "Graminoid", "Barren")

# Pivot wider
prs_table <- prs_long %>%
  pivot_wider(names_from = functional_group, values_from = mean_value)

prs_table <- prs_table %>%
  select(nutrient, all_of(fg_order))

# Create GT table with per-column coloring
prs_table %>%
  gt(rowname_col = "nutrient") %>%
  tab_header(
    title = "PRS Nutrient Content by Functional Group",
    subtitle = "Mean concentration across 40 plots"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 2
  ) 

light_model_nuts <- c("NO₃⁻", "NH₄⁺", "Cu²⁺", "H₂PO₄⁻")
dark_model_nuts  <- c( "SO₄²⁻", "Mn²⁺")

#  table
prs_table %>%
  gt(rowname_col = "nutrient") %>%
  tab_header(
    title = "",
    subtitle = "* Bold indicates nutrient used in Models; grey indicates Dark Model-only additions"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 2
  ) %>%
  # Highlight Light Model nutrients
  tab_style(
    style = cell_fill(color = "grey95"),
    locations = cells_body(
      rows = nutrient %in% light_model_nuts
    )
  ) %>%
  # Highlight Dark Model nutrients (different color)
  tab_style(
    style = cell_fill(color = "grey85"),
    locations = cells_body(
      rows = nutrient %in% dark_model_nuts
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = nutrient %in% union(light_model_nuts, dark_model_nuts)
    )
  )


prs_table %>%
  gt(rowname_col = "nutrient") %>%
  tab_header(
    title = ""
  ) %>%
  tab_source_note(  # subtitle at the bottom
    source_note = md("*Grey indicates nutrient used in Models; Dark grey indicates Dark Model-only additions*")
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 2
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_fill(color = "grey95"),
    locations = cells_body(
      rows = nutrient %in% light_model_nuts
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "grey85"),
    locations = cells_body(
      rows = nutrient %in% setdiff(dark_model_nuts, light_model_nuts)
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = nutrient %in% union(light_model_nuts, dark_model_nuts)
    )
  )

prs_gt <- prs_table %>%
  gt(rowname_col = "nutrient") %>%
  tab_header(
    title = ""  
  ) %>%
  fmt_number(columns = everything(), decimals = 2) %>%
  # Highlight Light Model 
  tab_style(
    style = cell_fill(color = "grey95"),
    locations = cells_body(rows = nutrient %in% light_model_nuts)
  ) %>%
  # Highlight Dark Model 
  tab_style(
    style = cell_fill(color = "grey85"),
    locations = cells_body(rows = nutrient %in% setdiff(dark_model_nuts, light_model_nuts))
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = nutrient %in% union(light_model_nuts, dark_model_nuts))
  ) %>%
  # Center the column headers
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  { 
    tmp <- .
    for (fg in names(grp_cols)) {
      tmp <- tmp %>%
        tab_style(
          style = cell_fill(color = alpha(grp_cols[[fg]], 0.2)),
          locations = cells_column_labels(columns = fg)
        )
    }
    tmp
  } %>%
  tab_source_note(
    source_note = "* Bold indicates nutrient used in Models; grey indicates Dark Model-only additions"
  )

prs_gt

# gtsave(prs_gt, "prs_table.png")


# ------------------------ PH panel plot (kenzie) --------------------

vars_box <- c(
  "C_N_ratio",
    "total_n", "cu", 
  "p",
  "p_h",
  "soil_moisture_avg",
  "x12cm_soil_temp",
  "gcc_mean"
)


df_box <- dfmL %>%
  select(functional_group, all_of(vars_box)) %>%
  pivot_longer(
    cols = -functional_group,
    names_to = "variable",
    values_to = "value"
  )

var_labs <- c(
  C_N_ratio = "Soil C/N",
  total_n = "Total N (µg 10 cm⁻²)",
  cu = "Cu (µg 10 cm⁻²)",
  p = "P (µg 10 cm⁻²)",
  p_h = "Soil pH",
  soil_moisture_avg = "Soil Moisture (%)",
  x12cm_soil_temp = "Soil Temp (°C, 12cm)",
  gcc_mean = "GCC"
)


df_box <- df_box %>%
  mutate(variable_lab = factor(variable, levels = vars_box, labels = var_labs))

df_box <- df_box %>%
  mutate(
    functional_group = factor(
      functional_group,
      levels = c(
        "Erect Shrub",
        "Prostrate Shrub",
        "Forb",
        "Graminoid",
        "Barren"
      )
    )
  )


p_box <- ggplot(
  df_box,
  aes(x = functional_group, y = value, fill = functional_group)
) +
  
  geom_boxplot(
    outlier.alpha = 0.15,
    linewidth = 0.4,
    alpha = 0.65
  ) +
  
  geom_jitter(
    aes(color = functional_group),
    width = 0.12,
    shape = 16,
    size = 1.5,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  
  stat_summary(
    fun = median,
    geom = "point",
    size = 1.4,
    color = "black"
  ) +
  
  scale_fill_manual(values = grp_cols, name = "Vegetation Type") +
  scale_color_manual(values = grp_cols, guide = "none") +
  
  facet_wrap(
    ~ variable_lab,
    scales = "free_y",
    ncol = 4
  ) +
  
  ylab(NULL) +
  labs(x = "Vegetation type") +
  guides(
    fill = guide_legend(override.aes = list(size = 1))
  ) +

  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "grey92"),
    axis.title = element_text(size =18),
    axis.text = element_text(size = 15),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    strip.text = element_text(size = 14), face = "bold")
    
  
# ggsave(
#   "raw_summary_by_vegetation.png",
#   p_box,
#   width = 17,
#   height = 9,
#   dpi = 300
# )


# -------------- Summaries pretty table (paper) --------

df_flux_summary <- df_master %>%
  mutate(light_dark = recode(light_dark, L = "light", D = "dark")) %>%
  group_by(plot_key, functional_group) %>%
  summarise(
    CH4_light = mean(ch4_mgCH4_m2_d[light_dark=="light"], na.rm = TRUE),
    CH4_dark  = mean(ch4_mgCH4_m2_d[light_dark=="dark"], na.rm = TRUE),
    CH4_total = sum(ch4_mgCH4_m2_d, na.rm = TRUE),
    
    GPP_mean = mean(GPP_gC_m2_d, na.rm = TRUE),
    ER_mean  = mean(ER_gC_m2_d, na.rm = TRUE),
    NEE      = mean(co2_gC_m2_d, na.rm = TRUE),
    
    total_C  = mean(total_C, na.rm = TRUE),
    total_N  = mean(total_N, na.rm = TRUE),
    
    C_N_ratio = mean(total_C, na.rm = TRUE) / mean(total_N, na.rm = TRUE),
    inorg_C = mean(inorg_C, na.rm = TRUE),
    LOI_mean = mean(LOI_mean, na.rm = TRUE),
    soil_depth_avg = mean(soil_depth_avg, na.rm = TRUE),
    x12cm_soil_temp = mean(x12cm_soil_temp, na.rm = TRUE),
    soil_moisture_avg = mean(soil_moisture_avg, na.rm = TRUE),
    volumetric_water_content_percent = mean(volumetric_water_content_percent, na.rm = TRUE),
    
    
    across(c(total_n, no3_n, nh4_n, ca, mg, k, p, fe, mn, cu, zn, b, s, cd, al),
           ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    
    .groups = "drop"
  )


vars_for_test <- df_flux_summary %>%
  select(-plot_key, -functional_group) %>%
  names()

kw_pvals <- map_dfr(vars_for_test, function(v){
  
  f <- as.formula(paste(v, "~ functional_group"))
  
  res <- kruskal.test(f, data = df_flux_summary)
  
  tibble(
    variable = v,
    p_value  = res$p.value
  )
})


df_fg_summary <- df_flux_summary %>%
  group_by(functional_group) %>%
  summarise(
    across(-plot_key,
           ~ paste0(round(mean(.x, na.rm = TRUE), 2), " ± ", round(sd(.x, na.rm = TRUE), 2))),
    .groups = "drop"
  )


df_fg_summary_1 <- df_flux_summary %>%
  # group_by(functional_group) %>%
  summarise(
    across(-plot_key,
           ~ paste0(round(mean(.x, na.rm = TRUE), 2), " ± ", round(sd(.x, na.rm = TRUE), 2))),
    .groups = "drop"
  )


chem_labels <- c(
  volumetric_water_content_percent = "VWC (cm³/cm³)",
  soil_depth_avg = "Soil Depth (cm)",
  x12cm_soil_temp = " Soil Temp (12cm)",
  soil_moisture_avg = "Soil Moisture (%)",
  total_n = "Total N",
  no3_n = "NO₃⁻",
  nh4_n = "NH₄⁺",
  ca    = "Ca²⁺",
  mg    = "Mg²⁺",
  k     = "K⁺",
  p     = "H₂PO₄⁻",
  fe    = "Fe³⁺",
  mn    = "Mn²⁺",
  cu    = "Cu²⁺",
  zn    = "Zn²⁺",
  b     = "B³⁺",
  s     = "SO₄²⁻",
  cd    = "Pb²⁺",
  al    = "Al³⁺",
  total_C = "%C",
  total_N = "%N",
  C_N_ratio = "C/N",
  inorg_C = "Inorganic C (%)",
  LOI_mean = "LOI (%)",
  CH4_light = "CH₄ (Light, mg m⁻² d⁻¹)",
  CH4_dark  = "CH₄ (Dark, mg m⁻² d⁻¹)",
  CH4_total = "CH₄ (Total, mg m⁻² d⁻¹)",
  GPP_mean  = "GPP (g C m⁻² d⁻¹)",
  ER_mean   = "ER (g C m⁻² d⁻¹)",
  NEE       = "NEE (g C m⁻² d⁻¹)"
)

kw_pvals <- kw_pvals %>%
  mutate(variable = recode(variable, !!!chem_labels))


df_long <- df_fg_summary %>%
  pivot_longer(-functional_group, names_to = "variable", values_to = "mean_value") %>%
  mutate(variable = recode(variable, !!!chem_labels))


fg_order <- c("Erect Shrub", "Prostrate Shrub", "Forb", "Graminoid", "Barren")

df_table <- df_long %>%
  pivot_wider(names_from = functional_group, values_from = mean_value) %>%
  dplyr::select(variable, all_of(fg_order))

df_table <- df_table %>%
  mutate(`Erect Shrub` = ifelse(variable == "LOI (%)" & `Erect Shrub` == "NaN ± NA", "", `Erect Shrub`))
df_table <- df_table %>%
  mutate(`Erect Shrub` = ifelse(variable == "VWC (cm³/cm³)" & `Erect Shrub` == "NaN ± NA", "", `Erect Shrub`))


light_model_nuts <- c("Total N", "Cu²⁺", "H₂PO₄⁻")
dark_model_nuts  <- c("SO₄²⁻", "Mn²⁺")

df_table <- df_table %>%
  left_join(kw_pvals, by = "variable") %>%
  mutate(
    p_value = ifelse(is.na(p_value), "",
                     format.pval(p_value, digits = 3, eps = 0.001))
  )


gt_table <- df_table %>%
  gt(rowname_col = "variable") %>%
  cols_label(p_value = "p-value") %>%
  # tab_header(title = "Functional Group Summary") %>%
  fmt_number(columns = everything(), decimals = 2) %>%
  # Highlight model nutrients
  # tab_style(style = cell_fill(color = "grey95"),
  #           locations = cells_body(rows = variable %in% light_model_nuts)) %>%
  # tab_style(style = cell_fill(color = "grey85"),
  #           locations = cells_body(rows = variable %in% setdiff(dark_model_nuts, light_model_nuts))) %>%
  # tab_style(style = cell_text(weight = "bold"),
  #           locations = cells_body(rows = variable %in% union(light_model_nuts, dark_model_nuts))) %>%
  # Add borders after methane
  tab_style(
    style = cell_borders(sides = "bottom", weight = px(2), color = "black"),
    locations = cells_body(rows = variable == "CH₄ (Total, mg m⁻² d⁻¹)")
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", weight = px(2), color = "black"),
    locations = cells_body(rows = variable == "NEE (g C m⁻² d⁻¹)")
  ) %>%
  # Add borders after LOI/C/N/nutrients
  tab_style(
    style = cell_borders(sides = "bottom", weight = px(2), color = "black"),
    locations = cells_body(rows = variable == "LOI (%)")
  ) %>%
  
  tab_style(
    style = cell_borders(sides = "bottom", weight = px(2), color = "black"),
    locations = cells_body(rows = variable == "VWC (cm³/cm³)")
  ) %>%
  # color column labels
  {
    tmp <- .
    for (fg in names(grp_cols)) {
      tmp <- tmp %>%
        tab_style(
          style = cell_fill(color = alpha(grp_cols[[fg]], 0.2)),
          locations = cells_column_labels(columns = fg)
        )
    }
    tmp
  } 
# %>%
#   tab_source_note(source_note = "* Bold indicates nutriens used in models; dark grey indicates Dark Model-only additions")


  

gt_table
# gtsave(gt_table, filename = "flux_summary_table_p_val_all.png")
# gtsave(gt_table, filename = "flux_summary_table.docx")

