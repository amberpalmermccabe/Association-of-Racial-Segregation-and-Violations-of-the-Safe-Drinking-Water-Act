library(sf)               
library(spdep)             
library(MASS)
library(pscl)
library(spatialreg)          
library(tidyverse)
library(writexl)
library(knitr)
library(kableExtra)
library(broom)
library(modelsummary)
library(viridis)
library(gtsummary)
###########################################################################
### HOLC Analysis

######## Descriptive Stats ######
# Calculate key statistics
violation_stats <- data.frame(
  Statistic = c("Mean", "Variance", "Variance/Mean Ratio", 
                "% Zeros", "Min", "Max", "Median"),
  Value = c(
    round(mean(HOLC_model_dat$violations), 3),
    round(var(HOLC_model_dat$violations), 3),
    round(var(HOLC_model_dat$violations) / mean(HOLC_model_dat$violations), 3),
    round(mean(HOLC_model_dat$violations == 0) * 100, 1),
    min(HOLC_model_dat$violations),
    max(HOLC_model_dat$violations),
    median(HOLC_model_dat$violations)
  )
)

print(violation_stats)
# mean 0.980, variance 22.47, mean/variance ratio 22.935
# 74.8% zeros 


##############################################################
# Creating Spatial Weights Matrix

# Get unique water systems
unique_holc_systems <- HOLC_model_dat %>%
  group_by(PWSID) %>%
  slice(1) %>%
  ungroup()

cat("Unique water systems:", nrow(unique_holc_systems), "\n") #642
cat("Total PWSID-year combinations:", nrow(HOLC_model_dat), "\n") #7062

# Create spatial weights on UNIQUE systems
centroids <- st_centroid(unique_holc_systems)
coords <- st_coordinates(centroids)

# K-nearest neighbors
nb_knn <- knn2nb(knearneigh(coords, k = 8))

# Contiguity (if feasible)
system.time({
  nb_queen <- poly2nb(unique_holc_systems, queen = TRUE)
})

# Hybrid approach
nb_hybrid_holc <- union.nb(nb_queen, nb_knn)

cat("\nHybrid neighbors summary:\n")
summary(nb_hybrid_holc)


# Convert to weights
listw_holc <- nb2listw(nb_hybrid_holc, style = "W", zero.policy = TRUE)

# Expand to full panel data
pwsid_lookup <- data.frame(
  PWSID = unique_holc_systems$PWSID,
  index = 1:nrow(unique_holc_systems)
)

HOLC_model_dat <- HOLC_model_dat %>%
  left_join(pwsid_lookup, by = "PWSID")

# Create expanded neighbor list for full panel
nb_panel_holc <- lapply(HOLC_model_dat$index, function(i) {
  nb_hybrid_holc[[i]]
})
class(nb_panel_holc) <- "nb"
attr(nb_panel_holc, "region.id") <- as.character(1:nrow(HOLC_model_dat))

listw_panel_holc <- nb2listw(nb_panel_holc, style = "W", zero.policy = TRUE)

# Test for spatial autocorrelation
moran_holc <- moran.test(HOLC_model_dat$violations, listw_panel_holc)
print(moran_holc)

# Create spatial lag
HOLC_model_dat <- HOLC_model_dat %>%
  mutate(violations_lag = lag.listw(listw_panel_holc, violations))

########################################
##### Run Models
library(MASS)
# Model: Zero-Inflated Poisson
#zip_holc <- zeroinfl(
# violations ~ wg_HOLC + ownrshp + pop_dns_1000 + 
#   ws_size + GW_SW_C + year_c |
#  wg_HOLC + ownrshp + ws_size,
#data = HOLC_model_dat,
#dist = "poisson"
#)

# Model 1: Negative Binomial (non-spatial)
nb_holc <- glm.nb(
  violations ~ wg_HOLC + ownrshp + pop_dns_1000  + GW_SW_C + ws_size ,
  data = HOLC_model_dat
)

#nb_holc_offset <- glm.nb(
#  violations ~ wg_HOLC + ownrshp + pop_dns_1000 + GW_SW_C +
#  as.factor(year_c) +  offset(log(POP_SERVED)),
# data = HOLC_model_dat
#)

# Model 2:Negative Binomial (with lag)
spatial_nb_holc <- glm.nb(
  violations ~ wg_HOLC + ownrshp + pop_dns_1000 + 
    ws_size + GW_SW_C  + violations_lag,
  data = HOLC_model_dat
)

#spatial_nb_holc_offset <- glm.nb(
# violations ~ wg_HOLC + ownrshp + pop_dns_1000 + 
#  ws_size + GW_SW_C + year_c + violations_lag + offset(log(POP_SERVED)),
#data = HOLC_model_dat
#)

# Check summaries
#summary(zip_holc)
summary(nb_holc)
#summary(nb_holc_offset)
summary(spatial_nb_holc)
#summary((spatial_nb_holc_offset))

#############################
## Stratified Analysis
strata <- levels(HOLC_model_dat$ws_size)

stratified_models <- lapply(strata, function(s) {
  dat <- HOLC_model_dat[HOLC_model_dat$ws_size == s, ]
  glm.nb(
    violations ~ wg_HOLC + ownrshp + pop_dns_1000 + GW_SW_C +
      as.factor(year),
    data = dat          # <-- comma after formula, and use dat not HOLC_model_dat
  )
})

names(stratified_models) <- strata

library(broom)

results <- lapply(strata, function(s) {
  tidy(stratified_models[[s]], exponentiate = TRUE, conf.int = TRUE) |>
    mutate(stratum = s)
}) |> bind_rows()

# Format the table
table_out <- results |>
  mutate(
    IRR     = round(estimate, 3),
    CI      = paste0("(", round(conf.low, 3), ", ", round(conf.high, 3), ")"),
    p.value = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(p.value, 3))
    )
  ) |>
  select(term, stratum, IRR, CI, p.value) |>
  tidyr::pivot_wider(
    names_from  = stratum,
    values_from = c(IRR, CI, p.value)
  )

write_xlsx(table_out, "HOLC_stratified_FINAL.xlsx")


########################
# Model Comparison

# Extract key statistics
comparison_table <- data.frame(
  Model = c("Negative Binomial",
            "Negative Binomial with Offset", 
            "Negative Binomial with Spatial Lag",
            "Negative Binomial w/ Lag and Offset"),
  
  AIC = c(AIC(nb_holc), 
          AIC(nb_holc_offset),
          AIC(spatial_nb_holc),
          AIC(spatial_nb_holc_offset)),
  
  BIC = c(BIC(nb_holc), 
          BIC(nb_holc_offset),
          BIC(spatial_nb_holc),
          BIC(spatial_nb_holc_offset)),
  
  LogLik = c(logLik(nb_holc)[1], 
             logLik(nb_holc_offset)[1],
             logLik(spatial_nb_holc)[1],
             logLik(spatial_nb_holc_offset)),
  
  N_Parameters = c(attr(logLik(nb_holc), "df"),
                   attr(logLik(nb_holc_offset), "df"),
                   attr(logLik(spatial_nb_holc), "df"),
                   attr(logLik(spatial_nb_holc_offset), "df")),
  
  N_Observations = c(length(nb_holc$y),
                     nobs(nb_holc_offset),
                     nobs(spatial_nb_holc),
                     nobs(spatial_nb_holc_offset))
)

# Calculate delta AIC (relative to best model)
comparison_table <- comparison_table %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Delta_BIC = BIC - min(BIC),
    Best_AIC = ifelse(Delta_AIC == 0, "✓", ""),
    Best_BIC = ifelse(Delta_BIC == 0, "✓", "")
  )

# Round for readability
comparison_table <- comparison_table %>%
  mutate(
    AIC = round(AIC, 1),
    BIC = round(BIC, 1),
    LogLik = round(LogLik, 1),
    Delta_AIC = round(Delta_AIC, 1),
    Delta_BIC = round(Delta_BIC, 1)
  )

print(comparison_table)

####### Violations by type
viol_type <- model_data %>% 
  group_by(PWSID) %>%
  summarise(
    violations       = n_distinct(VIOLATI),
    n_viol_TT        = n_distinct(VIOLATI[viol_ct == "TT"]),
    n_viol_MRDL      = n_distinct(VIOLATI[viol_ct == "MRDL"]),
    n_viol_MR        = n_distinct(VIOLATI[viol_ct == "MR"]),
    n_viol_OTHER     = n_distinct(VIOLATI[viol_ct == "OTHER"]),
    n_viol_MCL       = n_distinct(VIOLATI[viol_ct == "MCL"]),
    n_viol_MON       = n_distinct(VIOLATI[viol_ct == "MON"]),
    n_viol_RPT       = n_distinct(VIOLATI[viol_ct == "RPT"]),
    wg_HOLC          = first(wg_HOLC),
    ownrshp          = first(ownrshp),
    ws_size          = first(ws_size),
    GW_SW_C          = first(GW_SW_C),
    POP_SERVED       = first(POPULAT),
    pop_dns          = first(pop_dns),
    .groups = "drop"
  )

viol_type <- viol_type %>% 
  mutate(ownership_recoded = case_when( 
    ownrshp %in% c("F", "L", "S") ~ "public", 
    ownrshp %in% c("P", "M") ~ "private", 
    TRUE ~ NA_character_)) # in case there are other values ))

viol_model_dat <- viol_type %>%
  select(-ownrshp)

viol_model_dat <- viol_model_dat %>% 
  rename(ownrshp = ownership_recoded)

# Set reference category for a single factor
viol_model_dat <- viol_model_dat %>%
  mutate(ownrshp = relevel(factor(ownrshp), ref = "private"))

# For water system size
viol_model_dat <- viol_model_dat %>%
  mutate(ws_size = relevel(factor(ws_size), ref = "Large"))

# For groundwater/surface water
viol_model_dat <- viol_model_dat %>%
  mutate(GW_SW_C = relevel(factor(GW_SW_C), ref = "GW"))


# Total violations 
nb_total <- MASS::glm.nb(
  violations ~ wg_HOLC + pop_dns + GW_SW_C + ws_size + ownrshp,
  data = viol_model_dat
)
summary(nb_total)

# Incidence Rate Ratios
irr_total <- exp(cbind(IRR = coef(nb_total), confint(nb_total)))
print(round(irr_total, 3))


#  violation type 
# Health-based violations (regulatory failures - most substantively important)
nb_MCL <- MASS::glm.nb(
  n_viol_MCL ~ wg_HOLC + pop_dns + GW_SW_C + ws_size + ownrshp,
  data = viol_model_dat
)

nb_TT <- MASS::glm.nb(
  n_viol_TT ~ wg_HOLC + pop_dns + GW_SW_C + ws_size + ownrshp,
  data = viol_model_dat
)

# No water systems with this violation in this subset
#nb_MRDL <- MASS::glm.nb(
# n_viol_MRDL ~ wg_HOLC + pop_dns + GW_SW_C + ws_size + ownrshp,
#data = viol_model_dat
#)

# Monitoring/reporting violations (administrative capacity)
nb_MON <- MASS::glm.nb(
  n_viol_MON ~ wg_HOLC + pop_dns + GW_SW_C + ws_size + ownrshp,
  data = viol_model_dat
)

nb_MR <- MASS::glm.nb(
  n_viol_MR ~ wg_HOLC + pop_dns + GW_SW_C + ws_size + ownrshp,
  data = viol_model_dat
)

nb_RPT <- MASS::glm.nb(
  n_viol_RPT ~ wg_HOLC + pop_dns + GW_SW_C + ws_size + ownrshp,
  data = viol_model_dat
)

# --- Compile IRRs for all models ---
models <- list(
  Total  = nb_total,
  MCL    = nb_MCL,
  TT     = nb_TT,
  MON    = nb_MON,
  MR     = nb_MR,
  RPT    = nb_RPT
)

irr_table <- lapply(names(models), function(m) {
  irr <- exp(cbind(IRR = coef(models[[m]]), confint(models[[m]])))
  as.data.frame(round(irr, 3)) %>%
    rownames_to_column("term") %>%
    mutate(outcome = m)
}) %>%
  bind_rows() %>%
  filter(term == "wg_HOLC") %>%   # pull just the exposure coefficient
  dplyr::select(outcome, IRR, `2.5 %`, `97.5 %`)

print(irr_table)

############################################################################
### ICE Analysis

########################################
# Checking mean and variance to determine if poisson or ZIP is more appropriate
# Basic descriptives
mean(model_data_recode$violations) #  1.6625
var(model_data_recode$violations) #  79.6586

# Dispersion ratio
var(model_data_recode$violations) / mean(model_data_recode$violations)  #   47.914
# If ≈ 1: equidispersion (Poisson ok)
# If > 1: overdispersion (consider Negative Binomial)
# If < 1: underdispersion (rare)

# Check zero inflation
observed_zeros <- mean(model_data_recode$violations == 0) #  Observed proportion of zeros: 0.6916
expected_zeros <- exp(-mean(model_data_recode$violations))  #  Expected zeros under Poisson: 0.1896

cat("Observed proportion of zeros:", observed_zeros, "\n")
cat("Expected zeros under Poisson:", expected_zeros, "\n")
cat("Ratio (observed/expected):", observed_zeros / expected_zeros, "\n")
# If ratio > 1.5: zero-inflation present
## Ratio (observed/expected): 3.646


#####
######################################
# Create quintiles of ICErace
model_data_recode <- model_data_recode %>%
  mutate(
    ICErace_quint = cut(ICErace, 
                        breaks = quantile(ICErace, probs = seq(0, 1, 0.2), na.rm = TRUE),
                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                        include.lowest = TRUE),
    ICErace_quint = relevel(factor(ICErace_quint), ref = "Q5")  
  )

# Check distribution
table(model_data_recode$ICErace_quint)

# Negative Binomial
nb_model <- glm.nb(
  violations ~ ICErace_quint + ownership_recoded + pop_dns_1000 + ws_size + GW_SW_CODE + as.factor(year),
  data = model_data_recode
)
nb_model_offset <- glm.nb(
  violations ~ ICErace_quint + ownership_recoded + pop_dns_1000 + ws_size + GW_SW_CODE + as.factor(year) + offset(log(POP_SERVED)),
  data = model_data_recode
)

########## Spatial Lag ##################

########## Create spatial weights

# Get centroids
centroids <- st_centroid(model_data_recode)
coords <- st_coordinates(centroids)

# Contiguity-based neighbors (touching polygons)
nb_queen <- poly2nb(model_data_recode, queen = TRUE)
cat("Contiguity neighbors:\n")
summary(nb_queen)

# K-nearest neighbors (backup for isolated polygons)
nb_knn <- knn2nb(knearneigh(coords, k = 5))
cat("\nK-nearest neighbors (k=5):\n")
summary(nb_knn)

# Combine: union ensures no isolated observations
nb_hybrid <- union.nb(nb_queen, nb_knn)
cat("\nHybrid neighbors (contiguity + KNN):\n")
summary(nb_hybrid)
table(card(nb_hybrid))  # Check all have neighbors

# Convert to spatial weights list
listw_hybrid <- nb2listw(nb_hybrid, style = "W", zero.policy = TRUE)

model_data_recode <- model_data_recode %>%
  mutate(violations_lag = lag.listw(listw_hybrid, violations))


################ Test for spatial autocorrelation

moran_ice_before <- moran.test(model_data_recode$violations, listw_hybrid)
print(moran_ice_before)

# Global Moran's I plot
moran.plot(model_data_recode$violations, listw_hybrid,
           xlab = "Violations", ylab = "Spatially Lagged Violations",
           main = "Moran's I Scatterplot")


#######################################################
##### Main model (negative binomial with spatial lag)

# See the actual ICE ranges
model_data_recode %>%
  group_by(ICErace_quint) %>%
  summarise(
    n = n(),
    min_ICE = min(ICErace, na.rm = TRUE),
    mean_ICE = mean(ICErace, na.rm = TRUE),
    max_ICE = max(ICErace, na.rm = TRUE)
  )

# Main model with quintiles
main_nonspatial_quint <- glm.nb(
  violations ~ ICErace_quint + ownership_recoded + pop_dns_1000 + 
    ws_size + GW_SW_CODE + as.factor(year),
  data = model_data_recode
)

main_spatial_quint <- glm.nb(
  violations ~ ICErace_quint + ownership_recoded + pop_dns_1000 + 
    ws_size + GW_SW_CODE + as.factor(year) + violations_lag,
  data = model_data_recode
)

main_spatial_quint_offset <- glm.nb(
  violations ~ ICErace_quint + ownership_recoded + pop_dns_1000 + 
    ws_size + GW_SW_CODE + as.factor(year) + violations_lag + offset(log(POP_SERVED)),
  data = model_data_recode
)

summary(main_spatial_quint)
summary(main_spatial_quint_offset)

# Interaction model with quintiles
interact_spatial_quint <- glm.nb(
  violations ~ ICErace_quint * as.factor(year) + ownership_recoded + pop_dns_1000 + 
    ws_size + GW_SW_CODE + violations_lag,
  data = model_data_recode
)

interact_spatial_quint_offset <- glm.nb(
  violations ~ ICErace_quint * as.factor(year) + ownership_recoded + pop_dns_1000 + 
    ws_size + GW_SW_CODE + violations_lag + offset(log(POP_SERVED)),
  data = model_data_recode
)
summary(interact_spatial_quint)
summary(interact_spatial_quint_offset)

# Compare
AIC(main_nonspatial_quint, main_spatial_quint, interact_spatial_quint)
AIC(main_spatial_quint_offset, interact_spatial_quint_offset)

# Testing an interaction between ICE and ws_size
interact_size_ICE <- glm.nb(
  violations ~ ICErace_quint * ws_size + ownership_recoded + pop_dns_1000  + GW_SW_CODE 
  + violations_lag + as.factor(year),
  data = model_data_recode
)

summary (interact_size_ICE)
############## test residuals for autocorrelation
# Get row indices actually used in the model
model_indices <- as.numeric(names(residuals(main_spatial_quint)))

# Subset to only rows used in model
model_data_used <- model_data_recode[model_indices, ]

cat("\n========== SPATIAL WEIGHTS CHECK ==========\n")
cat("Original data rows:", nrow(model_data_recode), "\n")
cat("Model used rows:", nrow(model_data_used), "\n")
cat("Residuals length:", length(residuals(main_spatial_quint)), "\n")
cat("Match:", nrow(model_data_used) == length(residuals(main_spatial_quint)), "\n")

# Create spatial weights for ONLY the observations used in model
coords_used <- st_coordinates(st_centroid(model_data_used))
nb_used <- knn2nb(knearneigh(coords_used, k = 8))
listw_used <- nb2listw(nb_used, style = "W")

# Verify lengths match
cat("\nSpatial weights length:", length(listw_used$neighbours), "\n")

# Test AFTER model (on residuals)
cat("\n========== SPATIAL AUTOCORRELATION TEST ==========\n")
moran_ice_after <- moran.test(residuals(main_spatial_quint), listw_used)
print(moran_ice_after)

# Summary
cat("\n========== SUMMARY ==========\n")
cat("BEFORE model - Moran's I:", round(moran_ice_before$estimate[1], 4), 
    "p-value:", format.pval(moran_ice_before$p.value), "\n")
cat("AFTER model - Moran's I:", round(moran_ice_after$estimate[1], 4), 
    "p-value:", format.pval(moran_ice_after$p.value), "\n")

### Moran's I before: 0.1056 p < 2.22e-16, after 0.0325 p < 2.22e-16

################################
######## Stratified Analysis###

strata <- unique(model_data_recode$ws_size)

ice_stratified_no_offset <- lapply(strata, function(s) {
  dat <- model_data_recode[model_data_recode$ws_size == s, ]
  glm.nb(
    violations ~ ICErace_quint + ownership_recoded + pop_dns_1000 + GW_SW_CODE + as.factor(year),
    data = dat)
})
names(ice_stratified_no_offset) <- strata

####### ICE ethnicity 
# Create quintiles of ICEeth
model_data_recode <- model_data_recode %>%
  mutate(
    ICEeth_quint = cut(ICEeth, 
                       breaks = quantile(ICEeth, probs = seq(0, 1, 0.2), na.rm = TRUE),
                       labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                       include.lowest = TRUE),
    ICEeth_quint = relevel(factor(ICEeth_quint), ref = "Q5")  
  )

# Check distribution
table(model_data_recode$ICEeth_quint)


########################
# Negative Binomial
nb_model_eth <- glm.nb(
  violations ~ ICEeth_quint + ownership_recoded + pop_dns_1000 + ws_size + GW_SW_CODE + as.factor(year),
  data = model_data_recode
)

# NB with spatial lag
nb_lag_eth <- glm.nb(
  violations ~ ICEeth_quint + ownership_recoded + pop_dns_1000 + 
    ws_size + GW_SW_CODE + as.factor(year) + violations_lag,
  data = model_data_recode
)


# Interaction model with quintiles
nb_lag_interaction <- glm.nb(
  violations ~ ICEeth_quint * as.factor(year) + ownership_recoded + pop_dns_1000 + 
    ws_size + GW_SW_CODE + violations_lag,
  data = model_data_recode
)


# Compare
AIC(nb_lag_eth,  nb_lag_interaction)
BIC(nb_lag_eth, nb_lag_interaction)

####################################################
# Sensitivity Analysis Restricting ICE to HOLC area

# Filter model_data_recode to only include PWSIDs that are in HOLC_model_dat
combined_data <- model_data_recode %>%
  filter(PWSID %in% HOLC_model_dat$PWSID)

########## Create spatial weights

# Get centroids
centroids <- st_centroid(combined_data)
coords <- st_coordinates(centroids)

# Contiguity-based neighbors (touching polygons)
nb_queen <- poly2nb(combined_data, queen = TRUE)
cat("Contiguity neighbors:\n")
summary(nb_queen)

# K-nearest neighbors (backup for isolated polygons)
nb_knn <- knn2nb(knearneigh(coords, k = 5))
cat("\nK-nearest neighbors (k=5):\n")
summary(nb_knn)

# Combine: union ensures no isolated observations
nb_hybrid <- union.nb(nb_queen, nb_knn)
cat("\nHybrid neighbors (contiguity + KNN):\n")
summary(nb_hybrid)
table(card(nb_hybrid))  # Check all have neighbors

# Convert to spatial weights list
listw_hybrid <- nb2listw(nb_hybrid, style = "W", zero.policy = TRUE)

combined_data <- combined_data %>%
  mutate(violations_lag = lag.listw(listw_hybrid, violations))

# Negative Binomial
nb_sens <- glm.nb(
  violations ~ ICErace_quint + ownership_recoded + pop_dns_1000 + ws_size + GW_SW_CODE + as.factor(year),
  data = combined_data
)

# NB with spatial lag
nb_lag_sens <- glm.nb(
  violations ~ ICErace_quint + ownership_recoded + pop_dns_1000 + 
    ws_size + GW_SW_CODE + as.factor(year) + violations_lag,
  data = combined_data
)

##############################################################################
### Violations by type

# Fit nb models — ICErace

models_ICE_race <- list(
  Total = MASS::glm.nb(violations  ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  MCL   = MASS::glm.nb(n_viol_MCL  ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  TT    = MASS::glm.nb(n_viol_TT   ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  MON   = MASS::glm.nb(n_viol_MON  ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  MR    = MASS::glm.nb(n_viol_MR   ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  RPT   = MASS::glm.nb(n_viol_RPT  ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  OTHER = MASS::glm.nb(n_viol_OTHER ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  MRDL  = MASS::glm.nb(n_viol_MRDL ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean)
)

# Fit Spatial lag nb models — ICErace
models_ICE_race_slag <- list(
  Total = MASS::glm.nb(violations   ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_violations, data = dat_clean),
  MCL   = MASS::glm.nb(n_viol_MCL   ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_MCL,        data = dat_clean),
  TT    = MASS::glm.nb(n_viol_TT    ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_TT,         data = dat_clean),
  MON   = MASS::glm.nb(n_viol_MON   ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_MON,        data = dat_clean),
  MR    = MASS::glm.nb(n_viol_MR    ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_MR,         data = dat_clean),
  RPT   = MASS::glm.nb(n_viol_RPT   ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_RPT,        data = dat_clean),
  OTHER = MASS::glm.nb(n_viol_OTHER ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_OTHER,      data = dat_clean),
  MRDL  = MASS::glm.nb(n_viol_MRDL  ~ ICErace + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_MRDL,       data = dat_clean)
)

# Fit nb models — ICEeth
models_ICE_eth <- list(
  Total = MASS::glm.nb(violations   ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  MCL   = MASS::glm.nb(n_viol_MCL   ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  TT    = MASS::glm.nb(n_viol_TT    ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  MON   = MASS::glm.nb(n_viol_MON   ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  MR    = MASS::glm.nb(n_viol_MR    ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  RPT   = MASS::glm.nb(n_viol_RPT   ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  OTHER = MASS::glm.nb(n_viol_OTHER ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean),
  MRDL  = MASS::glm.nb(n_viol_MRDL  ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership, data = dat_clean)
)

# Fit Spatial lag nb models- ICEeth
models_ICE_eth_slag <- list(
  Total = MASS::glm.nb(violations   ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_violations, data = dat_clean),
  MCL   = MASS::glm.nb(n_viol_MCL   ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_MCL,        data = dat_clean),
  TT    = MASS::glm.nb(n_viol_TT    ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_TT,         data = dat_clean),
  MON   = MASS::glm.nb(n_viol_MON   ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_MON,        data = dat_clean),
  MR    = MASS::glm.nb(n_viol_MR    ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_MR,         data = dat_clean),
  RPT   = MASS::glm.nb(n_viol_RPT   ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_RPT,        data = dat_clean),
  OTHER = MASS::glm.nb(n_viol_OTHER ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_OTHER,      data = dat_clean),
  MRDL  = MASS::glm.nb(n_viol_MRDL  ~ ICEeth + pop_dns + GW_SW_CODE + ws_size + ownership + log_slag_MRDL,       data = dat_clean)
)

# Check convergence for all model sets
for (model_set_name in c("models_ICE_race", "models_ICE_race_slag",
                         "models_ICE_eth",  "models_ICE_eth_slag")) {
  model_set <- get(model_set_name)
  results <- sapply(names(model_set), function(m) {
    tryCatch({ nobs(model_set[[m]]); "converged" }, error = function(e) "FAILED")
  })
  cat("\n", model_set_name, ":\n"); print(results)
}


extract_irr <- function(model, exposure) {
  irr   <- exp(cbind(IRR = coef(model), confint.default(model)))
  pvals <- summary(model)$coefficients[, "Pr(>|z|)"]
  as.data.frame(round(irr, 3)) %>%
    rownames_to_column("term") %>%
    mutate(p_value = round(pvals, 4)) %>%
    filter(term == exposure)
}

get_obs <- function(models, dataset) {
  data.frame(
    outcome = names(models),
    n_systems = sapply(names(models), function(m) nobs(models[[m]])),
    n_nonzero = sapply(names(models), function(m) {
      v <- switch(m,
                  Total = "violations", MCL = "n_viol_MCL", TT    = "n_viol_TT",
                  MON   = "n_viol_MON", MR  = "n_viol_MR",  RPT   = "n_viol_RPT",
                  OTHER = "n_viol_OTHER", MRDL = "n_viol_MRDL"
      )
      sum(dataset[[v]] > 0, na.rm = TRUE)
    }),
    total_violations = sapply(names(models), function(m) {
      v <- switch(m,
                  Total = "violations", MCL = "n_viol_MCL", TT    = "n_viol_TT",
                  MON   = "n_viol_MON", MR  = "n_viol_MR",  RPT   = "n_viol_RPT",
                  OTHER = "n_viol_OTHER", MRDL = "n_viol_MRDL"
      )
      sum(dataset[[v]], na.rm = TRUE)
    })
  )
}

build_comparison_table <- function(models_base, models_slag, dataset, exposure) {
  obs <- get_obs(models_base, dataset)
  
  nb_results <- lapply(names(models_base), function(m) {
    extract_irr(models_base[[m]], exposure) %>%
      mutate(outcome = m, model = "NB",
             AIC = round(AIC(models_base[[m]]), 1))
  }) %>% bind_rows()
  
  slag_results <- lapply(names(models_slag), function(m) {
    extract_irr(models_slag[[m]], exposure) %>%
      mutate(outcome = m, model = "Spatial Lag NB",
             AIC = round(AIC(models_slag[[m]]), 1))
  }) %>% bind_rows()
  
  bind_rows(nb_results, slag_results) %>%
    left_join(obs, by = "outcome") %>%
    mutate(CI = paste0("(", `2.5 %`, ", ", `97.5 %`, ")")) %>%
    dplyr::select(outcome, model, n_systems, n_nonzero,
                  total_violations, IRR, CI, p_value, AIC) %>%
    arrange(outcome, model)
}

full_coefs <- function(models, label) {
  lapply(names(models), function(m) {
    irr   <- exp(cbind(IRR = coef(models[[m]]), confint.default(models[[m]])))
    pvals <- summary(models[[m]])$coefficients[, "Pr(>|z|)"]
    as.data.frame(round(irr, 3)) %>%
      rownames_to_column("term") %>%
      mutate(p_value = round(pvals, 4), outcome = m, model = label)
  }) %>%
    bind_rows() %>%
    dplyr::select(model, outcome, term, IRR, `2.5 %`, `97.5 %`, p_value)
}

model_fit_table <- function(models, exposure, model_label) {
  lapply(names(models), function(m) {
    data.frame(
      exposure  = exposure,
      model     = model_label,
      outcome   = m,
      AIC       = round(AIC(models[[m]]), 1),
      theta     = round(models[[m]]$theta, 3),
      null_dev  = round(models[[m]]$null.deviance, 1),
      resid_dev = round(models[[m]]$deviance, 1),
      df        = models[[m]]$df.residual,
      n         = nobs(models[[m]])
    )
  }) %>% bind_rows()
}
