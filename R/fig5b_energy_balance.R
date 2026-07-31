# Figure 5B — leaf energy balance across stomatal conductance
# No input files; the model is parameterised in the constants below.

library(ggplot2)
library(rlang)

# --- physical constants ------------------------------------------------------
emissivity       <- 0.96
stefan_boltzmann <- 5.67e-8   # W m^-2 K^-4
lambda           <- 44000     # latent heat of vaporisation, J mol^-1
Cp               <- 29.3      # specific heat of air, J mol^-1 °C^-1

# --- conditions ---------------------------------------------------------------
# Chamber conditions at the 60 °C treatment, with leaf properties measured for
# T. oblongifolia.
T_AIR       <- 60      # air temperature, °C
RH          <- 0.3     # relative humidity, 0-1
P           <- 101     # atmospheric pressure, kPa
SUN_DIRECT  <- 275     # direct shortwave radiation, W m^-2
SUN_DIFFUSE <- 20      # diffuse shortwave radiation, W m^-2
SKY_TEMP    <- 55      # sky (chamber surface) temperature, °C
LEAF_ANGLE  <- 0       # from horizontal, degrees
LEAF_LENGTH <- 0.017   # characteristic length, m
WIND_SPEED  <- 0.1     # m s^-1
RHO         <- 0.65    # leaf shortwave reflectance

# observed range of minimum ΔT at 60 °C (Fig. 5A)
OBS_DELTA <- c(-13, -10)

# --- energy balance ----------------------------------------------------------
# boundary layer conductance to heat, mol m^-2 s^-1
calc_gH <- function(wind_speed, leaf_length)
  2 * 0.135 * sqrt(wind_speed / leaf_length)

# stomatal and boundary layer conductances act in series
harmonic_mean_conductance <- function(g_vapor, gH) 1 / ((1 / g_vapor) + (1 / gH))

# slope of the saturation vapour pressure curve, kPa °C^-1
slope_vapor_pressure <- function(T_air)
  (4098.17 / (T_air + 237.31)^2) * 0.611 * exp((17.27 * T_air) / (T_air + 237.31))

saturated_vapor_pressure <- function(T_air)
  0.611 * exp((17.27 * T_air) / (T_air + 237.31))

# absorbed shortwave plus incoming longwave, W m^-2
calc_total_radiation <- function(leaf_angle, sun_direct, sun_diffuse, sky_temp, rho) {
  longwave <- emissivity * stefan_boltzmann * (sky_temp + 273.15)^4
  absorbed <- (1 - rho) * (cos(leaf_angle * pi / 180) * sun_direct + sun_diffuse)
  absorbed + longwave
}

calc_rad_emission <- function(T_C)
  emissivity * stefan_boltzmann * (T_C + 273.15)^4

# leaf-air temperature difference from the energy balance, °C
calc_delta_T <- function(radiation, rad_emission, gvtotal, e_air_sat, RH, P, gH, s) {
  ((radiation - rad_emission) - (lambda * gvtotal * e_air_sat * (1 - RH) / P)) /
    ((gH * Cp) + (lambda * s * gvtotal / P))
}

# transpiration, mmol m^-2 s^-1
calc_transpiration_E <- function(gvtotal, P, e_air_sat, delta_T, s, RH)
  1000 * (gvtotal / P) * ((e_air_sat + delta_T * s) - e_air_sat * RH)

# --- sweep over stomatal conductance -----------------------------------------
# Radiative emission depends on leaf temperature, which is what we are solving
# for, so ΔT is computed twice: once assuming the leaf is at air temperature,
# then again using the emission implied by that first estimate.
conductance_sweep <- function(g_min = 0.01, g_max = 1.5, interval = 0.01) {
  gH        <- calc_gH(WIND_SPEED, LEAF_LENGTH)
  s         <- slope_vapor_pressure(T_AIR)
  e_air_sat <- saturated_vapor_pressure(T_AIR)
  radiation <- calc_total_radiation(LEAF_ANGLE, SUN_DIRECT, SUN_DIFFUSE,
                                    SKY_TEMP, RHO)
  
  do.call(rbind, lapply(seq(g_min, g_max, by = interval), function(g_vapor) {
    gvtotal <- harmonic_mean_conductance(g_vapor, gH)
    
    delta_T_0 <- calc_delta_T(radiation, calc_rad_emission(T_AIR),
                              gvtotal, e_air_sat, RH, P, gH, s)
    emission  <- calc_rad_emission(T_AIR + delta_T_0)
    delta_T   <- calc_delta_T(radiation, emission,
                              gvtotal, e_air_sat, RH, P, gH, s)
    
    data.frame(g_vapor = g_vapor,
               delta_T = delta_T,
               E       = calc_transpiration_E(gvtotal, P, e_air_sat,
                                              delta_T, s, RH))
  }))
}

sweep_df <- conductance_sweep()

# --- plot --------------------------------------------------------------------
p <- ggplot(sweep_df, aes(x = g_vapor)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = OBS_DELTA[1], ymax = OBS_DELTA[2],
           fill = "grey80", alpha = 0.4) +
  geom_line(aes(y = delta_T, color = "Delta T"), linewidth = 2) +
  geom_line(aes(y = E, color = "E"), linewidth = 2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text", x = Inf, y = OBS_DELTA[2],
           label = expression("Obs Min " * Delta * italic(T)),
           hjust = 1.1, vjust = -0.5, size = 4, family = "Helvetica",
           color = "grey40") +
  scale_color_manual(
    values = c("Delta T" = "#db7b9b", "E" = "#475643"),
    labels = c("Delta T" = expression(Delta * italic(T) ~ "(°C)"),
               "E" = expression(italic(E) ~ "(" * mmol ~ m^{-2} ~ s^{-1} * ")"))
  ) +
  scale_y_continuous(breaks = seq(-13, 11, by = 3),
                     labels = function(x) gsub("-", "\u2212", x)) +
  guides(color = guide_legend(keywidth = 0.8, keyheight = 0.8)) +
  labs(x = expression(italic(g)[sw] ~ "(mol m"^{-2} ~ s^{-1} * ")"),
       y = expression(Delta * italic(T) ~ "and" ~ italic(E))) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.title      = element_blank(),
    legend.text       = element_text(size = 11),
    legend.position   = c(0.70, 0.75),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key        = element_blank(),
    axis.text         = element_text(size = 12),
    axis.title.x      = element_text(size = 13),
    axis.title.y      = element_text(size = 14),
    axis.title.y.left = element_text(margin = margin(r = -5)),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA)
  )

pdf("fig5b_energy_balance.pdf", width = 3.6, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p); dev.off()