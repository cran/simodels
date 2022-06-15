## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(simodels)
library(dplyr)
library(ggplot2)

## ----import-------------------------------------------------------------------
u1 = "https://github.com/Robinlovelace/simodels/releases/download/0.0.1/zones_aus.geojson"
zones_aus = sf::read_sf(u1)
u2 = "https://www.dropbox.com/s/wi3zxlq5pff1yda/AusMig2011.csv?raw=1"
od_aus = read.csv(u2)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  u1 = "https://www.dropbox.com/s/0fg80nzcxcsybii/GCCSA_2016_AUST_New.geojson?raw=1" # original
#  zones_original = sf::read_sf(u1)
#  sf::write_sf(zones_original, "GCCSA_2016_AUST_New.geojson")
#  lobstr::obj_size(zones_original$geometry) / lobstr::obj_size(zones_original)
#  # 99% of size taken up with geometry, no need for spatial resolution -> simplify
#  piggyback::pb_upload("GCCSA_2016_AUST_New.geojson")
#  readr::write_csv(od_aus, "od_migration_original.csv")
#  piggyback::pb_upload("od_migration_original.csv")
#  library(sf)
#  zones_aus_simplified = rmapshaper::ms_simplify(zones_original, keep = 0.001)
#  lobstr::obj_size(zones_aus_simplified) / lobstr::obj_size(zones_original)
#  # 0.2% of the size
#  plot(zones_aus_simplified$geometry[1]) # still too detailed?
#  sf::write_sf(zones_aus_simplified, "zones_aus.geojson")
#  piggyback::pb_upload("zones_aus.geojson")

## ----clean--------------------------------------------------------------------
dim(zones_aus)
names(zones_aus)
key_zone_names = c("GCCSA_CODE", "GCCSA_NAME", "AREA_SQKM")
zones = zones_aus[key_zone_names]
head(zones, 2)
dim(od_aus)
names(od_aus)
key_od_names = c("Orig_code", "Dest_code", "Flow")
od = od_aus[key_od_names]
head(od, 2)

## -----------------------------------------------------------------------------
summary(od[[1]] %in% zones[[1]])
summary(od[[2]] %in% zones[[1]])

## -----------------------------------------------------------------------------
od_sim = si_to_od(origins = zones, destinations = zones)
names(od_sim)

## ----unconstrained1-----------------------------------------------------------
si_power = function(d, beta) (d / 1000)^beta
od_calculated = si_calculate(
  od_sim,
  fun = si_power,
  d = distance_euclidean,
  beta = -0.8
  )
plot(od_calculated["interaction"], logz = TRUE)

## -----------------------------------------------------------------------------
od_interzonal = od %>%
  filter(Orig_code != Dest_code)
od_calculated_interzonal = od_calculated %>%
  filter(O != D) 
scale_factor = sum(od_interzonal$Flow) /
  sum(od_calculated_interzonal$interaction)
od_calculated_interzonal = od_calculated_interzonal %>% 
  mutate(interaction_scaled = interaction * scale_factor)
od_joined = inner_join(
  od_calculated_interzonal,
  od %>% rename(O = Orig_code, D = Dest_code)
  )
od_joined %>% 
  ggplot() +
  geom_point(aes(Flow, interaction_scaled))
cor(od_joined$Flow, od_joined$interaction_scaled)^2

## -----------------------------------------------------------------------------
od_joined %>% 
  mutate(decay = distance_euclidean^-0.8) %>% 
  mutate(decay = decay * (sum(Flow) / sum(decay))) %>% 
  ggplot() +
  geom_point(aes(distance_euclidean, Flow)) +
  geom_line(aes(distance_euclidean, decay), colour = "red") 

## -----------------------------------------------------------------------------
od_originating = od_joined %>% 
  group_by(O) %>% 
  mutate(originating_per_zone = sum(Flow)) %>% 
  ungroup()

## -----------------------------------------------------------------------------
od_constrained_p = si_calculate(
  od_originating,
  fun = si_power,
  d = distance_euclidean,
  beta = -0.8,
  constraint_production = originating_per_zone
  )
od_constrained_p %>% 
  ggplot() +
  geom_point(aes(Flow, interaction))
cor(od_constrained_p$Flow, od_constrained_p$interaction)^2

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  f = Flow ~ a * (distance_euclidean)^b
#  m = nls(
#    formula = f,
#    data = od_originating,
#    start = list(b = 0.8, a = 0.001),
#    upper = c(5, 1e5), lower = c(-5, 0.00001),
#    algorithm = "port"
#    )
#  m

## -----------------------------------------------------------------------------
library(minpack.lm)
f = Flow ~ a * (distance_euclidean)^b
m = nlsLM(
  formula = f,
  data = od_originating,
  )
m
# Nonlinear regression model
#   model: Flow ~ a * (distance_euclidean)^b
#    data: od_originating
#          a          b 
#  2.182e+07 -5.801e-01 

## -----------------------------------------------------------------------------
od_joined %>% 
  mutate(decay = distance_euclidean^-5.801e-01) %>% 
  mutate(decay = decay * 2.182e+07) %>% 
  ggplot() +
  geom_point(aes(distance_euclidean, Flow)) +
  geom_line(aes(distance_euclidean, decay), colour = "red") 

## -----------------------------------------------------------------------------
od_pred = si_predict(od_originating, model = m)
cor(od_pred$Flow, od_pred$interaction)^2
od_pred_const = si_predict(od_originating, model = m,
  constraint_production = originating_per_zone)
cor(od_pred_const$Flow, od_pred_const$interaction)^2

## -----------------------------------------------------------------------------
library(tmap)
ttm()
tm_shape(od_pred_const) +
  tm_lines("interaction_scaled", palette = "viridis")

