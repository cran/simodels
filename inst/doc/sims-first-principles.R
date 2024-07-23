## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE------------------------------------------------------------
library(tmap)
library(dplyr)
library(ggplot2)

## ----inputs-------------------------------------------------------------------
zones = simodels::si_zones
centroids = simodels::si_centroids
od = simodels::si_od_census
tm_shape(zones) + tm_polygons("all", palette = "viridis")

## -----------------------------------------------------------------------------
od_df = od::points_to_od(centroids)
od_sfc = od::odc_to_sfc(od_df[3:6])
sf::st_crs(od_sfc) = 4326
od_df$length = sf::st_length(od_sfc)
od_df = od_df %>% transmute(
  O, D, length = as.numeric(length) / 1000,
  flow = NA, fc = NA
  )
od_df = sf::st_sf(od_df, geometry = od_sfc, crs = 4326)

## ----unconstrained------------------------------------------------------------
beta = 0.3
i = 1
j = 2
for(i in seq(nrow(zones))) {
  for(j in seq(nrow(zones))) {
    O = zones$all[i]
    n = zones$all[j]
    ij = which(od_df$O == zones$geo_code[i] & od_df$D == zones$geo_code[j])
    od_df$fc[ij] = exp(-beta * od_df$length[ij])
    od_df$flow[ij] = O * n * od_df$fc[ij]
  }
}
od_top = od_df %>% 
  filter(O != D) %>% 
  top_n(n = 2000, wt = flow)

tm_shape(zones) +
  tm_borders() +
  tm_shape(od_top) +
  tm_lines("flow")

## ----distance_decay-----------------------------------------------------------
summary(od_df$fc)
od_df %>% 
  ggplot() +
  geom_point(aes(length, fc))

## ----constrained--------------------------------------------------------------
od_dfj = left_join(
  od_df,
  zones %>% select(O = geo_code, all) %>% sf::st_drop_geometry()
)
od_dfj = od_dfj %>% 
  group_by(O) %>% 
  mutate(flow_constrained = flow / sum(flow) * first(all)) %>%
  ungroup()
sum(od_dfj$flow_constrained) == sum(zones$all)
od_top = od_dfj %>% 
  filter(O != D) %>% 
  top_n(n = 2000, wt = flow_constrained)

tm_shape(zones) +
  tm_borders() +
  tm_shape(od_top) +
  tm_lines("flow_constrained")

## ----validation---------------------------------------------------------------
od_dfjc = inner_join(od_dfj %>% select(-all), od)
od_dfjc %>% 
  ggplot() +
  geom_point(aes(all, flow_constrained))
cor(od_dfjc$all, od_dfjc$flow_constrained)^2

