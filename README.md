# Svalbard Reindeer Movement Analysis

This project analyses GPS movement data from the **Svalbard Reindeer Project – Nordenskiöld Land** (Movebank Study ID: 2608802883) as part of a scientific visualisation assignment. The analysis focuses on daily movement distances and spatial behaviour of three individual reindeer (B244, B245, B251) tracked across the Svalbard archipelago.

## Visualisation

![Svalbard Reindeer Movement](svalbard_reindeer.png)

The figure shows three GPS-tagged individuals plotted on a satellite basemap, with fixes coloured by date (purple = earlier, yellow = more recent). Below each map, a histogram shows the distribution of daily distances travelled per individual. The colour scale and legendry-based histogram guide are shared across both panels, linking spatial location to movement intensity.

Key patterns visible in the data:
- **B244** shows the most concentrated home range with the majority of days spent travelling less than 2 km
- **B245** and **B251** show broader spatial spread and more variable daily movement, with occasional long-distance days exceeding 10 km
- All three individuals show a right-skewed distribution of daily distances, typical of central-place foraging behaviour

## Data

- **Source:** Svalbard Reindeer Project – Nordenskiöld Land, accessed via [Movebank](https://www.movebank.org)
- **Study ID:** 2608802883
- **Sensor type:** GPS
- **Individuals:** B244, B245, B251
- **Temporal resolution:** Even-hour fixes only

## Methods

- Data downloaded and processed using the `move2` R package
- Fixes filtered to even hours for data availability
- Step lengths calculated per fix in UTM coordinates (zone 33N)
- Daily total distances aggregated per individual
- Visualisation built with `ggplot2`, `tidyterra`, `maptiles`, `patchwork`, and `legendry`

## Requirements

```r
install.packages(c("move2", "sf", "terra", "ggplot2", "lubridate",
                   "dplyr", "legendry", "tidyterra", "maptiles",
                   "patchwork", "showtext", "rnaturalearth", "mapview"))
```

Access to Movebank requires a registered account and study download permission.

## Author

Marlene Sehrbrock — University of Würzburg, MSc EAGLE  
Created for Scientific Graphs Assignment, Semester 2
