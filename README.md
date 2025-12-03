# wstelemetrythreats-ms
Derived data and analysis code supporting the manuscript "Habitat use differences mediate population threat exposure in white sturgeon" by Walter and colleagues.

Raw telemetry data are not made available here due to file size constraints. These detection data have been publicly archived through the PATH database (https://fishdb.wfcb.ucdavis.edu), but were reformatted to match Ocean Tracking Network (OTN) standards in the process. Exact versions of raw detection data used in these analyses are available by request to lead author Jon Walter (jawalter@ucdavis.edu). In the interest of openness, code files include data munging steps that will not be fully reproducible without obtaining exact data versions, but can be sidestepped by using the provided derived (summarized) data products.

Code files included are:
1) ws_occupancy_analyses_wateryear_allfish.R: analysis script testing for differences in habitat occupancy; how they depend on lifestage, etc.; and how they translate to differences in exposure to angling.
2) ws_vs_hab.R: analysis script investigating timing of occupancy of areas where a severe HAB occurred in late Summer 2022.

Data files included are:
1) occupancyMatrix.csv: matrix of total occupancy (days) of unique fish-water year combinations at 15 sites. This is the main derived data product used in analyses.
2) fish_to_use.csv: vector of fish meeting selection criteria.
3) fishSummary.csv: miscellaneous attributes of fish selected for analysis.
4) generalArea_coords.csv: geographic coordinates of locations with telemetry receives.
5) sites_used.csv: geographic coordinates of locations used for analyses.
6) statesp020.shp: used in plotting
7) background_layer.shp: used in plotting
