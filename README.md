gluc
====

Generalized Land Use Change (gluc) is a probabilistic cellular automata
program used by the Landuse Evolution and Assessment Model (LEAM). Briefly,
gluc reads a configuration file that provides the names and locations 
of population and employment growth projections, an initial land use
raster map, residential and commercial probabilities raster maps, an
optional no growth raster map, and optional residential and commercial
density maps.

Based on the desired growth projecttions and map inputs, gluc will change
the state of the land use map to meet the growth demand. After execution
the model will produce a land use change raster map (only modified
cells will have have values), a people per cell raster map (total 
number of people assigned to a cell), a employment per cell raster map (
total number of jobs assigned to a cell), and optionally a final
probility map.

Note
====
The gluc program may be used independently of LEAM but certain historic
anomalies may be confusing.  Hopefully these will be removed over time.
