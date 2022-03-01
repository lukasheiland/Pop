# README for data adaption

## Tables
Two tables are assumed: ``Stages``, ``Stages_transitions``. Unfortunately, the heuristics for transitions are very sampling design-specific, so that ``Stages_transitions`` has to be generated independently.

- ``Stages``
	- format: **sf** data.frame with POINT geometries
	- find column specifications in ``spec_Stages.csv``

- ``Stages_transitions``
	- format: data.frame
	- find column specifications in ``spec_Stages_transitions.csv``
	- This table is used to estimate the transition rates per year ``g`` (from J to A) and ``h`` (from A to B) on a plot. There is a quite complicated data prep involved to fit a simple model here. Instead of providing this table, you can also use other means to estimate prior distributions of the above transition rates and feed them in manually (just set target weakpriors).

## Other essential settings and data to provide
These settings are provided as targets in the scripts. Here, they are specified following the scheme "name : type : default : explanation"
- "threshold_dbh"	:	``numeric``		:	200			:	dbh threshold between A and B in [mm]
- "offsetname"		:	``character``	:	"offset"	:	column name of the value used as an offset in the JAB model
- "crs"				:	``character``	: 				: 	if "Stages" is not provided as sf object but in some table format, please provide the crs as an EPSG code


---------------


## Notes for fitting a new data set

### Target adaptions
Change targets:
- ``dir_publish``
- ``dir_fit``
- In addition, rename all targets.

### Target interfaces
- The newly included table Stages should be treated like the target ``Stages_env``, and split by ``taxon_s``, then fed into  ``tar_target(fits_s, fitS(BA_s), pattern = map(BA_s), iteration = "list")``, and then into ``tar_target(Stages_s, predictS(fits_s, Stages_env), iteration = "list")``
- Both, the above and the newly included table ``Stages_transitions`` can be fed into ``tar_target(data_stan, formatStanData(Stages_scaled, Stages_transitions, taxon_s, threshold_dbh, timestep = 1, parfactor = 1))``
