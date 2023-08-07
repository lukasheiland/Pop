# Pop – development of a simple forest population model

This is a project for developing a forest population model that resulted in the JAB model. The JAB model is stored in `Model_2023-02_bb/` (for more see Section 'File structure').


## File structure

* Directories beginning with `Model_` contain versions of the population model, `Model_2023-02_bb/` being the published JAB model.
* The root directory contains all the R scripts for processing the data, test simulations, and fitting the models--including for old and discarded model versions.
* In particular, the files `Main….R`, call the targets pipelines formulated in `_targets.R`
	* `Main.R` calls the pipeline published in Heiland et al. (2023).	
	* `Main_env.R` calls the pipeline with environmenntal response of the rates.
	* `Main_test.R` calls a test pipeline.
* Files `…_functions.R` contain the actual functions that are called in the targets pipelines.
* `Inventory.nosync/` contains the main data, including the NFIs and environmental data. The directory is a subset of the data organized in the 'Inventories' project, git branch `pop`.
* `Data/` contains data processed or downloaded by the pipeline.
* `Range.nosync/` contains only data used in comparing the environmental distributions of the German NFI and Fagus sylvatica in the EAFTS.
* `Data model/` contains the data model.
* `Fits.nosync/`, `Publish.nosync/` contain the **generated results**.
	* `Fits.nosync/` the MCMC draws plus some infos on the fit including the offset used
		* (The info on which run corresponds to which offsets is stored in file names of empty txt files within "Fits.nosync")
	+ `Publish.nosync/` figures, tables



-------


## Inventory data

The data in `Inventory.nosync/` is a curation of the German and the Slovakian national forest inventories (NFI) for use with the JAB model, a size-structured dynamic forest population model.

### Description of the data and file structure

The directory `Inventory.nosync` contains all data needed for projects with the JAB model. Details are presented below for the different data files.

The main software project folder, e.g., `Pop_bb_publish` at https://doi.org/10.5281/zenodo.8032461, is intended to contain the data folder Inventory.nosync (currently contains a placeholder).

The sub-directory `Fits.nosync` will contain the generated fits, while `Publish.nosync` will contain the generated plots and tables.


#### Details for: `DE BWI/Data/DE_BWI_big_abund.rds`
* Description: Abundances of sample trees ('large trees') in the German NFI from constant angle sampling 
    * The preprocessed inventory data has the following format:
        * all forest plots are represented with at least one row per species,
        * plots without observation have one row per plot with NAs.
* Columns:
    * taxid: id of the taxon
    * tax: short clear text taxon
    * treeid: id of the tree
    * dbh: diameter at breast height in mm
    * count: count with respect to the given countarea
    * countarea: area to which the count refers
    * age: age of a tree, guesstimate from the German NFI
    * height: height of a tree
    * forestrydamage: records of damage by forestry derived from different categories of the German NFI
    * dbhLower: lower dbh threshold of the inventory
    * plotid: id of the forest plot
    * plotobsid: id of forest plot observations, interaction of NFI, year, and plot
    * methodid: id of a consistent method, can include multiple surveys (obsid)
    * envjoinid: id to join spatial data, usually plot or cluster
    * inventory: inventory, e.g. a certain NFI
    * time: date of observation


#### Details for: `DE BWI/Data/`
* Description: Additional status information of sample trees ('large trees') in the German NFI from constant angle sampling. Status info includes whether dead or alive and the distance to the plot center.
* Columns:
    * plotid: id of the forest plot
    * treeid: id of the tree
    * obsid: id of the survey, interaction of NFI and year
    * distance: distance to the plot center of the constant angle count method in cm
    * dbh: diameter at breast height in mm
    * count_ha: count per hectare (usually 1/area)
    * kennziffer: NFI internal id
    * alive: logical, whether tree is dead
    * dead: logical, whether tree is alive
    * harvested: logical, whether tree was harvested
    * samplingissue: logical, was there an issue with sampling?
    * time: date of observation
    * clusterid: id of a cluster of plots (in the German NFI 'Trakt')
    * taxid: id of the taxon
    * tax: short clear text taxon



#### Details for: `DE BWI/Data/DE_BWI_small_abund.rds`
* Description: Counts of small trees ('saplings') in size classes with corresponding sample areas in the German NFI.
    * The preprocessed inventory data has the following format:
        * all forest plots are represented with at least one row per species,
        * plots without observation have one row per plot with NAs.
* Columns:
    * plotid: id of the forest plot
    * taxid: id of the taxon
    * tax: short clear text taxon
    * regclass: size class of saplings
    * regclasstype: "h" height, "d" diameter, or "hd" lower height and upper diameter
    * regclassLower: lower threshold in cm
    * regclassUpper: upper threshold in cm
    * count: count with respect to 'countarea'
    * countarea: area for the counts within a certain size class on a forest plot
    * obsid: id of the survey, interaction of NFI and year
    * plotobsid: id of forest plot observations, interaction of NFI, year, and plot
    * methodid: id of a consistent method, can include multiple surveys (obsid)
    * envjoinid: id to join spatial data, usually plot or cluster
    * inventory: inventory, e.g. a certain NFI
    * time: date of observation



#### Details for: `DE BWI/Data/DE_BWI_Env_sf.rds`
* Description: Simple feature (sf) data frame of spatial information including environmental variables for each forest plot of the German NFI.
* Columns:
    * plot: NFI internal id of a forest plot within a cluster
    * plotid: id of the forest plot
    * clusterid: id of a cluster of forest plots. Clusters include up to 4 plots in the German NFI
    * slope_loc: slope data from the German NFI
    * aspect_loc: aspect data from the German NFI
    * alt_loc: elevation data from the German NFI in m
    * time_predictors: the time the environmental variables refer to; empty
    * nyears_predictors: the number of years the environmental variables refer to; empty
    * geometry: coordinates in sf format
    * inventory: inventory, e.g. a certain NFI
    * envjoinid: id to join spatial data, usually plot or cluster
    
    * phCaCl_esdacc:
    * waterLevel_loc:
    
    * regType_DE_BWI_2: regeneration type in the second survey (BWI 2)
    * allNatRegen_DE_BWI_2: plot-level aggregated logical whether all regeneration type was natural in the second survey (BWI 2)
    * anyUnnatRegen_DE_BWI_2: plot-level aggregated logical whether any regeneration type was unnatural in the second survey (BWI 2)
    * anyHarvested_DE_BWI_2: plot-level aggregated logical whether there was any harvested tree in the second survey (BWI 2)
    * anyForestryDamage_DE_BWI_2: plot-level aggregated logical whether there was any damage from forestry in the second survey (BWI 2)
    * regType_DE_BWI_3: as above but for the third survey (BWI 3)
    * allNatRegen_DE_BWI_3: as above but for the third survey (BWI 3)
    * anyUnnatRegen_DE_BWI_3: as above but for the third survey (BWI 3)
    * anyHarvested_DE_BWI_3: as above but for the third survey (BWI 3)
    * anyForestryDamage_DE_BWI_3: as above but for the third survey (BWI 3)

    

#### Details for: `DE BWI/Data/DE_BWI_geo.rds`
* Description: Simple feature (sf) data frame of geographical information for clusters of plots in the German NFI. This includes information to which grid of plots a cluster ('Trakt') belongs. See German NFI for more information.

* Columns:
    * Tnr: NFI internal cluster id
    * Netz: which grid does a cluster belong to?
    * geometry: coordinates in sf format
    
    
    
#### Details for: `DE BWI/Data/Waterlevels.csv`
* Description: A pseudo-ratio ranking of hydromorphic soil water levels on plots of the German NFI. The ranking is intended to translate levels as recorded by https://www.openagrar.de/receive/openagrar_mods_00049873?cmdf=benning+2019+umweltdatenbank+bodenprofile into the ratio-scale.
* Columns:
    * Icode: id in the German NFI database
    * Acode: factor level in the German NFI database, labels correspond to hydromorphic soil water levels
    * KurzD: factor level in the German NFI database, labels correspond to hydromorphic soil water levels
    * KurzE: factor level in the German NFI database
    * LangD: long description of the hydromorphic soil water levels in German (translations of relevant labels in column Ellenberg_EN)
    * zu.whhform: id in the German NFI database
    * waterlevel_loc: pseudo-ratio ranking
    * Ellenberg: levels corresponding to the Ellenberg ecogram
    * Ellenberg_EN: English translation of levels corresponding to the Ellenberg ecogram
    
    
    
#### Details for: `SK NIML/Data/SK_NIML_complete.rds`
* Description: Tree abundances from the Slovakian NFI. The data includes large sample tree measurements on semi-constant plot areas and counts of seedlings and saplings of one survey in the years 2015–2017.
* Columns:
    * plotid: id of the forest plot
    * year: year of the survey
    * taxon: long taxon namme
    * sizeclass: "big" measured sample trees with dbh, or "small" counted seedlings and saplings
    * height: tree height
    * area: area to which the counts refer
    * count: tree counts within taxon, sizeclass, plotid, year with respect to 'area'
    * dbh: diameter at breast height
    * regeneration: plot level variable ('Artificial', 'Natural - coppice', or 'Natural - seed')
    * count_ha: counts per hectare
    * treeid: id of a tree
    * ba: basal area of a tree
    * ba_ha: basal area per hectare of a tree (varying plot area)
    * status: 'Damage stem', 'Dead tree', or 'Live tree' 
    * dbh_lower: lower dbh threshold for measured sample trees of the NFI
    * management_plot: management type at plot level ('Clear cut', 'No management', 'Selection', or 'Shelter')
    * WGS_E: easting
    * WGS_N: northing
    * isforest: whether a plot is forest. This is important for `SK NIML/Data/SK_NIML_complete_fullgrid.rds`, where plot grid locations without forest are included for a complete regular grid.


#### Details for: `SK NIML/Data/SK_NIML_complete_fullgrid.rds`
* Columns:
    * columns are identical to `SK NIML/Data/SK_NIML_complete.rds`

