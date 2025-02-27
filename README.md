# The Nitrogen Risk Assessment Model for Scotland (NIRAMS II)

NIRAMS II is a simple, spatially distributed nitrate leaching model designed to estimate nitrate loadings to groundwater from the base of the soil profile. 
The model was developed by the [James Hutton Institute](https://www.hutton.ac.uk/) and has been used to support the Scottish Government's designations and reviews of [Nitrate Vulnerable Zones (NVZ)](https://www.gov.scot/policies/agriculture-and-the-environment/nvz/) in line with the EU Nitrates Directive.


## Model documentation

The model is coded in Python and can be found [here](https://github.com/JamesSample/nirams_ii), where the background to the model's development, its uses and some basic documentation are also provided.

## Deriving input files for NIRAMS

This repository contains **R scripts** for processing input data and deriving the required input HDF5 file for NIRAMS. To run the NIRAMS II model, the user must first obtain such input HDF5 file containing gridded datasets of meterological time series, soil properites and agricultural activities. Unfortunatey, this file cannot be made freely available, as several of the datasets included in it were obtained under a license that precludes redistribution.

NIRAMS II makes use of the following spatial input datasets:
* **Climate grids** from the UK Met Office [HadUK-Grid](https://www.metoffice.gov.uk/research/climate/maps-and-data/data/haduk-grid/haduk-grid), with potential evapotranspiration estimated using the [FAO56 modified Penman-Monteith methodology](https://www.fao.org/4/x0490e/x0490e00.htm). The gridded HadUK climate data are available as netCDF files and need to be downloaded prior to the processing here. For the work here, climate data have been downloaded at 5km resolution and then re-sampled to 1km resolution.
* **Soil properties** from the [Hydrology of Soil Types (HOST) database](https://nora.nerc.ac.uk/id/eprint/7369/). 
* **Crop type, land use and livestock** data from [EDINA AgCensus the Scottish](https://agcensus.edina.ac.uk/). The gridded AgCensus data (at 2km resolution) are available as raster files and need to be downloaded prior to the processing here. The AgCensus data are re-samled to 1km resolution.
* **Land cover** data from [UKCEH Land Cover Maps (LCM)](https://www.ceh.ac.uk/data/ukceh-land-cover-maps). The LCM data are available as raster files and have for the work here been downloaded at 1km resolution.
* **Atmospheric deposition estimates** from the [European Monitoring and Evaluation Programme (EMEP)](https://www.emep.int/). The deposition data from EMEP is available at 0.1x0.1 degree (long,lat) resolution and available as NetCDF files. The files need to be downloaded prior to the processing here. The data is re-sampled to 1km resolution.

The AgCensus data are used (together with LCM) to derive annual estimates of organic and inorganic nitrogen inputs as well as nitrogen plant uptake. AgCensus contain gridded data (at 2 km resolution) on the distribution of land uses/crop types as well as the number, age, and type of livestock. The amount of organic nitrogen excreted annually by each animal class was taken from manure planning documentation issued to farmers within Scotland's NVZs. Annual estimates for the application rate of inorganic nitrogen fertiliser are informed by the results of the British Survey of Fertiliser Practice. The annual estimates of organic and inorganic nitrogen application as well as of nitrogen uptake are distributed spatially at 1 km resolution based on appropriate land classes from LCM. These annual estimates are then distributed temporally using a set of idealised time series that define, for a variety of crop classes, the length of the growing season, the amount of nitrogen uptake and the timing of fertiliser application (see details in Model documentation above). 

Note that previous versions of NIRAMS (before 2013) used field-scale land use and crop data from IACS and business-scale data on livestock from June Agricultural Census (JAC) to derive the annual estimates of organic N, inorganic N, and N uptake. This part has since been replaced with AgCensus data, partly because of ease (the AgCensus data are already gridded and include detailed information on both crops and livestock numbers) and partly because of issues with accessibility to IACS data. 

