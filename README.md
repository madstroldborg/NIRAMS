# The Nitrogen Risk Assessment Model for Scotland (NIRAMS II)

NIRAMS II is a simple, spatially distributed nitrate leaching model designed to estimate nitrate loadings to groundwater from the base of the soil profile. 
The model was developed by the [James Hutton Institute](https://www.hutton.ac.uk/) and has been used to support the Scottish Government's designations and reviews of [Nitrate Vulnerable Zones (NVZ)](https://www.gov.scot/policies/agriculture-and-the-environment/nvz/) in line with the EU Nitrates Directive.


## Model documentation

The model is coded in Python and can be found [here](https://github.com/JamesSample/nirams_ii), where the background to the model's development, its uses and some basic documentation are also provided.

## Deriving input files for NIRAMS

This repository contains R scripts for processing all the necessary input data and deriving the required input HDF5 file for NIRAMS. To run the NIRAMS II model directly, the user must first obtain the input HDF5 file containing gridded datasets of meterological time series, soil properites and agricultural activities. Unfortunatey this file cannot be made freely available, as several of the datasets included in it were obtained under a license that precludes redistribution.
