# Gapless Daily Mean LST on Google Earth Engine

## Overview

This repository provides Google Earth Engine (GEE) JavaScript code for estimating global gapless daily mean land surface temperature (LST) and generating associated pixel-level quality information.

## Gapless Daily Mean LST

The script `reconstruct_gapless_daily_mean_lst_with_eatc.js` first estimates gapped daily mean LST using linear combinations of two to four valid Terra and Aqua MODIS instantaneous LST observations. Temporal gaps are subsequently reconstructed using the enhanced annual temperature cycle (EATC) model. For tropical regions, outlier ATC parameters and noise were filtered, and the remaining gaps were filled using spatial interpolation.

## Pixel-Level Quality Information

The script `encode_modis_observation_temporal_upscaling_qa.js` encodes information on the four nominal daily MODIS observations, the number of valid daily observations, and the validity of the linear combinations performed in Step 1. For each pixel-day, this information is stored as a single 16-bit unsigned integer using the following encoding scheme:

| Bits | Field | Encoding |
|:---:|---|---|
| 2–0 | Terra daytime flag | `000` = invalid observation<br>`001` = average LST error ≤ 1 K<br>`010` = 1 K < average LST error ≤ 2 K<br>`011` = 2 K < average LST error ≤ 3 K<br>`100` = average LST error > 3 K |
| 5–3 | Terra nighttime flag | Same encoding as above |
| 8–6 | Aqua daytime flag | Same encoding as above |
| 11–9 | Aqua nighttime flag | Same encoding as above |
| 14–12 | Number of valid daily MODIS observations | `000` = 0 valid observations<br>`001` = 1 valid observation<br>`010` = 2 valid observations<br>`011` = 3 valid observations<br>`100` = 4 valid observations |
| 15 | Step 1 temporal-upscaling flag | `0` = invalid<br>`1` = valid |

The script `calculate_eatc_fitting_rmse.js` calculates the pixel-level root-mean-square error (RMSE) between the input LST values used for EATC fitting and the corresponding fitted values. 



## References

Xing, Z., et al. (2021). Estimation of daily mean land surface temperature at global scale using pairs of daytime and nighttime MODIS instantaneous observations. *ISPRS Journal of Photogrammetry and Remote Sensing, 178*, 51–67. https://doi.org/10.1016/j.isprsjprs.2021.05.017

Zou, Z., et al. (2018). Enhanced modeling of annual temperature cycles with temporally discrete remotely sensed thermal observations. *Remote Sensing, 10*, 650. https://doi.org/10.3390/rs10040650
