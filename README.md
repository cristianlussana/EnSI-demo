# EnSI-demo

This repository contains the supporting information for the scientific article: 
Ensemble-based statistical interpolation of atmospheric variables near the surface C Lussana, TN Nipen, B Menetrier and IA Seierstad (Almost) Ready for submission to Quarterly Journal of the Royal Meteorological Society

## Features
- EnSI (Ensemble-based statistical interpolation) for post-processing numerical model output by combining it with in-situ observations using statistical techniques. EnSI extends Optimal Interpolation (OI) and Ensemble OI—traditionally used in spatial analysis for individual variables and in a time-independent manner—by incorporating consecutive observation times and cross-correlations between variables.
- started Box-Cox data transformation for precipitation to enhance spatial analysis performance.
- EnSI is demonstrated through the reconstruction of hourly precipitation and temperature fields over:
  - Oslo (ensi-metno_0d);
  - One-dimensional cross section (ensi-metno_1d);
  - Two-dimensional regular grid over Scandinavia (ensi-metno_2d).
- EnSI configuration files for all the experiments described in the scientific article.
  

## Resources
For information on how to use EnSI, check out the wiki at https://github.com/cristianlussana/EnSI-demo/wiki.

## Copyright and license
Copyright © 2024-2024 Norwegian Meteorological Institute. EnSI-deom is licensed under the GNU General Public License v3.0. See LICENSE file.
