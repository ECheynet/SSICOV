# Operational modal analysis with automated SSI-COV algorithm
The modal parameters of a line-like structure are automatically identified using an SSI-COV algorithm applied to ambient vibration data

[![View Operational modal analysis with automated SSI-COV algorithm on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/69030-operational-modal-analysis-with-automated-ssi-cov-algorithm)
[![DOI](https://zenodo.org/badge/248938005.svg)](https://zenodo.org/badge/latestdoi/248938005)
<a href="https://www.buymeacoffee.com/echeynet" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 25px !important;width: 120px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>

## Summary

The function SSICOV.m aims to automatically identify the eigenfrequencies, mode shapes and damping ratios of a line-like structure using ambient vibrations only. The covariance-driven stochastic subspace identification method (SSI-COV) is used in combination with a clustering algorithm to automatically analyse the stabilization diagrams. 

The algorithm is inspired by the one used by Magalhaes et al. [1]. It has been applied for ambient vibration monitoring of the Lysefjord Bridge [2] and was compared to the frequency domain decomposition technique [3]. Finally, the algorithm was found accurate enough to visualise the evolution of the bridge eigenfrequencies with the temperature [4].

## content

The submission file contains:
- A data file BridgeData.mat
- A Matlab Live Script Example1.mlx that illustrates the application of the algorithm.
- A Matlab Live Script Example1_noToolbox.mlx that reproduce Example1 but using the function SSICOV_noToolbox.
- The function SSICOV which is the automated SSI-COV algorithm.
- The function SSICOV_noToolbox which is the automated SSI-COV algorithm but does not use the Statistics and Machine Learning Toolbox. The Linkage algorithm is replaced by the function "PHA_Clustering" by [5] and the function "cluster" is replaced by the function "Cluster2", which is derived from [6].
- The function plotStabDiag.m, which plot the stabilization diagram.

Any question, suggestion or comment is welcomed.

## References

[1] Magalhaes, F., Cunha, A., & Caetano, E. (2009). Online automatic identification of the modal parameters of a long span arch bridge. Mechanical Systems and Signal Processing, 23(2), 316-329.

[2] Cheynet, E., Jakobsen, J. B., & Snæbjörnsson, J. (2016).Buffeting response of a suspension bridge in complex terrain. Engineering Structures, 128, 474-487.

[3] Cheynet, E., Jakobsen, J. B., & Snæbjörnsson, J. (2017).Damping estimation of large wind-sensitive structures.Procedia Engineering, 199, 2047-2053.

[4] Cheynet, E., Snæbjörnsson, J., & Jakobsen, J. B. (2017).Temperature Effects on the Modal Properties of a Suspension Bridge.In Dynamics of Civil Structures, Volume 2 (pp. 87-93). Springer.

[5] Yonggang (2021). Fast hierarchical clustering method - PHA (https://www.mathworks.com/matlabcentral/fileexchange/46134-fast-hierarchical-clustering-method-pha), MATLAB Central File Exchange. Retrieved February 4, 2021. 

[6]  Eric Ogier (2021). Hierarchical clustering (https://www.mathworks.com/matlabcentral/fileexchange/56844-hierarchical-clustering), MATLAB Central File Exchange. Retrieved February 4, 2021. 
