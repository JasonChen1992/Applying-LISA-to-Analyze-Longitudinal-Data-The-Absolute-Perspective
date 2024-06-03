# Applying Local Indicators of Spatial Association to Analyze Longitudinal Data: The Absolute Perspective

The traditional LISA methods are enhanced in this project to handle longitudinal data, allowing for the analysis of changes in spatial patterns over time. This approach helps to reveal absolute distributional dynamics, providing a more comprehensive view of spatial processes by Ran Tao and Yuzhou Chen.

## Abstract
Local Indicators of Spatial Association (LISA) are a class of spatial statistical methods that have been widely applied in various scientific fields. When applying LISA to make longitudinal comparisons of spatial data, a common way is to run LISA analysis at each time point, then compare the results to infer the distributional dynamics of spatial processes. Given that LISA hinges on the global mean value that often varies across time, the LISA result generated at time Ti reflects the spatial patterns strictly with respect to Ti. Therefore, the typical comparative cross-sectional analysis with LISA can only characterize the relative distributional dynamics. However, the relative perspective alone is inadequate to comprehend the full picture, as the patterns are not directly associated with the changes of the spatial process’s intensity. We argue that it is important to obtain the absolute distribution dynamics to complement the relative perspective, especially for tracking how spatial processes evolve across time at the local level. We develop a solution that modifies the significance test when implementing LISA analysis of longitudinal data to reveal and visualize the absolute distribution dynamics. Experiments were conducted with Mongolian livestock data and Rwanda population data.

## Methodology

### Local G Statistic & Local Moran's I

Longitudinal Analysis
The traditional LISA methods are enhanced in this project to handle longitudinal data, allowing for the analysis of changes in spatial patterns over time. This approach helps to reveal absolute distributional dynamics, providing a more comprehensive view of spatial processes.

## How to Use
1. Clone the repository to your local machine.
2. Ensure you have the necessary dependencies installed: `clusterpy`, `numpy`, `shapefile`, and `pandas`.
3. Place your input data files in the appropriate directory.
4. Run the scripts to calculate absolute perspective of Local G(*) and Local Moran's I
5. The results will be saved in a CSV file.


References
For more detailed information on the methodology and applications, refer to the following paper:

Ran Tao, Yuzhou Chen. "Applying Local Indicators of Spatial Association to Analyze Longitudinal Data: The Absolute Perspective." Geographical Analysis (2022). DOI: 10.1111/gean.12323

Contact
For any questions or issues, please contact Yuzhou Chen at yuzhouchen@usf.edu or Ran Tao at rtao@usf.edu.

