2019-thunderstorms-analysis

# Lagrangian Analysis of Thunderstorms in Switzerland
## Analysis Tools for Investigation of Thunderstorm Initiation Conditions

Author:

Thomas M. Lanz | 14-119-564 | MSc in Climate Sciences | OCCR - University of Bern | thomas.lanz@students.unibe.ch

## Introduction and Research Question
The aim of the analysis tools is the investigation of atmospheric conditions for thunderstorm initiation. Such atmospheric conditions for the evolution of deep convective cells are an unstable stratification of the atmosphere, a substantial amount of ground level moisture and an instability triggering process like convergence or lifting (Wallace & Hobbs, 2006). These conditions represent the basic ingredients for thunderstorm initiation. In general, a thunderstorm is defined by the American Meteorological Society (AMS) as '[...]a local storm, invariably produced by a cumulonimbus cloud and always accompanied by lightning and thunder, usually with strong gusts of wind, heavy rain, and sometimes with hail' (AMS, 2012). Despite their frequent occurrence and damage potential (MeteoSwiss, 2018; Nisi et al., 2016; Trefalt et al., 2018), the initiation of thunderstorms is still incompletely understood. The goal of the developed analysis tools is to fill the above mentioned research gap and to answer the following research question: 

How far can the analysis tools unveil the atmospheric conditions and processes responsible for thunderstorm initiation?

This coding project is conducted in the framework of the seminar 'Geodata Analysis and Modelling' (FS2019-438745) and is settled within a master thesis project.

## About the Analysis Tools
For reaching the prediscribed aim of this project, the analysis tools consists of five different tools for investigating the condition and processes of the atmosphere: horizontal maps, soundings, vertical cross-section, maps of trajectories and temporal evolution along trajectories. The figures created with these tools serve as basis for the analysis of the atmosphere and the answering of the research question.

The programming was done in IDE Jupyter Notebook (v6.0.3) with the IPython Kernel (v7.13.0) and with package and environment management by Conda. For each of the five analysis tools, a notebook document was created (also python files are provided). In the following section the workflow of the analysis tools is described.

## What the Analysis Tools Do and How to Use
For the analysis tools data is used from the Weather Research and Forecasting (WRF) model (v4 ARW) (Powers et al., 2017) and from the Lagranto program (Sprenger & Wernli, 2015), which takes input data from the WRF model. The WRF data is in netCDF and the Lagranto output data in ASCII file format. 

The following subsections provide information regarding the useage of the codes. For background information of the scientific methods, see the file scientific_methods.pdf. Example figures of each analysis tool are presented in the subsequent Results section.

### Horizontal Maps
This analysis tool produces horizontal (2D) maps of different variables. Supported variables for the plotting function (horizontal_map) are updraft, reflectivity, helicity, pw, cape, cin, ctt, temperature_surface, wind_shear, rh, omega, pvo, avo, theta_e, water_vapor, uv_wind, updraft_reflectivity and divergence. Besides the selection of the desired variable_name for the plotting function, more input parameters need to be defined like date, start_hour, end_hour, pressure_level, subset, initiation, save, gif. The first five parameters (variable_name included) need a specific input value (e.g. variable_name="divergence", date="2018-05-30", start_hour=16, end_hour=17, pressure_level=850), where the remaining parameters need Boolean values (True or False). Further, the definition of some more variables in first section of the code is necessary, before using the plotting function (data_directory, save_dir, subset_extent and initiation_location). The horizontal_map function iterates with a 5 minutes time step over a list of files from start_hour until end_hour. If the gif parameter is set True, a gif is created from all the generated figures and the single files are deleted.

### Soundings
The sounding is an analysis tool for investigating the vertical distribution of atmospheric physical properties (e.g. temperature, pressure, wind, etc.) and represents the WRF model data in a similar way (thermodynamic diagram) like measurement data from an real world atmospheric sounding (e.g. balloon sounding). For a selected location (lat, lon) and time (date), the analysis tool generates a skew T-log p diagram, based on the variables derived from WRF data. Because the WRF model data lacks some of these specific variables (e.g. pressure, dewpoint temperature or wind speed), these required variables need to be computed (by wrf-python function getvar()) and added back to the WRF dataset. Afterwards, the variables are selected for a specific location and some variables need to be calculated by MetPy functions. Finally, the code generates a figure according to the specific layout of a skew T-log p diagram (see Results).

### Vertical Cross-Section
Vertical cross-sections show a vertical slice of the atmosphere along a line with specified start (start_lat, start_lon) and end points (end_lat, end_lon). The analysis tool is represented by a plotting fucntion (cross_section), which supports the following variables: vertical_velocity, rh, omega, absolute_vorticity, theta_e and reflectivity. The only input parameters left to define are date, time and save, if a saving of the figure is desired (default save=False). Before the cross_section function can be used, the data directory need to be adjusted. After interpolation and removing of white space between terrain height and contour of the variable, a vertical cross-section figure is generated with filled mountain area.

### Maps of Trajectories
The analysis tool for mapping trajectories uses the output data (trajectories) from the Lagranto program. With the help of the lagranto_plotting function, the desired variables (water_vapor, updraft and height) along the trajectories can be plotted on the background of terrain height contours (greys). The start_time and end_time of the calculated trajectories need to be indicated as input parameters, as well as Booloean value True, if a subset or a saving of the figure is requested. Furthermore, a bunch of trajectories can be selected according to their height level (pbl, 5 or 10) or otherwise, all available trajectrories will be plotted. For defining the Planteary Boundary Layer (PBL) height, some more variables need to be set. The number of plotted trajectories in the figure, can be adjusted by specifying number_trajs_plot. Now all necessary input parameters for the function are defined and a horizontal map of trajectories can be made.

### Temporal Evolution along Trajectories
This analysis tool creates a kind of cross-section, but in this case with time since initiation on the x-axis and the values of the chosen variable along the trajectories on the y-axis. As already mentioned in the section before, the temporal_evolution_trajectories is also capable of seperation the trajectories in bunches of different heights (trajs_bunch=', default = 'all'). Besides the data directory of the trajectory data and the desired saving directory, time delta needs to be specified too and the number of plotted trajectories in the figure can be adjusted as well. Further, some variables for getting the PBL height need to be set. For the analysis and comparison between different vertical trajectory bunches, for each bunch of trajectories the 10th and 90th percentile and mean is computed and inidicated with colored lines on the figure and on the added legend. With all needed variables defined, the function is ready for creating a figure of the temporal evolution along trajectories for the selected variable.

## Results
This section shows selected results of the respective analysis tools.

### Horizontal Maps

### Soundings

### Vertical Cross-Section

### Maps of Trajectories

### Temporal Evolution along Trajectories

## Thanks

## References
AMS. (2012, April 25). Thunderstorm - AMS glossary. Retrieved April 2, 2019, from http://glossary.ametsoc.org/wiki/Thunderstorm

MeteoSwiss. (2018, December 21). 2018: Rekordwärme und massive regenarmut. Retrieved March 23, 2020, from https://www.meteoschweiz. admin.ch/home/suche.subpage.html/de/data/blogs/2018/12/ 2018-waermstes-jahr-seit-messbeginn-.html?query=sommer+2018+ schweiz&topic=0&pageIndex=0&tab=search tab

Nisi, L., Martius, O., Hering, A., Kunz, M., & Germann, U. (2016). Spatial and temporal distribution of hailstorms in the alpine region: A long-term, high resolution, radar-based analysis. Quarterly Journal of the Royal Meteorological Society, 142(697), 1590–1604. https://doi.org/10.1002/ qj.2771

Powers, J. G., Klemp, J. B., Skamarock, W. C., Davis, C. A., Dudhia, J., Gill, D. O., Coen, J. L., Gochis, D. J., Ahmadov, R., Peckham, S. E., Grell, G. A., Michalakes, J., Trahan, S., Benjamin, S. G., Alexander, C. R., Dimego, G. J., Wang, W., Schwartz, C. S., Romine, G. S., . . . Duda, M. G. (2017). The weather research and forecasting model: Overview, system efforts, and future directions. Bulletin of the American Meteo- rological Society, 98(8), 1717–1737. https://doi.org/10.1175/BAMS-D- 15-00308.1

Sprenger, M., & Wernli, H. (2015). The LAGRANTO lagrangian analysis tool – version 2.0. Geoscientific Model Development, 8(8), 2569–2586.

Trefalt, S., Martynov, A., Barras, H., Besic, N., Hering, A. M., Lenggenhager, S., Noti, P., Röthlisberger, M., Schemm, S., Germann, U., & Martius, O. (2018). A severe hail storm in complex topography in switzerland - observations and processes. Atmospheric Research, 209, 76–94. https: //doi.org/10.1016/j.atmosres.2018.03.007

Wallace, J. M., & Hobbs, P. V. (2006, March 24). Atmospheric science: An introductory survey [Google-Books-ID: HZ2wNtDOU0oC]. Elsevier.
