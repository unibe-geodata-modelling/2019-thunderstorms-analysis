2019-thunderstorms-analysis

# Lagrangian Analysis of Thunderstorms in Switzerland
## Analysis Tools for Investigation of Thunderstorm Initiation Conditions

Author:

Thomas M. Lanz | 14-119-564 | MSc in Climate Sciences | OCCR - University of Bern | thomas.lanz@students.unibe.ch

## Introduction and Research Question
The aim of these tools is the analysis of atmospheric conditions for thunderstorm initiation. Such atmospheric conditions for the evolution of deep convective cells are an unstable stratification of the atmosphere, a substantial amount of ground level moisture and an instability triggering process like convergence or lifting (Wallace & Hobbs, 2006). These conditions represent the basic ingredients for thunderstorm initiation. In general, a thunderstorm is defined by the American Meteorological Society (AMS) as '[...]a local storm, invariably produced by a cumulonimbus cloud and always accompanied by lightning and thunder, usually with strong gusts of wind, heavy rain, and sometimes with hail' (AMS, 2012).

During the summer months in Switzerland, thunderstorms occur frequently and cause a lot of damage (MeteoSwiss, 2018; Nisi et al., 2016; Trefalt et al., 2018). Although, they appear in Switzerland on small spatial scales (compared to e.g. a foehn storm), they have a high damage potential by hail, lightning, wind gusts and high precipitation amounts in short time periods (Doswell et al., 1996; García-Ortega et al., 2007; Trefalt et al., 2018). This reveals that in-depth understanding of these destructive weather events is of great importance and an accurate forecast a desirable target. Despite their frequent occurrence and damage potential, the initiation of thunderstorms is still incompletely understood. The aim of developing the analysis tools is to fill the above mentioned research gap and to answer the following research question: 

Which atmospheric conditions and processes are responsible for thunderstorm initiation in May 2018 in Switzerland?

This coding project is conducted in the framework of the seminar 'Geodata Analysis and Modelling' (FS2019-438745) and is settled within a master thesis project.

## About the Analysis Tools
For reaching the prediscribed aim of this project, the analysis tools consists of five different tools for investigating the state of the atmosphere: Horizontal maps, soundings, vertical cross-section, maps of trajectories and temporal evolution along trajectories. Based on these tools, figures are created, which serve as basis for analysis of the atmospheric conditions.

The programming was done in IDE Jupyter Notebook (v6.0.3) with the IPython Kernel (v7.13.0) and with package and environment management by Conda. For each of the five analysis tools, a notebook document was generated, but also a python file is provided. In the following section the workflow of the analysis tools is described.

## What the Analysis Tools Do
For the analysis tools is data used from the Weather Research and Forecasting (WRF) model (v4 ARW) (Powers et al., 2017) and from the Lagranto program (Sprenger & Wernli, 2015), which input data is also from the WRF model. The WRF data is in netCDF and the Lagranto output data in ASCII file format. 

In the following subsections, the codes of the analysis tools are presented and information regarding how to use the codes is provided. For background information of the scientific methods, see the file xy.

### Horizontal Maps
This analysis tool produces horizontal (2D) maps of different variables. Supported variables for plotting procedure are: updraft, reflectivity, helicity, pw, cape, cin, ctt, temperature_surface, wind_shear, rh, omega, pvo, avo, theta_e, water_vapor, uv_wind, updraft_reflectivity and divergence. Besides the selection of the desired variable_name for the plotting function 'horizontal_map', there are more input parameters to be defined: date, start_hour, end_hour, pressure_level=False, subset=False, initiation=False, save=False, gif=False. The first five parameters (variable_name included) need a specific input value (e.g. variable_name="divergence", date="2018-05-30", start_hour=16, end_hour=17, pressure_level=850), where the remaining paramteres need Boolean values (True or False). Further, some more variables need to be defined before using the plotting function (e.g. data_directory, initiation_location or variables depending on the chosen variable_name, see first section of function). Afterwards, the horizontal_map plotting function is ready to usage and generates figures according to the chosen input parameters. For some examples of different types, see Results section.

### Soundings
The sounding is an analysis tool for investigating the vertical structure of the atmosphere for a selected location and date. At this point in space and time, the analysis tool generates a skew T-log p diagram. First, the required variables (e.g. pressure or temperature) need to be computed (by wrf-python function getvar()) and then added back to the WRF dataset. Afterwards, the variables are selected for a specific location and then plotted according to the specific layout of a skew T-log p diagram (see Results).

### Vertical Cross-Section
Vertical cross-sections show a vertical slice of the atmosphere along the specified start (start_lat, start_lon) and end points (end_lat, end_lon). The analysis tool is represented by a plotting fucntion (cross_section), which supports the following variables: vertical_velocity, rh, omega, absolute_vorticity, theta_e and reflectivity. The only input parameters left to define are date, time and save (only necessary if save=True, default is save=False). Based on the defined input parameters, a vertical cross-section is generated with filled mountain area.

### Maps of Trajectories
The analysis tool for mapping trajectories uses the output data (trajectories) from the Lagranto program. With the help of a function (lagranto_plotting) the desired variables (water_vapor, updraft and height) along the trajectories can be plotted on the background of terrain height contours (grey). The start_time and end_time of the calculated trajectories need to be indicated as input parameters, as well as Booloean value True if a subset or a svaing of the figure is required. A bunch of trajectories can be selected according to their height level (pbl, 5 or 10). Now all necessary input parameters for the function are available and a horizontal map of trajectories can be made. 

### Temporal Evolution along Trajectories
This analysis tool creates a kind of cross-section, but in this case with time on the x-axis and the chosen variable along the trajectories on the y-axis. As mentioned in the section before, also this function (temporal_evolution_trajectories) can seperate the trajectories in bunches of different heights (trajs_bunch=', default = 'all'). Besides the data directory of the trajectory data and the desired save directory, also the frequency of trajectories can be adjusted and the time delta delta needs to be specified. Further, some variables for getting the PBL height need to be set. With all needed variables defined, the function is ready for creating a figure of the temporal evolution along trajectories for the selected variable.

## Results

## Thanks

## References
AMS. (2012, April 25). Thunderstorm - AMS glossary. Retrieved April 2, 2019, from http://glossary.ametsoc.org/wiki/Thunderstorm

Doswell, C. A., Brooks, H. E., & Maddox, R. A. (1996). Flash flood fore- casting: An ingredients-based methodology. Weather and Forecasting, 11(4), 560–581. https://doi.org/10.1175/1520-0434(1996)011⟨0560: FFFAIB⟩2.0.CO;2

García-Ortega, E., Fita, L., Romero, R., López, L., Ramis, C., & Sánchez, J. L. (2007). Numerical simulation and sensitivity study of a severe hailstorm in northeast spain. Atmospheric Research, 83 (2), 225–241. https://doi. org/10.1016/j.atmosres.2005.08.004

MeteoSwiss. (2018, December 21). 2018: Rekordwärme und massive regenarmut. Retrieved March 23, 2020, from https://www.meteoschweiz. admin.ch/home/suche.subpage.html/de/data/blogs/2018/12/ 2018-waermstes-jahr-seit-messbeginn-.html?query=sommer+2018+ schweiz&topic=0&pageIndex=0&tab=search tab

Nisi, L., Martius, O., Hering, A., Kunz, M., & Germann, U. (2016). Spatial and temporal distribution of hailstorms in the alpine region: A long-term, high resolution, radar-based analysis. Quarterly Journal of the Royal Meteorological Society, 142(697), 1590–1604. https://doi.org/10.1002/ qj.2771

Powers, J. G., Klemp, J. B., Skamarock, W. C., Davis, C. A., Dudhia, J., Gill, D. O., Coen, J. L., Gochis, D. J., Ahmadov, R., Peckham, S. E., Grell, G. A., Michalakes, J., Trahan, S., Benjamin, S. G., Alexander, C. R., Dimego, G. J., Wang, W., Schwartz, C. S., Romine, G. S., . . . Duda, M. G. (2017). The weather research and forecasting model: Overview, system efforts, and future directions. Bulletin of the American Meteo- rological Society, 98(8), 1717–1737. https://doi.org/10.1175/BAMS-D- 15-00308.1

Sprenger, M., & Wernli, H. (2015). The LAGRANTO lagrangian analysis tool – version 2.0. Geoscientific Model Development, 8(8), 2569–2586.

Trefalt, S., Martynov, A., Barras, H., Besic, N., Hering, A. M., Lenggenhager, S., Noti, P., Röthlisberger, M., Schemm, S., Germann, U., & Martius, O. (2018). A severe hail storm in complex topography in switzerland - observations and processes. Atmospheric Research, 209, 76–94. https: //doi.org/10.1016/j.atmosres.2018.03.007

Wallace, J. M., & Hobbs, P. V. (2006, March 24). Atmospheric science: An introductory survey [Google-Books-ID: HZ2wNtDOU0oC]. Elsevier.
