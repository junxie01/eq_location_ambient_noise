# eq_location_ambient_noise
# citations:
1, Comparison of ground truth location of earthquake from InSAR and from ambient seismic noise: A case study of the 1998 Zhangbei earthquake, J Xie, X Zeng, W Chen, Z Zhan, 2011 Earthquake Science.
2, Ground truth location of earthquakes by use of ambient seismic noise from a sparse seismic network: A case study in Western Australia, X Zeng, J Xie, S Ni, Pure and Applied Geophysics 172 (6), 1397-1407, 2015

# this code is used to locate the earthquake by grid search, 
# taking the nearest station (~dozens km) as reference
# please prepare the group velocity dispersion curve from earthquake data and NCF between the reference station and remote stations 
# the code will read the dspersion curve file names by default. Please check the source code.
# unlike these two papers above, the original time of the earthquake is estimated by least square method. 
# 

# how to compile:
  make 
# how to use
  ./locat loc_par
# parameter file:
  loc_par
# which has the form as:
loc_station.lst    # station list, with: name, stlo, stla
NUUG               # reference station name
-52.3440 71.6400   # first guess earthquake location, usually given by default in extracting group velocity, such as USGS/IRIS.
output             # output file name
-53 -51 0.01       # lantitude range and step
71.2  72  0.01     # longtitude range and step
