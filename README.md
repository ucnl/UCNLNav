# UCNLNav

## WARNING!

If this repository achieves:

25 stars - we'll add partial implementation in MATLAB  
50 stars - we'll add AoA (Angle of arrival) estimation routines  
100 stars - we'll add C implementation of the library  
150 stars - we'll add complete documentation for the library  

This library contains routines for:
* Solving geodetic problems (Vincenty equations, haversine)
* Solving navigation & positioning problems (TOA & TDOA)

Basic example of usage is in UCNLNav_Tests demo application.


30-NOV-2019 Update:
* Added partial implementation of the library in Rust:
- All the functionality from Algorithms.cs (TOA/TDOA not covered with tests yet)


22-NOV-2019 Update:  
* In GeoPoints.cs new classes for metric point (MPoint and MPoint3D)  
* In Navigation.cs new methods for calculating centroids of clouds of MPoint and MPoint3D, 
converting between local and geographic coordinate systems, methods for calculating statistics (CEP, SEP, STD, MRSE, DRMS)


26-AUG-2019 Update:  
* Routines for VLBL (Virtual long baseline) positioning
* TDOA solution in 3D


### Please, let us know that our work is useful for you: star this repository =)
