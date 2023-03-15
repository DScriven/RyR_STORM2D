# RyR_STORM2D
Blink and Cluster Analyzer for surface ryanodine receptors.

This program analyzes dSTORM blinks emitted by fluorophores attached to ryanodine 
receptors on the surface of cardiomyocytes. The program uses Delaunay triangulation 
to divide the blinks into clusters according to the maximum distance that Ca2+ diffusion 
is assumed to excite other receptors, a parameter we called the neighbourhood limit. 

The program uses a two-step approach to define the clusters and how they are grouped 
according to the neighbourhood limit.

To compile this program you would need Qt5 (https://www.qt.io/ - Open Source ver 5.15.2) 
and access to CGAL (https://www.cgal.org/  ver 5.5 or later) and libtiff ver 4 or later 
(http://www.simplesystems.org/libtiff/). OpenGL ver 4.3 or later is also required.
 
To install this program on your Windows computer download RyR_STORM2D_Installer.exe. 
This program works on both Windows 10 & 11. The documents, installed as 
a subdirectory of where the program is installed, have a full description of its functions.
These documents can also be accessed from within the program using the F1 key.
       
