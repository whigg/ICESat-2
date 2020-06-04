Written by Eduard Heijkoop, University of Colorado
Eduard.Heijkoop@colorado.edu

Download subsetted ICESat-2 data and mask data based on a given shapefile.
Output is a csv file with longitude (deg), latitude (deg), height (m above WGS84), UTC time
Two options:
1. GCP: ATL03 high confidence photons over land (e.g. for ground control points)
2. OCEAN: ATL03 high & medium confidence photons over water

Requirements:
- Run Get_Token_NSIDC.py (with your NASA EarthData username) to generate a valid token. Tokens remain valid for 30 days, after which you need to run this script again.
- Shapefile to mask photons; (1) will keep everything inside it, (2) will keep everything outside it. Example is the OpenStreetMap dataset (https://osmdata.openstreetmap.de/data/land-polygons.html)
- Input file with: name, lon_min, lon_max, lat_min, lat_max. Set "header=None" for no header and "header=0" for a single header line in this file
- Change files to the right username, files and folders
- For OCEAN a geotiff of the DTU18 Mean Sea Surface is needed, available here: https://drive.google.com/file/d/1Jog6s0kfQ9ipFjX3srRo4mXVn_-TMyA6/view?usp=sharing

Run either one of the main scripts and ICESat-2 data will be downloaded into the folder specified by the name in your input file. The .h5 files will be loaded into memory, processed and masked, then written to a file.
