def strip_gdalinfo_lonlat(input_str):
	input_str = input_str.split("(")
	input_str = input_str[len(input_str)-1]
	input_str = input_str.rstrip(")\n")
	input_str = input_str.split(",")
	lon = input_str[0]
	lat = input_str[1]
	lon = lon.lstrip(" ")
	lat = lat.lstrip(" ")
	
	lon_degrees = lon.split('d')
	lon_remain = lon_degrees[1]
	lon_degrees = lon_degrees[0]
	lon_minutes = lon_remain.split("'")
	lon_remain = lon_minutes[1]
	lon_minutes = lon_minutes[0]
	lon_seconds = lon_remain.split('"')
	lon_hemisphere = lon_seconds[1]
	lon_seconds = lon_seconds[0]
	lon = float(lon_degrees) + float(lon_minutes)/60 + float(lon_seconds)/3600
	if lon_hemisphere == "W":
		lon = -1*lon
	

	lat_degrees = lat.split('d')
	lat_remain = lat_degrees[1]
	lat_degrees = lat_degrees[0]
	lat_minutes = lat_remain.split("'")
	lat_remain = lat_minutes[1]
	lat_minutes = lat_minutes[0]
	lat_seconds = lat_remain.split('"')
	lat_hemisphere = lat_seconds[1]
	lat_seconds = lat_seconds[0]    
	lat = float(lat_degrees) + float(lat_minutes)/60 + float(lat_seconds)/3600
	if lat_hemisphere == "S":
		lat = -1*lat
	
	return lon, lat

def cat_str_API(beam):
	beam_command = '/gt'+beam+'/heights/h_ph,/gt'+beam+'/heights/lon_ph,/gt'+beam+'/heights/lat_ph,/gt'+beam+'/heights/delta_time,' \
		'/gt'+beam+'/heights/signal_conf_ph,/gt'+beam+'/geolocation/sigma_h,' \
		'/gt'+beam+'/geolocation/reference_photon_lon,/gt'+beam+'/geolocation/reference_photon_lat,/gt'+beam+'/geophys_corr/tide_ocean,/gt'+beam+'/geophys_corr/dac,' \
		'/gt'+beam+'/geolocation/delta_time,/gt'+beam+'/geophys_corr/delta_time,'
	return beam_command

def find_extents_shp(shp):
        import geopandas as gpd
        import json
        import numpy as np

        df_shp = gpd.read_file(shp)
        g_shp = json.loads(df_shp.to_json())
        lon = np.empty([0,1],dtype=float)
        lat = np.empty([0,1],dtype=float)

        for i in range(len(g_shp['features'])):
                tmp = g_shp['features'][i]['geometry']['coordinates']
                tmp2 = np.array(tmp[0])
                tmp_lon = tmp2[:,0]
                tmp_lat = tmp2[:,1]
                lon = np.append(lon,tmp_lon)
                lat = np.append(lat,tmp_lat)

        lon_min = np.min(lon)
        lon_max = np.max(lon)
        lat_min = np.min(lat)
        lat_max = np.max(lat)

        return lon_min, lon_max, lat_min, lat_max

def SRTM_filter(lon,lat,h):
        import numpy as np
        import os, sys
        from osgeo import gdal, gdalconst
        import getpass

        user = 'EHeijkoop'
        pw = getpass.getpass()
        
        lon_min = np.min(lon)
        lon_max = np.max(lon)
        lat_min = np.min(lat)
        lat_max = np.max(lat)



        SRTM_list = []
        lon_range = range(int(np.floor(lon_min)),int(np.floor(lon_max))+1)
        lat_range = range(int(np.floor(lat_min)),int(np.floor(lat_max))+1)

        for jj in range(len(lon_range)):
                for kk in range(len(lat_range)):
                        if lon_range[jj] >= 0:
                                lonLetter = 'E'
                        else:
                                lonLetter = 'W'
                        if lat_range[kk] >= 0:
                                latLetter = 'N'
                        else:
                                latLetter = 'S'
                        lonCode = f"{int(np.abs(np.floor(lon_range[jj]))):03d}"
                        latCode = f"{int(np.abs(np.floor(lat_range[kk]))):02d}"
                        SRTMID = latLetter + latCode + lonLetter + lonCode
                        SRTM_list.append(SRTMID)
        
        for jj in range(len(SRTM_list)):
                DL_command = 'wget --user=' + user + ' --password=' + pw + ' https://e4ftl01.cr.usgs.gov//MODV6_Dal_D/SRTM/SRTMGL1.003/2000.02.11/' + SRTM_list[jj] + '.SRTMGL1.hgt.zip'
                os.system(DL_command)
                exists = os.path.isfile(SRTM_list[jj] + '.SRTMGL1.hgt.zip')
                if exists:
                        unzip_command = 'unzip ' + SRTM_list[jj] + '.SRTMGL1.hgt.zip'
                        os.system(unzip_command)
                        delete_command = 'rm ' + SRTM_list[jj] + '.SRTMGL1.hgt.zip'
                        os.system(delete_command)

                        merge_command = merge_command + SRTM_list[jj] + '.hgt '

def landmask(lon,lat,h,shp,inverse_flag):
        #Old landmask using GMT Select
        import numpy as np
        import os, sys
        import pandas as pd

        shp_no_ext = shp.split('.shp')
        shp_no_ext = shp_no_ext[0]
        shp_gmt = shp_no_ext + '.gmt'

        f_tmp = open('tmp_lonlat.txt','w')
        np.savetxt(f_tmp,np.c_[lon,lat,h],fmt='%10.5f',delimiter=',')
        f_tmp.close()
        
        os.system('ogr2ogr -f "GMT" '+ shp_gmt + ' ' + shp + ' -mapFieldType Integer64=Real')
        if inverse_flag == 0: #in the shapefile, i.e. land
                os.system('gmt select tmp_lonlat.txt -F' + shp_gmt + ' > tmp_landmask.txt')
        elif inverse_flag == 1: #out of the shapefile, i.e. water
                os.system('gmt select tmp_lonlat.txt -If -F' + shp_gmt + ' > tmp_landmask.txt')
        
        df_tmp = pd.read_csv('tmp_landmask.txt',header=None,names=['lon','lat','h'],dtype={'lon':'float','lat':'float','h':'float'})

        lon_masked = df_tmp.lon
        lat_masked = df_tmp.lat
        h_masked = df_tmp.h
        os.system('rm tmp_lonlat.txt')
        os.system('rm tmp_landmask.txt')
        return lon_masked,lat_masked,h_masked

def gps2utc(gps_time):
        #Converts ICESat-2 GPS time to UTC time
        from datetime import date,time,timedelta,datetime
        import numpy as np

        t0 = datetime(1980,1,6,0,0,0,0)
        leap_seconds = -18 #applicable to everything after 2017-01-01, UTC is currently 18 s behind GPS
        dt = (gps_time + leap_seconds) * timedelta(seconds=1)
        utc_time = t0+dt
        utc_time_str = [str(x) for x in utc_time]
        return utc_time_str


def inpoly(lon,lat,shp_path):
    #Derived from Darren Engwirda's inpoly Ray casting algorithm
    #https://github.com/dengwirda/inpoly
    import numpy as np
    import geopandas as gpd

    ftol = 1e-15

    shp = gpd.read_file(shp_path)

    lon_coast = np.empty([0,1],dtype=float)
    lat_coast = np.empty([0,1],dtype=float)

    for ii in range(len(shp)):
        tmp_geom_type = shp.geometry[ii].geom_type
        if tmp_geom_type == 'Polygon':
            tmp = np.asarray(shp.geometry[ii].exterior.xy)
            lon_coast = np.append(lon_coast,tmp[0,:])
            lon_coast = np.append(lon_coast,np.nan)
            lat_coast = np.append(lat_coast,tmp[1,:])
            lat_coast = np.append(lat_coast,np.nan)
        elif tmp_geom_type == 'MultiPolygon':
            tmp_list = list(shp.geometry[ii])
            for jj in range(len(tmp_list)):
                tmp = np.asarray(tmp_list[jj].exterior.xy)
                lon_coast = np.append(lon_coast,tmp[0,:])
                lon_coast = np.append(lon_coast,np.nan)
                lat_coast = np.append(lat_coast,tmp[1,:])
                lat_coast = np.append(lat_coast,np.nan)
    
    #input lon/lat points
    vert = np.stack((lon,lat),axis=1)
    nvrt = vert.shape[0]
    #vertices of shapefile lon/lat
    node = np.stack((lon_coast,lat_coast),axis=1)
    nnod = node.shape[0]
    #indices of shapefile vertices, just ascending order
    edge = np.stack((np.array(range(0,nnod-1)),np.array(range(1,nnod))),axis=1)
    edge = np.concatenate((edge,[[nnod-1,0]]),axis=0)

    #flip to ensure the y-axis is the "long" axis
    vmin = np.min(vert,axis=0)
    vmax = np.max(vert,axis=0)
    ddxy = vmax - vmin
    lbar = np.sum(ddxy) / 2

    if ddxy[0] > ddxy[1]:
        vert = np.fliplr(vert)
        node = np.fliplr(node)
    
    #sort points via y-value
    with np.errstate(invalid='ignore'):
        swap = node[edge[:,1],1] < node[edge[:,0],1]
    edge[swap,:] = np.fliplr(edge[swap,:])
    idx_sort = np.argsort(vert[:,1])
    vert = vert[idx_sort]

    feps = ftol * lbar**2
    veps = ftol * lbar

    nvrt = vert.shape[0]
    stat = np.zeros([nvrt],dtype=bool)
    bnds = np.zeros([nvrt],dtype=bool)

    #loop over polygon edges
    for epos in range(edge.shape[0]):
        inod = edge[epos,0]
        jnod = edge[epos,1]

        yone = node[inod,1]
        ytwo = node[jnod,1]
        xone = node[inod,0]
        xtwo = node[jnod,0]

        xmin = np.minimum(xone,xtwo)
        xmax = np.maximum(xone,xtwo)

        xmax = xmax + veps

        ymin = yone - veps
        ymax = ytwo + veps

        ydel = ytwo - yone
        xdel = xtwo - xone

        #Binary search
        ilow = 0
        iupp = nvrt

        while (ilow < iupp - 1): #-2 instead of -1?
            imid = np.int(ilow + np.floor((iupp - ilow)/2))
            if (vert[imid,1] < ymin):
                ilow = imid
            else:
                iupp = imid
        
        if (vert[ilow,1] >= ymin):
            ilow = ilow - 1
        
        for jpos in range(ilow+1,nvrt):
            if bnds[jpos]:
                continue

            xpos = vert[jpos,0]
            ypos = vert[jpos,1]

            if (ypos <= ymax):
                if (xpos >= xmin):
                    if (xpos <= xmax):
                        mul1 = ydel * (xpos - xone)
                        mul2 = xdel * (ypos - yone)

                        if (feps >= np.abs(mul2 - mul1)):
                            bnds[jpos] = True
                            stat[jpos] = True
                        elif (ypos == yone and  xpos == xone ):
                            bnds[jpos] = True
                            stat[jpos] = True
                        elif (ypos == ytwo and  xpos == xtwo ):
                            bnds[jpos] = True
                            stat[jpos] = True
                        elif (mul1 < mul2):
                            if (ypos >= yone and ypos < ytwo):
                                stat[jpos] = not stat[jpos]

                else:
                    if (ypos >= yone and ypos < ytwo):
                        stat[jpos] = not stat[jpos]

            else:
                break


    landmask = np.empty(stat.shape)
    landmask[idx_sort] = stat
    #returns landmask True/False array
    return landmask
