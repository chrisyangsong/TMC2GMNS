import os
import math
import datetime
import numpy as np
import pandas as pd

import utm
import time
import os.path

import shapely
import shapely.geometry as geom
import geopandas as gpd

from shapely.geometry import LineString
from shapely.geometry import MultiLineString
from shapely.geometry import Point
from shapely.wkt import loads

def getDegree(latA, lonA, latB, lonB):
    """
    Args:
        point p1(latA, lonA)
        point p2(latB, lonB)
    Returns:
        bearing between the two GPS points,
        default: the basis of heading direction is north
    """
    radLatA = math.radians(latA)
    radLonA = math.radians(lonA)
    radLatB = math.radians(latB)
    radLonB = math.radians(lonB)
    dLon = radLonB - radLonA
    y = math.sin(dLon) * math.cos(radLatB)
    x = math.cos(radLatA) * math.sin(radLatB) - math.sin(radLatA) * math.cos(radLatB) * math.cos(dLon)
    brng = np.degrees(math.atan2(y, x))
    brng = (brng + 360) % 360
    # brng = 360 - brng
    return brng

def bearing(geometry):
   x, y = shapely.wkt.loads(geometry).coords.xy
   xy = pd.DataFrame({'LON':x,'LAT':y})
   return getDegree(xy.head(1)['LAT'], xy.head(1)['LON'], xy.tail(1)['LAT'], xy.tail(1)['LON'])



def MatchTMC2BASENetwork(link_tmc,link_base,link_measurement_TMC):
    link_base = pd.read_csv(link_base, low_memory=False)
    link_base['bearing_angle'] = link_base['geometry'].apply(bearing)
    
    multiline_string_base_list = []
    multiline_string_base_list_sub = []
    for j in link_base.index:
        link_base_geometry_list = link_base.loc[j,'geometry'][12:-1].split(", ")
        for link_base_geometry in link_base_geometry_list:
            multiline_string_base_list_sub.append((float(link_base_geometry.split(" ")[0]),float(link_base_geometry.split(" ")[1])))
        multiline_string_base_list_sub = tuple(multiline_string_base_list_sub)
        multiline_string_base_list.append(multiline_string_base_list_sub)
        multiline_string_base_list_sub = []

    line_base = MultiLineString(multiline_string_base_list)



    link_measurement_TMC = pd.read_csv(link_measurement_TMC, low_memory=False)
    link_TMC = pd.read_csv(link_tmc, low_memory=False)
    link_TMC = link_TMC[link_TMC['link_id'].isin(set(link_measurement_TMC['link_id']))]
    link_TMC = link_TMC.reset_index()
    link_TMC = link_TMC.drop(['index'], 1)

    in_bbox_index_list = []
    for i in link_TMC.index:
        link_tmc_geometry_list = link_TMC.loc[i,'geometry'][12:-1].split(", ")
        start_longitude = float(link_tmc_geometry_list[0].split(" ")[0])
        start_latitude = float(link_tmc_geometry_list[0].split(" ")[1])
        end_longitude = float(link_tmc_geometry_list[1].split(" ")[0])
        end_latitude = float(link_tmc_geometry_list[1].split(" ")[1])
        if (start_longitude > line_base.bounds[0]) & (start_longitude < line_base.bounds[2]) & \
            (end_longitude > line_base.bounds[0]) & (end_longitude < line_base.bounds[2]) & \
                (start_latitude > line_base.bounds[1]) & (start_latitude < line_base.bounds[3]) & \
            (end_latitude > line_base.bounds[1]) & (end_latitude < line_base.bounds[3]):
            in_bbox_index_list.append(i)

    link_TMC = link_TMC.loc[in_bbox_index_list]
    link_TMC = link_TMC.reset_index()
    link_TMC = link_TMC.drop(['index'], 1)

    link_TMC['bearing_angle'] = link_TMC['geometry'].apply(bearing)



    matching_tmc2base_dict = {}
    k = 0
    p = 1
    for i in link_TMC.index:
        angle_list = list(abs(link_TMC.loc[i,'bearing_angle']-list(link_base['bearing_angle'])))
        small_angle_list = [i for i, value in enumerate(angle_list) if value < 45]
        for dividing_rate in range(1,11):
            start_point_tmc = loads(link_TMC['geometry'][i]).interpolate((dividing_rate-1)/10, normalized=True)
            end_point_tmc = loads(link_TMC['geometry'][i]).interpolate(dividing_rate/10, normalized=True)
            distance_list = []
            
            for j in small_angle_list:
                distance_list.append(LineString([start_point_tmc, end_point_tmc]).distance(loads(link_base['geometry'][j])))
            df_distance = pd.DataFrame({'index':small_angle_list,'distance':distance_list})
            nearest_index = int(df_distance.loc[df_distance['distance'].idxmin()]['index'])

            matching_tmc2base_dict[k] = {'name_tmc':link_TMC.loc[i]['link_id'],\
                                            'from_node_id_tmc':link_TMC.loc[i]['from_node_id'],\
                                            'to_node_id_tmc':link_TMC.loc[i]['to_node_id'],\
                                            'bearing_angle_tmc':link_TMC.loc[i]['bearing_angle'],\
                                            'geometry_tmc':link_TMC.loc[i]['geometry'],\
                                            'name_base':link_base['name'][nearest_index],\
                                            'link_id_base':link_base['link_id'][nearest_index],\
                                            'from_node_id_base':link_base['from_node_id'][nearest_index],\
                                            'to_node_id_base':link_base['to_node_id'][nearest_index],\
                                            'geometry_base':link_base['geometry'][nearest_index],\
                                            'distance':min(distance_list),\
                                            'geometry_tmc_base':'MULTILINESTRING ('+ link_TMC.loc[i]['geometry'][11:] + \
                                                                ', ' + link_base['geometry'][nearest_index][11:]+')'}
            k += 1

            
        if link_TMC.index.get_loc(i) > p/10 * len(link_TMC): 
            print(str(p*10)+"%"+' matching completed!')
            p = p + 1
                

    matching_tmc2base = pd.DataFrame(matching_tmc2base_dict).transpose()
    matching_tmc2base_drop_duplicates = matching_tmc2base.drop(columns=['distance']).drop_duplicates()
    matching_tmc2base_drop_duplicates['distance'] = matching_tmc2base.loc[matching_tmc2base_drop_duplicates.index]['distance']
    matching_tmc2base_drop_duplicates = matching_tmc2base_drop_duplicates.reset_index()
    matching_tmc2base_drop_duplicates = matching_tmc2base_drop_duplicates.drop(['index'], 1)
    matching_tmc2base = matching_tmc2base_drop_duplicates

    matching_tmc2base.to_csv('matching_tmc2base.csv',index = False)
    print('matching_tmc2base.csv generated!')
 