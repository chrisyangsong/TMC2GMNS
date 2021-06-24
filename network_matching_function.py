import os
import math
import datetime
import numpy as np
import pandas as pd
import time
import os.path


def LLs2Dist(lon1, lat1, lon2, lat2): #WGS84 transfer coordinate system to distance(meter) #xy #credit to xtHuang0927
    R = 6371
    dLat = (lat2 - lat1) * math.pi / 180.0
    dLon = (lon2 - lon1) * math.pi / 180.0

    a = math.sin(dLat / 2) * math.sin(dLat/2) + math.cos(lat1 * math.pi / 180.0) * math.cos(lat2 * math.pi / 180.0) * math.sin(dLon/2) * math.sin(dLon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    distance = R * c * 1000
    return distance


def osm2tmc_matching(osm_network,TMC_Identification):
    "Get the OSM Network"
    network_type = 'auto'
    import osm2gmns as og

    net = og.getNetFromOSMFile(osm_network,network_type=(network_type), default_lanes=True, default_speed=True)


    og.consolidateComplexIntersections(net)
    og.outputNetToCSV(net, output_folder=os.getcwd())
    link_road = pd.read_csv('link.csv', low_memory=False)

    "Convert TMC Data into GMNS Format"

    tmc = pd.read_csv(TMC_Identification)

    
    
    '''build node.csv'''
    print('converting tmc data into gmns format...')


    node_tmc = pd.DataFrame()
    node_tmc['name'] = None
    node_tmc['x_coord'] = None
    node_tmc['y_coord'] = None
    node_tmc['z_coord'] = None
    node_tmc['node_type'] = None
    node_tmc['ctrl_type'] = None
    node_tmc['zone_id'] = None
    node_tmc['parent_node_id'] = None
    node_tmc['geometry'] = None

    for i in range(0,len(tmc)-1):
        if tmc.loc[i+1,'road_order'] > tmc.loc[i,'road_order']:
            node_tmc = node_tmc.append({'name': tmc.loc[i,'tmc'],\
                                        'x_coord': tmc.loc[i,'start_longitude'], \
                                        'y_coord': tmc.loc[i,'start_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_start',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'parent_node_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i,'start_longitude'].astype(str) + " " + tmc.loc[i,'start_latitude'].astype(str) +")"}, ignore_index=True)
        else:
            node_tmc = node_tmc.append({'name': tmc.loc[i,'tmc'],\
                                        'x_coord': tmc.loc[i,'start_longitude'], \
                                        'y_coord': tmc.loc[i,'start_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_start',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'parent_node_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i,'start_longitude'].astype(str) + " " + tmc.loc[i,'start_latitude'].astype(str) +")"}, ignore_index=True)
            node_tmc = node_tmc.append({'name': tmc.loc[i,'tmc']+'END',\
                                        'x_coord': tmc.loc[i,'end_longitude'], \
                                        'y_coord': tmc.loc[i,'end_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_end',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'parent_node_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i,'end_longitude'].astype(str) + " " + tmc.loc[i,'end_latitude'].astype(str) +")"}, ignore_index=True)


    node_tmc = node_tmc.append({'name': tmc.loc[i+1,'tmc'],\
                                        'x_coord': tmc.loc[i+1,'start_longitude'], \
                                        'y_coord': tmc.loc[i+1,'start_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_start',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'parent_node_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i+1,'start_longitude'].astype(str) + " " + tmc.loc[i+1,'start_latitude'].astype(str) +")"}, ignore_index=True)

    node_tmc = node_tmc.append({'name': tmc.loc[i+1,'tmc']+'END',\
                                        'x_coord': tmc.loc[i+1,'end_longitude'], \
                                        'y_coord': tmc.loc[i+1,'end_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_end',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'parent_node_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i+1,'end_longitude'].astype(str) + " " + tmc.loc[i+1,'end_latitude'].astype(str) +")"}, ignore_index=True)

    node_tmc.index.name = 'node_id'

    node_tmc.index += 100000001 #index from 0

    node_tmc.to_csv('node_tmc.csv')
    print('node_tmc.csv generated!')


    '''build link_tmc.csv'''
    link_tmc = pd.DataFrame()
    link_tmc['name'] = None
    link_tmc['corridor_id'] = None
    link_tmc['corridor_link_order'] = None
    link_tmc['from_node_id'] = None
    link_tmc['to_node_id'] = None
    link_tmc['directed'] = None
    link_tmc['geometry_id'] = None
    link_tmc['geometry'] = None
    link_tmc['dir_flag'] = None
    link_tmc['parent_link_id'] = None
    link_tmc['length'] = None
    link_tmc['grade'] = None
    link_tmc['facility_type'] = None
    link_tmc['capacity'] = None
    link_tmc['free_speed'] = None
    link_tmc['lanes'] = None

    for i in range(0,len(tmc)):
        link_tmc = link_tmc.append({'name': tmc.loc[i,'tmc'],\
                                    'corridor_id': tmc.loc[i,'road']+'_'+tmc.loc[i,'direction'],\
                                    'corridor_link_order' : tmc.loc[i,'road_order'],\
                                    'from_node_id': node_tmc[(node_tmc['x_coord']==tmc.loc[i,'start_longitude']) & (node_tmc['y_coord']==tmc.loc[i,'start_latitude'])].index.values[0], \
                                    'to_node_id': node_tmc[(node_tmc['x_coord']==tmc.loc[i,'end_longitude']) & (node_tmc['y_coord']==tmc.loc[i,'end_latitude'])].index.values[0],\
                                    'directed': 1,\
                                    'geometry_id': None,\
                                    'geometry': "LINESTRING (" + tmc.loc[i,'start_longitude'].astype(str) + " " + tmc.loc[i,'start_latitude'].astype(str) + "," +\
                                        tmc.loc[i,'end_longitude'].astype(str) +" "+ tmc.loc[i,'end_latitude'].astype(str) + ")",\
                                    'dir_flag': 1,\
                                    'parent_link_id': None,\
                                    'length': tmc.loc[i,'miles'],\
                                    'grade': None,\
                                    'facility_type': 'interstate' if tmc.loc[i,'road'][0] == 'I'else None ,\
                                    'capacity':None,\
                                    'free_speed':None,\
                                    'lanes': None}, ignore_index=True)
    link_tmc.index.name = 'link_id'
    link_tmc.index += 100000001


    link_tmc.to_csv('link_tmc.csv')
    print('link_tmc.csv generated!')


    for i in range(len(link_road)):
        lon_list = []
        lat_list = [] 
        link_geometry_list = link_road.loc[i,'geometry'][12:-1].split(", ")
        for link_geometry in link_geometry_list:
            lon_list.append(float(link_geometry.split(" ")[0]))
            lat_list.append(float(link_geometry.split(" ")[1]))
        center_lon = np.mean(lon_list)
        center_lat = np.mean(lat_list)


        distance_list = []
        for j in link_tmc.index:
            lon_tmc_list = []
            lat_tmc_list = []
            link_tmc_geometry_list = link_tmc.loc[j,'geometry'][12:-1].split(",")
            for link_tmc_geometry in link_tmc_geometry_list:
                lon_tmc_list.append(float(link_tmc_geometry.split(" ")[0]))
                lat_tmc_list.append(float(link_tmc_geometry.split(" ")[1]))
            center_tmc_lon = np.mean(lon_tmc_list)
            center_tmc_lat = np.mean(lat_tmc_list)
            
            distance_list.append(LLs2Dist(center_lon, center_lat, center_tmc_lon, center_tmc_lat))
        nearest_index = distance_list.index(min(distance_list))

        link_road.loc[i,'tmc_name'] = link_tmc.iloc[nearest_index]['name']
        link_road.loc[i,'tmc_link_id'] = int(link_tmc.iloc[[nearest_index]].index.values[0])
        link_road.loc[i,'matched_from_node_id'] = int(link_tmc.iloc[nearest_index]['from_node_id'])
        link_road.loc[i,'matched_to_node_id'] = int(link_tmc.iloc[nearest_index]['to_node_id'])
        if i % 5000 == 0 : print("{:.0%}".format(i/len(link_road))+' matching completed!')
    link_road.to_csv('network_matching.csv')
    print('osm2tmc_network_matching.csv generated!')



def osm2utd_matching(osm_network,utd_links,utd_detectors):
    "Get the OSM Network"
    network_type = 'auto'
    import osm2gmns as og

    net = og.getNetFromOSMFile(osm_network,network_type=(network_type), default_lanes=True, default_speed=True)


    og.consolidateComplexIntersections(net)
    og.outputNetToCSV(net, output_folder=os.getcwd())
    link_road = pd.read_csv('link.csv', low_memory=False)

    "Convert utd Data into GMNS Format"
    utd = pd.read_csv(utd_links)
    detectors = pd.read_csv(utd_detectors)

    utd.rename(columns = {'Unnamed: 0':'utd','long':'utd_long','lat':'utd_lat','citycode':'utd_citycode'}, inplace = True)
    detectors.rename(columns = {'Unnamed: 0':'detectors','long':'detectors_long',\
                    'lat':'detectors_lat','citycode':'detectors_citycode','length':'detectors_length'}, inplace = True)
    detectors['linkid'] = detectors['linkid'].astype(int)
    detectors['lanes'] = detectors['lanes'].astype(int)

    detectors = detectors[-detectors['linkid'].duplicated()]
    utd_detectors = utd.merge(detectors, on='linkid', how='left')  


    num_link = 0
    '''build node.csv'''
    print('converting utd data into gmns format...')


    node_utd = pd.DataFrame()
    node_utd['name'] = None
    node_utd['x_coord'] = None
    node_utd['y_coord'] = None
    node_utd['z_coord'] = None
    node_utd['node_type'] = None
    node_utd['ctrl_type'] = None
    node_utd['zone_id'] = None
    node_utd['parent_node_id'] = None
    node_utd['geometry'] = None

    for i in range(len(utd)-1):
        if utd.loc[i+1,'order'] > utd.loc[i,'order']:
            node_utd = node_utd.append({'name': utd.loc[i,'utd'],\
                                        'x_coord': utd.loc[i,'utd_long'], \
                                        'y_coord': utd.loc[i,'utd_lat'],\
                                        'z_coord': None,\
                                        'node_type': utd.loc[i,'linkid'],\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'parent_node_id': None,\
                                        'geometry': "POINT (" + utd.loc[i,'utd_long'].astype(str) + " " + utd.loc[i,'utd_lat'].astype(str) +")"}, ignore_index=True)
            num_link += 1
        else:
            
            node_utd = node_utd.append({'name': utd.loc[i,'utd'],\
                                        'x_coord': utd.loc[i,'utd_long'], \
                                        'y_coord': utd.loc[i,'utd_lat'],\
                                        'z_coord': None,\
                                        'node_type': utd.loc[i,'linkid'],\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'parent_node_id': None,\
                                        'geometry': "POINT (" + utd.loc[i,'utd_long'].astype(str) + " " + utd.loc[i,'utd_lat'].astype(str) +")"}, ignore_index=True)

    node_utd = node_utd.append({'name': utd.loc[i+1,'utd'],\
                                        'x_coord': utd.loc[i+1,'utd_long'], \
                                        'y_coord': utd.loc[i+1,'utd_lat'],\
                                        'z_coord': None,\
                                        'node_type': utd.loc[i+1,'linkid'],\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'parent_node_id': None,\
                                        'geometry': "POINT (" + utd.loc[i+1,'utd_long'].astype(str) + " " + utd.loc[i+1,'utd_lat'].astype(str) +")"}, ignore_index=True)


    node_utd.index.name = 'node_id'
    node_utd.index += 100000001 #index from 0



    node_utd.to_csv('node_utd.csv')
    print('node_utd.csv generated!')


    '''build link_utd.csv'''
    link_utd = pd.DataFrame()
    link_utd['name'] = None
    link_utd['corridor_id'] = None
    link_utd['corridor_link_order'] = None
    link_utd['from_node_id'] = None
    link_utd['to_node_id'] = None
    link_utd['directed'] = None
    link_utd['geometry_id'] = None
    link_utd['geometry'] = None
    link_utd['dir_flag'] = None
    link_utd['parent_link_id'] = None
    link_utd['length'] = None
    link_utd['grade'] = None
    link_utd['facility_type'] = None
    link_utd['capacity'] = None
    link_utd['free_speed'] = None
    link_utd['lanes'] = None
    link_utd['detid'] = None

    for i in range(0,len(utd_detectors)-1):
        if utd_detectors.loc[i+1,'order'] > utd_detectors.loc[i,'order']:
            link_utd = link_utd.append({'name': utd_detectors.loc[i,'utd'],\
                                        'corridor_id': utd_detectors.loc[i,'road'],\
                                        'corridor_link_order' : utd_detectors.loc[i,'order'],\
                                        'from_node_id': node_utd[(node_utd['x_coord']==utd_detectors.loc[i,'utd_long']) & \
                                                        (node_utd['y_coord']==utd_detectors.loc[i,'utd_lat'])].index.values[0], \
                                        'to_node_id': node_utd[(node_utd['x_coord']==utd_detectors.loc[i+1,'utd_long']) & \
                                                        (node_utd['y_coord']==utd_detectors.loc[i+1,'utd_lat'])].index.values[0],\
                                        'directed': None,\
                                        'geometry_id': None,\
                                        'geometry': "LINESTRING (" + utd_detectors.loc[i,'utd_long'].astype(str) + " " + \
                                                    utd_detectors.loc[i,'utd_lat'].astype(str) + "," + \
                                                    utd_detectors.loc[i+1,'utd_long'].astype(str) +" "+ \
                                                    utd_detectors.loc[i+1,'utd_lat'].astype(str) + ")",\
                                        'dir_flag': None,\
                                        'parent_link_id': None,\
                                        'length': None,\
                                        'grade': utd_detectors.loc[i,'fclass'],\
                                        'facility_type': None ,\
                                        'capacity':None,\
                                        'free_speed':utd_detectors.loc[i,'limit'],\
                                        'lanes': utd_detectors.loc[i,'lanes'],\
                                        'detid': utd_detectors.loc[i,'detid']}, ignore_index=True)
    link_utd.index.name = 'link_id'
    link_utd.index += 100000001


    link_utd.to_csv('link_utd.csv')
    print('link_utd.csv generated!')


    for i in range(len(link_road)):
        lon_list = []
        lat_list = [] 
        link_geometry_list = link_road.loc[i,'geometry'][12:-1].split(", ")
        for link_geometry in link_geometry_list:
            lon_list.append(float(link_geometry.split(" ")[0]))
            lat_list.append(float(link_geometry.split(" ")[1]))
        center_lon = np.mean(lon_list)
        center_lat = np.mean(lat_list)


        distance_list = []
        for j in link_utd.index:
            lon_utd_list = []
            lat_utd_list = []
            link_utd_geometry_list = link_utd.loc[j,'geometry'][12:-1].split(",")
            for link_utd_geometry in link_utd_geometry_list:
                lon_utd_list.append(float(link_utd_geometry.split(" ")[0]))
                lat_utd_list.append(float(link_utd_geometry.split(" ")[1]))
            center_utd_lon = np.mean(lon_utd_list)
            center_utd_lat = np.mean(lat_utd_list)
            
            distance_list.append(LLs2Dist(center_lon, center_lat, center_utd_lon, center_utd_lat))
        nearest_index = distance_list.index(min(distance_list))

        link_road.loc[i,'detid'] = link_utd.iloc[nearest_index]['detid']
        link_road.loc[i,'utd_link_id'] = int(link_utd.iloc[[nearest_index]].index.values[0])
        link_road.loc[i,'matched_from_node_id'] = int(link_utd.iloc[nearest_index]['from_node_id'])
        link_road.loc[i,'matched_to_node_id'] = int(link_utd.iloc[nearest_index]['to_node_id'])
        if i % 5000 == 0 : print("{:.0%}".format(i/len(link_road))+' matching completed!')
    link_road.to_csv('network_matching.csv')
    print('osm2utd_network_matching.csv generated!')



