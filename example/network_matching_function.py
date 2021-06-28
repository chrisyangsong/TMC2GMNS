import os
import math
import datetime
import numpy as np
import pandas as pd
import time
import os.path


def TMCIdentification2GMNSNodeLinkFiles(TMC_file):
    "Convert TMC Data into GMNS Format"
    tmc = pd.read_csv(TMC_file)
    
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


def ConvertTMCReading2Measurement(Reading,link_tmc):

    ## reading by detid
    reading = pd.read_csv(Reading)
    # reading = reading.loc[0:500]
    link_tmc = pd.read_csv(link_tmc)

    reading_dict = {}
    gp = reading.groupby('tmc_code')
    for key, form in gp:
        reading_dict[key] = {
            'measurement_tstamp':form['measurement_tstamp'].tolist(),
            'speed':form['speed'].tolist()
            }


    '''build measurement_tmc.csv''' 

    measurement_tmc = pd.DataFrame()
    measurement_tmc['link_id_tmc'] = None
    measurement_tmc['corridor_id'] = None
    measurement_tmc['corridor_link_order'] = None
    measurement_tmc['from_node_id'] = None
    measurement_tmc['to_node_id'] = None
    measurement_tmc['time_period'] = None
    measurement_tmc['date'] = None
    measurement_tmc['geometry'] = None
    measurement_tmc['volume'] = None
    measurement_tmc['travel_time'] = None
    measurement_tmc['speed'] = None
    measurement_tmc['reference_speed'] = None
    measurement_tmc['density'] = None
    measurement_tmc['queue'] = None
    measurement_tmc['notes'] = None


    k=0
    measurement_tmc_dict = {}
    for i in link_tmc.index:
        try:
            reading_dict_selected = reading_dict[link_tmc['name'][i]]
            for j in range(0,len(reading_dict_selected['measurement_tstamp']),3):
                measurement_tmc_dict[k] = {'link_id_tmc': i,\
                                                'corridor_id': link_tmc['corridor_id'][i],\
                                                'corridor_link_order' : link_tmc['corridor_link_order'][i],\
                                                'from_node_id': link_tmc.loc[i,'from_node_id'], \
                                                'to_node_id': link_tmc.loc[i,'to_node_id'], \
                                                'time_period': reading_dict_selected['measurement_tstamp'][j][11:13]+\
                                                    reading_dict_selected['measurement_tstamp'][j][14:16]+'_'+\
                                                    reading_dict_selected['measurement_tstamp'][j+3][11:13]+\
                                                    reading_dict_selected['measurement_tstamp'][j+3][14:16],\
                                                'date': reading_dict_selected['measurement_tstamp'][j][:10],\
                                                'geometry': link_tmc['geometry'][i],\
                                                'volume': None,\
                                                'travel_time': None,\
                                                'speed':round(np.mean(reading_dict_selected['speed'][j:j+3])),\
                                                'reference_speed': None,\
                                                'density': None,\
                                                'queue': None,\
                                                'notes': None }
                k += 1

        except:
            measurement_tmc_dict[k] = {'link_id_tmc': i,\
                                            'corridor_id': link_tmc['corridor_id'][i],\
                                            'corridor_link_order' : link_tmc['corridor_link_order'][i],\
                                            'from_node_id': link_tmc.loc[i,'from_node_id'], \
                                            'to_node_id': link_tmc.loc[i,'to_node_id'], \
                                            'time_period': None,\
                                            'date': None,\
                                            'geometry': link_tmc['geometry'][i],\
                                            'volume': None,\
                                            'travel_time': None,\
                                            'speed': None,\
                                            'reference_speed': None,\
                                            'density': None,\
                                            'queue': None,\
                                            'notes': None }
            k += 1


        # if i % 20 == 0 : print(i)
    measurement_tmc = pd.DataFrame(measurement_tmc_dict).transpose()
    measurement_tmc.to_csv('measurement_tmc.csv',index = False)
    print('measurement_tmc.csv generated!')

def LLs2Dist(lon1, lat1, lon2, lat2): #WGS84 transfer coordinate system to distance(meter) #xy #credit to xtHuang0927
    R = 6371
    dLat = (lat2 - lat1) * math.pi / 180.0
    dLon = (lon2 - lon1) * math.pi / 180.0

    a = math.sin(dLat / 2) * math.sin(dLat/2) + math.cos(lat1 * math.pi / 180.0) * math.cos(lat2 * math.pi / 180.0) * math.sin(dLon/2) * math.sin(dLon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    distance = R * c * 1000
    return distance


def MatchTMC2GMNSNetwork(link_tmc,link_gmns):
    link_tmc = pd.read_csv(link_tmc)
    link_gmns = pd.read_csv(link_gmns)
    matching_tmc2gmns_dict = {}
    k = 0
    for i in range(len(link_gmns)):
        lon_list = []
        lat_list = [] 
        link_geometry_list = link_gmns.loc[i,'geometry'][12:-1].split(", ")
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

        matching_tmc2gmns_dict[k] = {'name_osm':link_gmns['name'][i],\
                        'link_id_osm':link_gmns['link_id'][i],\
                        'from_node_id_osm':link_gmns['from_node_id'][i],\
                        'to_node_id_osm':link_gmns['to_node_id'][i],\
                        'geometry_osm':link_gmns['geometry'][i],\
                        'name_tmc':link_tmc.iloc[nearest_index]['name'],\
                        'link_id_tmc':int(link_tmc.iloc[[nearest_index]].index.values[0]),\
                        'from_node_id_tmc':int(link_tmc.iloc[nearest_index]['from_node_id']),\
                        'to_node_id_tmc':int(link_tmc.iloc[nearest_index]['to_node_id']),\
                        'geometry_tmc':link_tmc.iloc[nearest_index]['geometry']}
        k += 1
        if i % 5000 == 0 : print("{:.0%}".format(i/len(link_gmns))+' matching completed!')

    matching_tmc2gmns = pd.DataFrame(matching_tmc2gmns_dict).transpose()
    matching_tmc2gmns.to_csv('matching_tmc2gmns.csv',index = False)
    print('matching_tmc2gmns.csv generated!')


def ConvertMeasurementsBasedOnMatching(link_gmns,matching_tmc2gmns,measurement_tmc):
    link_gmns = pd.read_csv(link_gmns)
    matching_tmc2gmns = pd.read_csv(matching_tmc2gmns)
    measurement_tmc = pd.read_csv(measurement_tmc)
    '''build measurement_osm.csv''' 
    measurement_osm = pd.DataFrame()
    measurement_osm['link_id'] = None
    measurement_osm['osm_way_id'] = None
    measurement_osm['from_node_id'] = None
    measurement_osm['to_node_id'] = None
    measurement_osm['lanes'] = None
    measurement_osm['length'] = None
    measurement_osm['time_period'] = None
    measurement_osm['date'] = None
    measurement_osm['geometry'] = None
    measurement_osm['volume'] = None
    measurement_osm['speed'] = None


    k=0
    measurement_osm_dict = {}
    for i in link_gmns.index:
        try:
            measurement_tmc_selected = measurement_tmc[measurement_tmc['link_id_tmc'] == matching_tmc2gmns['link_id_tmc'][i]]
            for j in measurement_tmc_selected.index:
                measurement_osm_dict[k] = {'link_id': link_gmns['link_id'][i],\
                                                'osm_way_id':link_gmns['osm_way_id'][i],\
                                                'from_node_id': link_gmns['from_node_id'][i],\
                                                'to_node_id': link_gmns['to_node_id'][i],\
                                                'lanes': link_gmns['lanes'][i], \
                                                'length': link_gmns['length'][i], \
                                                'time_period': measurement_tmc_selected['time_period'][j],\
                                                'date': measurement_tmc_selected['date'][j],\
                                                'geometry': link_gmns['geometry'][i],\
                                                'volume': None,\
                                                'speed': measurement_tmc_selected['speed'][j]}

                k += 1
        except:
            measurement_osm_dict[k] = {'link_id': link_gmns['link_id'][i],\
                                                'osm_way_id':link_gmns['osm_way_id'][i],\
                                                'from_node_id': link_gmns['from_node_id'][i],\
                                                'to_node_id': link_gmns['to_node_id'][i],\
                                                'lanes': link_gmns['lanes'][i], \
                                                'length': link_gmns['length'][i], \
                                                'time_period':None,\
                                                'date': None,\
                                                'geometry': link_gmns['geometry'][i],\
                                                'volume': None,\
                                                'speed': None}

            k += 1
        # print(i)
        # if i % 5000 == 0 : print(i)
    measurement_osm = pd.DataFrame(measurement_osm_dict).transpose()
    measurement_osm.to_csv('measurement_osm.csv',index = False)
    print('measurement_osm.csv generated!')