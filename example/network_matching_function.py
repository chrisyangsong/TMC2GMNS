import os
import math
import datetime
import numpy as np
import pandas as pd
import time
import os.path


def TMCIdentification2GMNSNodeLinkFiles(TMC_file):
    #output:node_tmc,link_tmc
    "Convert TMC Data into GMNS Format"
    tmc = pd.read_csv(TMC_file)
    tmc = tmc.drop_duplicates(subset=['direction','road_order']).sort_values(by=['direction','road_order'])
    tmc = tmc.reset_index()
    tmc = tmc.drop(['index'], 1)
    
    '''build node.csv'''
    print('converting tmc data into gmns format...')
    p=1

    node_tmc = pd.DataFrame()
    node_tmc['name'] = None
    node_tmc['x_coord'] = None
    node_tmc['y_coord'] = None
    node_tmc['z_coord'] = None
    node_tmc['node_type'] = None
    node_tmc['ctrl_type'] = None
    node_tmc['zone_id'] = None
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
                                        'geometry': "POINT (" + tmc.loc[i,'start_longitude'].astype(str) + " " + tmc.loc[i,'start_latitude'].astype(str) +")"}, ignore_index=True)
        else:
            node_tmc = node_tmc.append({'name': tmc.loc[i,'tmc'],\
                                        'x_coord': tmc.loc[i,'start_longitude'], \
                                        'y_coord': tmc.loc[i,'start_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_start',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i,'start_longitude'].astype(str) + " " + tmc.loc[i,'start_latitude'].astype(str) +")"}, ignore_index=True)
            node_tmc = node_tmc.append({'name': tmc.loc[i,'tmc']+'END',\
                                        'x_coord': tmc.loc[i,'end_longitude'], \
                                        'y_coord': tmc.loc[i,'end_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_end',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i,'end_longitude'].astype(str) + " " + tmc.loc[i,'end_latitude'].astype(str) +")"}, ignore_index=True)

        if i > p/10 * len(tmc): 
            print(str(p*10)+"%"+' nodes completed!')
            p = p + 1

    node_tmc = node_tmc.append({'name': tmc.loc[i+1,'tmc'],\
                                        'x_coord': tmc.loc[i+1,'start_longitude'], \
                                        'y_coord': tmc.loc[i+1,'start_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_start',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i+1,'start_longitude'].astype(str) + " " + tmc.loc[i+1,'start_latitude'].astype(str) +")"}, ignore_index=True)

    node_tmc = node_tmc.append({'name': tmc.loc[i+1,'tmc']+'END',\
                                        'x_coord': tmc.loc[i+1,'end_longitude'], \
                                        'y_coord': tmc.loc[i+1,'end_latitude'],\
                                        'z_coord': None,\
                                        'node_type': 'tmc_end',\
                                        'ctrl_type': None,\
                                        'zone_id': None,\
                                        'geometry': "POINT (" + tmc.loc[i+1,'end_longitude'].astype(str) + " " + tmc.loc[i+1,'end_latitude'].astype(str) +")"}, ignore_index=True)

    node_tmc.index.name = 'node_id'

    node_tmc.index += 100000001 #index from 0

    node_tmc.to_csv('node_tmc.csv')
    print('node_tmc.csv (' + str(len(node_tmc)) + ' nodes' + ') generated!')

    '''build link_tmc.csv'''
    p = 1
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
                                    'length': tmc.loc[i,'miles'],\
                                    'grade': None,\
                                    'facility_type': 'interstate' if tmc.loc[i,'road'][0] == 'I'else None ,\
                                    'capacity':None,\
                                    'free_speed':None,\
                                    'lanes': None}, ignore_index=True)

        if i > p/10 * len(tmc): 
            print(str(p*10)+"%"+' links completed!')
            p = p + 1
            
    link_tmc.index.name = 'link_id'
    link_tmc.index += 100000001
    link_tmc.to_csv('link_tmc.csv')

    print('link_tmc.csv (' + str(len(link_tmc)) + ' links' + ') generated!')


def ConvertTMCReading2Measurement(Reading,link_tmc):

    link_tmc = pd.read_csv('link_tmc.csv')
    ## reading by detid
    reading = pd.read_csv('Reading.csv')
    reading = reading.loc[0:5000]


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
    p=1
    measurement_tmc_dict = {}
    for i in link_tmc.index:
        try:
            reading_dict_selected = reading_dict[link_tmc['name'][i]]
            for j in range(0,len(reading_dict_selected['measurement_tstamp'])):
                measurement_tmc_dict[k] = {'link_id_tmc': i,\
                                                'corridor_id': link_tmc['corridor_id'][i],\
                                                'corridor_link_order' : link_tmc['corridor_link_order'][i],\
                                                'from_node_id': link_tmc.loc[i,'from_node_id'], \
                                                'to_node_id': link_tmc.loc[i,'to_node_id'], \
                                                'time_period': reading_dict_selected['measurement_tstamp'][j][11:13]+\
                                                    reading_dict_selected['measurement_tstamp'][j][14:16]+'_'+\
                                                    reading_dict_selected['measurement_tstamp'][j+1][11:13]+\
                                                    reading_dict_selected['measurement_tstamp'][j+1][14:16],\
                                                'date': reading_dict_selected['measurement_tstamp'][j][:10],\
                                                'geometry': link_tmc['geometry'][i],\
                                                'volume': None,\
                                                'travel_time': None,\
                                                'speed':round(np.mean(reading_dict_selected['speed'][j:j+1])),\
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

        if link_tmc.index.get_loc(i) > p/10 * len(link_tmc): 
            print(str(p*10)+"%"+' measurement_tmc completed!')
            p = p + 1

    measurement_tmc = pd.DataFrame(measurement_tmc_dict).transpose()
    measurement_tmc = measurement_tmc.dropna(subset=['time_period']) #remove na at the end of day of unrecorded ones
    measurement_tmc = measurement_tmc.reset_index()
    measurement_tmc = measurement_tmc.drop(['index'], 1)
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
    


def MatchTMC2GMNSNetwork(link_tmc,link_base):
    link_tmc = pd.read_csv('link_tmc.csv')
    link_base = pd.read_csv('link.csv', low_memory=False)
    # link_base = link_base[link_base['link_type_name'].isin(['motorway','trunk','primary','secondary'])]
    link_base = link_base[link_base['link_type_name'].isin(['motorway','trunk'])]
    link_base = link_base.reset_index()
    link_base = link_base.drop(['index'], 1)
    matching_tmc2gmns_dict = {}
    k = 0
    p = 1

    for j in link_tmc.index:
        lon_tmc_list = []
        lat_tmc_list = []
        link_tmc_geometry_list = link_tmc.loc[j,'geometry'][12:-1].split(",")
        for link_tmc_geometry in link_tmc_geometry_list:
            lon_tmc_list.append(float(link_tmc_geometry.split(" ")[0]))
            lat_tmc_list.append(float(link_tmc_geometry.split(" ")[1]))
        center_tmc_lon = np.mean(lon_tmc_list)
        center_tmc_lat = np.mean(lat_tmc_list)
        tmc_lon_1 = lon_tmc_list[0]
        tmc_lon_2 = lon_tmc_list[-1]
        tmc_lat_1 = lat_tmc_list[0]
        tmc_lat_2 = lat_tmc_list[-1]
        if getDegree(tmc_lat_1,tmc_lon_1,tmc_lat_2,tmc_lon_2)>180:
            angle_tmc = getDegree(tmc_lat_2,tmc_lon_2,tmc_lat_1,tmc_lon_1)
        else:
            angle_tmc = getDegree(tmc_lat_1,tmc_lon_1,tmc_lat_2,tmc_lon_2)


        distance_list = []
        angle_list = []
        for i in range(len(link_base)):
            lon_list = []
            lat_list = [] 
            link_geometry_list = link_base.loc[i,'geometry'][12:-1].split(", ")
            for link_geometry in link_geometry_list:
                lon_list.append(float(link_geometry.split(" ")[0]))
                lat_list.append(float(link_geometry.split(" ")[1]))
            '''distance'''
            center_lon = np.mean(lon_list)
            center_lat = np.mean(lat_list)
            distance_list.append(LLs2Dist(center_lon, center_lat, center_tmc_lon, center_tmc_lat))
            '''angle '''
            base_lon_1 = lon_list[0]
            base_lon_2 = lon_list[-1]
            base_lat_1 = lat_list[0]
            base_lat_2 = lat_list[-1]
            if getDegree(base_lat_1,base_lon_1,base_lat_2,base_lon_2)>180:
                angle_base = getDegree(base_lat_2,base_lon_2,base_lat_1,base_lon_1)
            else:
                angle_base = getDegree(base_lat_1,base_lon_1,base_lat_2,base_lon_2)
            
            if abs(angle_tmc - angle_base) >= 90:
                relative_angle = 180 - abs(angle_tmc - angle_base)
            else:
                relative_angle = abs(angle_tmc - angle_base)
            angle_list.append(relative_angle)

        small_angle_list = [i for i, value in enumerate(angle_list) if value < 20]
        df_distance = pd.DataFrame({'distance':distance_list})
        nearest_index = df_distance.loc[small_angle_list].idxmin().values[0]

        matching_tmc2gmns_dict[k] = {'name_tmc':link_tmc.loc[j]['name'],\
                                    'link_id_tmc':link_tmc.loc[[j]].index.values[0],\
                                    'from_node_id_tmc':link_tmc.loc[j]['from_node_id'],\
                                    'to_node_id_tmc':link_tmc.loc[j]['to_node_id'],\
                                    'category_id_tmc':link_tmc.index.get_loc(j)+1,\
                                    'geometry_tmc':link_tmc.loc[j]['geometry'],\
                                    'name_base':link_base['name'][nearest_index],\
                                    'link_id_base':link_base['link_id'][nearest_index],\
                                    'from_node_id_base':link_base['from_node_id'][nearest_index],\
                                    'to_node_id_base':link_base['to_node_id'][nearest_index],\
                                    'category_id_base':link_tmc.index.get_loc(j)+1,\
                                    'geometry_base':link_base['geometry'][nearest_index],\
                                    'distance':min(distance_list)}
        k += 1

        
        if link_tmc.index.get_loc(j) > p/10 * len(link_tmc): 
            print(str(p*10)+"%"+' matching completed!')
            p = p + 1
            

    matching_tmc2gmns = pd.DataFrame(matching_tmc2gmns_dict).transpose()

    matching_tmc2gmns.to_csv('matching_tmc2gmns.csv',index = False)
    print('matching_tmc2gmns.csv generated!')


def ConvertMeasurementBasedOnMatching(link_base,matching_tmc2gmns,measurement_tmc):
    link_base = pd.read_csv('link.csv')
    link_base = link_base[link_base['link_type_name'].isin(['motorway','trunk','primary','secondary'])]
    link_base = link_base.reset_index()
    link_base = link_base.drop(['index'], 1)
    matching_tmc2gmns = pd.read_csv('matching_tmc2gmns.csv')
    measurement_tmc = pd.read_csv('measurement_tmc.csv')

    '''build measurement_base.csv''' 
    measurement_base = pd.DataFrame()
    measurement_base['link_id'] = None
    measurement_base['osm_way_id'] = None
    measurement_base['from_node_id'] = None
    measurement_base['to_node_id'] = None
    measurement_base['lanes'] = None
    measurement_base['length'] = None
    measurement_base['time_period'] = None
    measurement_base['date'] = None
    measurement_base['geometry'] = None
    measurement_base['volume'] = None
    measurement_base['speed'] = None
    measurement_base['ip_address'] = None

    k=0
    p=1
    measurement_base_dict = {}
    for i in matching_tmc2gmns.index:
        try:
            measurement_tmc_selected = measurement_tmc[measurement_tmc['link_id_tmc'] == matching_tmc2gmns['link_id_tmc'][i]]
            link_base_selected = link_base[link_base['link_id'] == matching_tmc2gmns['link_id_base'][i]]
            for j in measurement_tmc_selected.index:
                measurement_base_dict[k] = {'link_id': link_base_selected['link_id'].values[0],\
                                                'osm_way_id':link_base_selected['osm_way_id'].values[0],\
                                                'from_node_id': link_base_selected['from_node_id'].values[0],\
                                                'to_node_id': link_base_selected['to_node_id'].values[0],\
                                                'lanes': link_base_selected['lanes'].values[0], \
                                                'length': link_base_selected['length'].values[0], \
                                                'link_type_name': link_base_selected['link_type_name'].values[0], \
                                                'time_period': measurement_tmc_selected['time_period'][j],\
                                                'date': measurement_tmc_selected['date'][j],\
                                                'geometry': link_base_selected['geometry'].values[0],\
                                                'volume': None,\
                                                'speed': measurement_tmc_selected['speed'][j],\
                                                'ip_address': 'www.openstreetmap.org/?way=' + str(link_base_selected['osm_way_id'].values[0])} 

                k += 1
        except:
            measurement_base_dict[k] = {'link_id': link_base_selected['link_id'].values[0],\
                                                'osm_way_id':link_base_selected['osm_way_id'].values[0],\
                                                'from_node_id': link_base_selected['from_node_id'].values[0],\
                                                'to_node_id': link_base_selected['to_node_id'].values[0],\
                                                'lanes': link_base_selected['lanes'].values[0], \
                                                'length': link_base_selected['length'].values[0], \
                                                'link_type_name': link_base_selected['link_type_name'].values[0], \
                                                'time_period':None,\
                                                'date': None,\
                                                'geometry': link_base_selected['geometry'].values[0],\
                                                'volume': None,\
                                                'speed': None,\
                                                'ip_address': 'www.openstreetmap.org/?way=' + str(link_base_selected['osm_way_id'].values[0])}

            k += 1
        
        if i+1 > p/10 * len(matching_tmc2gmns.index): 
            print(str(p*10)+"%"+' measurement_base completed!')
            p = p + 1

    measurement_base = pd.DataFrame(measurement_base_dict).transpose()
    measurement_base.to_csv('measurement_base.csv',index = False)
    print('measurement_base.csv generated!')


def OneFunction():
    TMCIdentification2GMNSNodeLinkFiles('TMC_Identification.csv')
    ConvertTMCReading2Measurement('Reading.csv','link_tmc.csv')
    MatchTMC2GMNSNetwork('link_tmc.csv','link.csv')
    ConvertMeasurementBasedOnMatching('link.csv','matching_tmc2gmns.csv','measurement_tmc.csv')