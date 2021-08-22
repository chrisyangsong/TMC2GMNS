import os
import math
import datetime
import numpy as np
import pandas as pd

import time
import os.path

import shapely.geometry as geom
import geopandas as gpd
from shapely.geometry import LineString
from shapely.geometry import MultiLineString
from shapely.geometry import Point



def TMCIdentification2GMNSNodeLinkFiles(TMC_file,link_base):#output:node_tmc,link_tmc
    
    '''reading link_base'''
    link_base = pd.read_csv(link_base, low_memory=False)
    # link_base = link_base[link_base['link_type_name'].isin(['motorway','trunk'])]
    # link_base = link_base[-link_base['name'].isna()]
    # link_base = link_base[link_base['link_type_name'].isin(['motorway','trunk','primary','secondary'])]
    # link_base = link_base.reset_index()
    # link_base = link_base.drop(['index'], 1)

    '''convert link_base to MultiLineString''' 
    multiline_string_base_list = []
    multiline_string_base_list_sub = []
    for j in link_base.index:
        link_base_geometry_list = link_base.loc[j,'geometry'][12:-1].split(", ")
        for link_base_geometry in link_base_geometry_list:
            multiline_string_base_list_sub.append((float(link_base_geometry.split(" ")[0]),float(link_base_geometry.split(" ")[1])))
        multiline_string_base_list_sub = tuple(multiline_string_base_list_sub)
        multiline_string_base_list.append(multiline_string_base_list_sub)
        multiline_string_base_list_sub = []

    from shapely.geometry import MultiLineString
    line_base = MultiLineString(multiline_string_base_list) 

    '''reading tmc'''
    tmc = pd.read_csv(TMC_file)
    tmc = tmc.drop_duplicates(subset=['direction','road_order']).sort_values(by=['direction','road_order'])
    tmc = tmc.reset_index()
    tmc = tmc.drop(['index'], 1)
    origin_tmc_num = len(tmc)

    '''remove out boundary tmc'''
    in_bbox_index_list = []
    for i in tmc.index:
        if (tmc['start_longitude'][i] > line_base.bounds[0]) & (tmc['start_longitude'][i] < line_base.bounds[2]) & \
            (tmc['end_longitude'][i] > line_base.bounds[0]) & (tmc['end_longitude'][i] < line_base.bounds[2]) & \
                (tmc['start_latitude'][i] > line_base.bounds[1]) & (tmc['start_latitude'][i] < line_base.bounds[3]) & \
            (tmc['end_latitude'][i] > line_base.bounds[1]) & (tmc['end_latitude'][i] < line_base.bounds[3]):
            in_bbox_index_list.append(i)

    tmc = tmc.loc[in_bbox_index_list]
    tmc = tmc.reset_index()
    tmc = tmc.drop(['index'], 1)
    if len(in_bbox_index_list) < origin_tmc_num:
        print('base map cannot cover all TMC nodes,' + str(origin_tmc_num-len(in_bbox_index_list)) + 'tmc nodes are out of boundary box, please use larger base map')

    form_tmc = pd.DataFrame([])
    gp = tmc.groupby(['road','direction'])
    for key, form in gp:
        form['end_latitude'][:-1] = form['start_latitude'][1:].tolist()
        form['end_longitude'][:-1] = form['start_longitude'][1:].tolist()
        form_tmc = pd.concat([form_tmc, form], axis=0)
    form_tmc = form_tmc.reset_index()
    form_tmc = form_tmc.drop(['index'], 1)
    tmc = form_tmc


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

    gp = tmc.groupby(['road','direction'])
    for key, form in gp:
        for i in form.index:
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



    # for i in range(0,len(tmc)-1):
    #     if tmc.loc[i+1,'road_order'] > tmc.loc[i,'road_order']:
    #         node_tmc = node_tmc.append({'name': tmc.loc[i,'tmc'],\
    #                                     'x_coord': tmc.loc[i,'start_longitude'], \
    #                                     'y_coord': tmc.loc[i,'start_latitude'],\
    #                                     'z_coord': None,\
    #                                     'node_type': 'tmc_start',\
    #                                     'ctrl_type': None,\
    #                                     'zone_id': None,\
    #                                     'geometry': "POINT (" + tmc.loc[i,'start_longitude'].astype(str) + " " + tmc.loc[i,'start_latitude'].astype(str) +")"}, ignore_index=True)
    #     else:
    #         node_tmc = node_tmc.append({'name': tmc.loc[i,'tmc'],\
    #                                     'x_coord': tmc.loc[i,'start_longitude'], \
    #                                     'y_coord': tmc.loc[i,'start_latitude'],\
    #                                     'z_coord': None,\
    #                                     'node_type': 'tmc_start',\
    #                                     'ctrl_type': None,\
    #                                     'zone_id': None,\
    #                                     'geometry': "POINT (" + tmc.loc[i,'start_longitude'].astype(str) + " " + tmc.loc[i,'start_latitude'].astype(str) +")"}, ignore_index=True)
    #         node_tmc = node_tmc.append({'name': tmc.loc[i,'tmc']+'END',\
    #                                     'x_coord': tmc.loc[i,'end_longitude'], \
    #                                     'y_coord': tmc.loc[i,'end_latitude'],\
    #                                     'z_coord': None,\
    #                                     'node_type': 'tmc_end',\
    #                                     'ctrl_type': None,\
    #                                     'zone_id': None,\
    #                                     'geometry': "POINT (" + tmc.loc[i,'end_longitude'].astype(str) + " " + tmc.loc[i,'end_latitude'].astype(str) +")"}, ignore_index=True)

    #     if i > p/10 * len(tmc): 
    #         print(str(p*10)+"%"+' nodes completed!')
    #         p = p + 1

    # node_tmc = node_tmc.append({'name': tmc.loc[i+1,'tmc'],\
    #                                     'x_coord': tmc.loc[i+1,'start_longitude'], \
    #                                     'y_coord': tmc.loc[i+1,'start_latitude'],\
    #                                     'z_coord': None,\
    #                                     'node_type': 'tmc_start',\
    #                                     'ctrl_type': None,\
    #                                     'zone_id': None,\
    #                                     'geometry': "POINT (" + tmc.loc[i+1,'start_longitude'].astype(str) + " " + tmc.loc[i+1,'start_latitude'].astype(str) +")"}, ignore_index=True)

    # node_tmc = node_tmc.append({'name': tmc.loc[i+1,'tmc']+'END',\
                                        # 'x_coord': tmc.loc[i+1,'end_longitude'], \
                                        # 'y_coord': tmc.loc[i+1,'end_latitude'],\
                                        # 'z_coord': None,\
                                        # 'node_type': 'tmc_end',\
                                        # 'ctrl_type': None,\
                                        # 'zone_id': None,\
                                        # 'geometry': "POINT (" + tmc.loc[i+1,'end_longitude'].astype(str) + " " + tmc.loc[i+1,'end_latitude'].astype(str) +")"}, ignore_index=True)

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
                                    # 'from_node_id': node_tmc[(node_tmc['x_coord']==tmc.loc[i,'start_longitude']) & (node_tmc['y_coord']==tmc.loc[i,'start_latitude'])].index.values[0], \
                                    # 'to_node_id': node_tmc[(node_tmc['x_coord']==tmc.loc[i,'end_longitude'])&(node_tmc['y_coord']==tmc.loc[i,'end_latitude'])].index.values[0],\
                                    'from_node_id': node_tmc[(node_tmc['y_coord']==tmc.loc[i,'start_latitude'])].index.values[0], \
                                    'to_node_id': node_tmc[(node_tmc['y_coord']==tmc.loc[i,'end_latitude'])].index.values[0],\
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


def volume_speed_func(speed,ffs=69,critical_density=37.5,mm=4): # fundamental diagram  (volume_delay fuction) 
    if speed > ffs:
        speed = ffs
    kernal=np.power(np.power(ffs/speed,mm),0.5)
    return speed*critical_density*np.power(kernal-1,1/mm)

def ConvertTMCReading2Measurement(Reading,link_tmc):

    link_tmc = pd.read_csv(link_tmc)
    ## reading by detid
    reading = pd.read_csv(Reading)
    reading = reading.loc[0:5000]
    reading_dict = {}
    gp = reading.groupby('tmc_code')
    for key, form in gp:
        reading_dict[key] = {
            'measurement_tstamp':form['measurement_tstamp'].tolist(),
            'speed':form['speed'].tolist()
            }
    reading_dict_selected = reading_dict[list(reading_dict.keys())[0]]
    delta_reading_time = divmod((pd.to_datetime(reading_dict_selected['measurement_tstamp'][1], format='%Y-%m-%d %H:%M:%S') - pd.to_datetime(reading_dict_selected['measurement_tstamp'][0], format='%Y-%m-%d %H:%M:%S')).seconds,60)[0]

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
                speed = round(np.mean(reading_dict_selected['speed'][j:j+2]))
                volume_rate = round(volume_speed_func(speed))
                volume = round(volume_rate/60*delta_reading_time)
                density = round(volume_rate/speed)
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
                                                'volume': volume,\
                                                'volume_rate': volume_rate,\
                                                'travel_time': None,\
                                                'speed':speed,\
                                                'reference_speed': None,\
                                                'density': density,\
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
                                            'volume_rate': volume_rate,\
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
    link_tmc = pd.read_csv(link_tmc)
    link_base = pd.read_csv(link_base, low_memory=False)
    # link_base = link_base[-link_base['name'].isna()]
    # link_base = link_base[link_base['link_type_name'].isin(['motorway','trunk','primary','secondary'])]
    # link_base = link_base[link_base['link_type_name'].isin(['motorway','trunk'])]
    # link_base = link_base.reset_index()
    # link_base = link_base.drop(['index'], 1)

    '''initial matching'''
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

        for dividing_rate in range(1,11):
            end_point_tmc = LineString([(tmc_lon_1, tmc_lat_1),  (tmc_lon_2, tmc_lat_2)]).interpolate(dividing_rate/10, normalized=True)
            start_point_tmc = LineString([(tmc_lon_1, tmc_lat_1),  (tmc_lon_2, tmc_lat_2)]).interpolate((dividing_rate-1)/10, normalized=True)

            distance_list = []
            angle_list = []
            for i in range(len(link_base)):
                lon_list = []
                lat_list = [] 
                point_base_list = []
                link_geometry_list = link_base.loc[i,'geometry'][12:-1].split(", ")
                for link_geometry in link_geometry_list:
                    lon_list.append(float(link_geometry.split(" ")[0]))
                    lat_list.append(float(link_geometry.split(" ")[1]))
                    point_base_list.append(tuple([float(link_geometry.split(" ")[0]),float(link_geometry.split(" ")[1])]))
                base_line = geom.LineString(tuple(point_base_list))
                '''distance'''
                distance_list.append(start_point_tmc.distance(base_line))
                
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

            small_angle_list = [i for i, value in enumerate(angle_list) if value < 45]
            df_distance = pd.DataFrame({'distance':distance_list})
            
            nearest_index = df_distance.loc[small_angle_list].idxmin().values[0]

            matching_tmc2gmns_dict[k] = {'name_tmc':link_tmc.loc[j]['name'],\
                                        'corridor_id_tmc':link_tmc.loc[j]['corridor_id'],\
                                        'link_id_tmc':link_tmc.loc[[j]].index.values[0],\
                                        'from_node_id_tmc':link_tmc.loc[j]['from_node_id'],\
                                        'to_node_id_tmc':link_tmc.loc[j]['to_node_id'],\
                                        'category_id_tmc':link_tmc.index.get_loc(j)+1,\
                                        'geometry_tmc':link_tmc.loc[j]['geometry'],\
                                        # 'name_base':link_base['name'][nearest_index],\
                                        'link_id_base':link_base['link_id'][nearest_index],\
                                        'from_node_id_base':link_base['from_node_id'][nearest_index],\
                                        'to_node_id_base':link_base['to_node_id'][nearest_index],\
                                        'category_id_base':link_tmc.index.get_loc(j)+1,\
                                        'geometry_base':link_base['geometry'][nearest_index],\
                                        'distance':min(distance_list),\
                                        'geometry_tmc_base':'MULTILINESTRING ('+ link_tmc.loc[j]['geometry'][11:] + \
                                                            ', ' + link_base['geometry'][nearest_index][11:]+')'}
            k += 1

        
        if link_tmc.index.get_loc(j) > p/10 * len(link_tmc): 
            print(str(p*10)+"%"+' initial matching completed!')
            p = p + 1
            

    matching_tmc2gmns = pd.DataFrame(matching_tmc2gmns_dict).transpose()
    matching_tmc2gmns_drop_duplicates = matching_tmc2gmns.drop(columns=['distance']).drop_duplicates()
    matching_tmc2gmns_drop_duplicates['distance'] = matching_tmc2gmns.loc[matching_tmc2gmns_drop_duplicates.index]['distance']
    matching_tmc2gmns_drop_duplicates = matching_tmc2gmns_drop_duplicates.reset_index()
    matching_tmc2gmns_drop_duplicates = matching_tmc2gmns_drop_duplicates.drop(['index'], 1)
    matching_tmc2gmns = matching_tmc2gmns_drop_duplicates

    # matching_tmc2gmns.to_csv('matching_tmc2gmns.csv',index = False)
    # print('matching_tmc2gmns.csv generated!')

    matching_name_base_counts = matching_tmc2gmns['name_base'].value_counts()
    link_base = link_base[link_base['name'].isin(matching_name_base_counts[matching_name_base_counts>30].index)]
    link_base = link_base.reset_index()
    link_base = link_base.drop(['index'], 1)
    # link_base.to_csv('link_base_small.csv',index = False)


    ''' extract base route '''
    link_base_iteration = link_base 
    base_route_dict = {}
    i = 0
    while len(link_base_iteration)>0:
        k = link_base_iteration.index[0]
        base_route_index_list = []
        base_route_list = []
        base_route_index_list.append(k)
        base_route_list.append(link_base.loc[k,'from_node_id'])
        base_route_list.append(link_base.loc[k,'to_node_id'])
        link_base_iteration = link_base_iteration.drop(k)

        selected_link_base_row = link_base_iteration[link_base_iteration['from_node_id'] == base_route_list[-1]]
        while len(selected_link_base_row) != 0:
            base_route_list.append(selected_link_base_row['to_node_id'].values[0])
            base_route_index_list.append(selected_link_base_row.index.values[0])
            link_base_iteration = link_base_iteration.drop(selected_link_base_row.index.values[0])
            selected_link_base_row = link_base_iteration[link_base_iteration['from_node_id'] == base_route_list[-1]]

        selected_link_base_row = link_base_iteration[link_base_iteration['to_node_id'] == base_route_list[0]]
        while len(selected_link_base_row) != 0:
            base_route_list.insert(0,selected_link_base_row['from_node_id'].values[0])
            base_route_index_list.insert(0,selected_link_base_row.index.values[0])
            link_base_iteration = link_base_iteration.drop(selected_link_base_row.index.values[0])
            selected_link_base_row = link_base_iteration[link_base_iteration['to_node_id'] == base_route_list[0]]

        base_route_dict[i] = {
            'base_route_list':base_route_list,
            'base_route_index_list':base_route_index_list
        }
        i+=1


    ''' extract tmc route'''
    tmc_route_dict = {}
    gp = matching_tmc2gmns.groupby('corridor_id_tmc')
    for key, form in gp:
        tmc_route_dict[key] = {
            'from_node_id_base':set(form['from_node_id_base'].tolist()),
            'tmc_index':list(form.index.values)
            }


    base_route_tmc_dict = {}
    k=0
    for i in range(len(base_route_dict)):
        for key, value in tmc_route_dict.items():
            base_route_tmc_dict[k] = {'base_route':i,
            'corridor': key,
            'tmc_index':value['tmc_index'],
            'number': int(len(set(base_route_dict[i]['base_route_list']).intersection(value['from_node_id_base'])))}
            k+=1

    base_route_tmc = pd.DataFrame(base_route_tmc_dict).transpose()     

    base_tmc_dict = {}
    gp = base_route_tmc.groupby('base_route')
    for key, form in gp:
        base_tmc_dict[key] = {
            'corridor':form.loc[pd.to_numeric(form['number']).idxmax()]['corridor']
            }
    base_tmc = pd.DataFrame(base_tmc_dict).transpose() 
    tmc_base_dict = {}
    gp = base_tmc.groupby('corridor')
    for key, form in gp:
        tmc_base_dict[key] = list(form.index.values)


    '''second matching'''
    matching_tmc2gmns_dict = {}
    k = 0
    p = 1
    tmc_count = 0
    for key, value in tmc_base_dict.items():
        link_tmc_sub = link_tmc[link_tmc['corridor_id'] == key]
        base_link_list_sub = []
        for base_route_id in value:
            base_link_list_sub += base_route_dict[base_route_id]['base_route_index_list']
        link_base_sub = link_base.loc[base_link_list_sub]


        for j in link_tmc_sub.index:
            lon_tmc_list = []
            lat_tmc_list = []
            link_tmc_geometry_list = link_tmc_sub.loc[j,'geometry'][12:-1].split(",")
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

            for dividing_rate in range(1,51):
                end_point_tmc = LineString([(tmc_lon_1, tmc_lat_1),  (tmc_lon_2, tmc_lat_2)]).interpolate(dividing_rate/50, normalized=True)
                start_point_tmc = LineString([(tmc_lon_1, tmc_lat_1),  (tmc_lon_2, tmc_lat_2)]).interpolate((dividing_rate-1)/50, normalized=True)

                distance_list = []
                angle_list = []
                for i in link_base_sub.index:
                    lon_list = []
                    lat_list = [] 
                    point_base_list = []
                    link_geometry_list = link_base_sub.loc[i,'geometry'][12:-1].split(", ")
                    for link_geometry in link_geometry_list:
                        lon_list.append(float(link_geometry.split(" ")[0]))
                        lat_list.append(float(link_geometry.split(" ")[1]))
                        point_base_list.append(tuple([float(link_geometry.split(" ")[0]),float(link_geometry.split(" ")[1])]))
                    base_line = geom.LineString(tuple(point_base_list))
                    '''distance'''
                    distance_list.append(start_point_tmc.distance(base_line))
                    
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

                small_angle_list = [i for i, value in enumerate(angle_list) if value < 45]
                df_distance = pd.DataFrame({'distance':distance_list})
                
                nearest_index = df_distance.loc[small_angle_list].idxmin().values[0]

                matching_tmc2gmns_dict[k] = {'name_tmc':link_tmc_sub.loc[j]['name'],\
                                            'corridor_id_tmc':link_tmc_sub.loc[j]['corridor_id'],\
                                            'link_id_tmc':link_tmc_sub.loc[[j]].index.values[0],\
                                            'from_node_id_tmc':link_tmc_sub.loc[j]['from_node_id'],\
                                            'to_node_id_tmc':link_tmc_sub.loc[j]['to_node_id'],\
                                            'category_id_tmc':link_tmc_sub.index.get_loc(j)+1,\
                                            'geometry_tmc':link_tmc_sub.loc[j]['geometry'],\
                                            # 'name_base':link_base_sub.iloc[nearest_index]['name'],\
                                            'link_id_base':link_base_sub.iloc[nearest_index]['link_id'],\
                                            'from_node_id_base':link_base_sub.iloc[nearest_index]['from_node_id'],\
                                            'to_node_id_base':link_base_sub.iloc[nearest_index]['to_node_id'],\
                                            'category_id_base':link_tmc_sub.index.get_loc(j)+1,\
                                            'geometry_base':link_base_sub.iloc[nearest_index]['geometry'],\
                                            'distance':min(distance_list),\
                                            'geometry_tmc_base':'MULTILINESTRING ('+ link_tmc_sub.loc[j]['geometry'][11:] + \
                                                                ', ' + link_base_sub.iloc[nearest_index]['geometry'][11:]+')'}
                k += 1

            tmc_count += 1
            if tmc_count > p/10 * len(link_tmc): 
                print(str(p*10)+"%"+' matching completed!')
                p = p + 1
            

    matching_tmc2gmns = pd.DataFrame(matching_tmc2gmns_dict).transpose()
    matching_tmc2gmns_drop_duplicates = matching_tmc2gmns.drop(columns=['distance']).drop_duplicates()
    matching_tmc2gmns_drop_duplicates['distance'] = matching_tmc2gmns.loc[matching_tmc2gmns_drop_duplicates.index]['distance']
    matching_tmc2gmns_drop_duplicates = matching_tmc2gmns_drop_duplicates.reset_index()
    matching_tmc2gmns_drop_duplicates = matching_tmc2gmns_drop_duplicates.drop(['index'], 1)
    matching_tmc2gmns = matching_tmc2gmns_drop_duplicates

    matching_tmc2gmns.to_csv('matching_tmc2gmns.csv',index = False)
    print('matching_tmc2gmns.csv generated!')


def ConvertMeasurementBasedOnMatching(link_base,matching_tmc2gmns,measurement_tmc):
    link_base = pd.read_csv(link_base)
    # link_base = link_base[link_base['link_type_name'].isin(['motorway','trunk'])]
    # link_base = link_base.reset_index()
    # link_base = link_base.drop(['index'], 1)
    matching_tmc2gmns = pd.read_csv(matching_tmc2gmns)
    measurement_tmc = pd.read_csv(measurement_tmc)

    '''build measurement_base.csv''' 

    matching_tmc2gmns_dict = {}
    gp = matching_tmc2gmns.groupby('link_id_base')
    for key, form in gp:
        matching_tmc2gmns_dict[key] = {
            'link_id_tmc':form['link_id_tmc'].tolist()
            }


    k=0
    p=1
    i=0
    measurement_base_dict = {}
    for key, value in matching_tmc2gmns_dict.items():
        i += 1
        try:
            link_base_selected = link_base[link_base['link_id'] == key]
            measurement_tmc_selected = measurement_tmc[measurement_tmc['link_id_tmc'].isin(value['link_id_tmc'])]
            gp_measurement_tmc_selected = measurement_tmc_selected.groupby(['time_period','date'])
            for key_measurement_tmc_selected, form_measurement_tmc_selected in gp_measurement_tmc_selected:
                measurement_base_dict[k] = {'link_id': link_base_selected['link_id'].values[0],\
                                                'osm_way_id':link_base_selected['osm_way_id'].values[0],\
                                                'from_node_id': link_base_selected['from_node_id'].values[0],\
                                                'to_node_id': link_base_selected['to_node_id'].values[0],\
                                                'lanes': link_base_selected['lanes'].values[0], \
                                                'length': link_base_selected['length'].values[0], \
                                                'link_type_name': link_base_selected['link_type_name'].values[0], \
                                                'corridor_id': form_measurement_tmc_selected['corridor_id'].tolist()[0],\
                                                'time_period': key_measurement_tmc_selected[0],\
                                                'date': key_measurement_tmc_selected[1],\
                                                'geometry': link_base_selected['geometry'].values[0],\
                                                'volume': round(form_measurement_tmc_selected['volume'].mean()),\
                                                'volume_rate': round(form_measurement_tmc_selected['volume_rate'].mean()),\
                                                'speed': round(form_measurement_tmc_selected['speed'].mean()),\
                                                'density': round(form_measurement_tmc_selected['density'].mean()),\
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
                                                'corridor_id': None,\
                                                'time_period':None,\
                                                'date': None,\
                                                'geometry': link_base_selected['geometry'].values[0],\
                                                'volume': None,\
                                                'volume_rate': None,\
                                                'speed': None,\
                                                'density': None,\
                                                'ip_address': 'www.openstreetmap.org/?way=' + str(link_base_selected['osm_way_id'].values[0])}

            k += 1
        
        if i+1 > p/10 * len(matching_tmc2gmns_dict): 
            print(str(p*10)+"%"+' measurement_base completed!')
            p = p + 1

    measurement_base = pd.DataFrame(measurement_base_dict).transpose()
    measurement_base.to_csv('measurement_base.csv',index = False)
    print('measurement_base.csv generated!')

def OneFunction():
    TMCIdentification2GMNSNodeLinkFiles('TMC_Identification.csv','link.csv')
    ConvertTMCReading2Measurement('Reading.csv','link_tmc.csv')
    MatchTMC2GMNSNetwork('link_tmc.csv','link.csv')
    ConvertMeasurementBasedOnMatching('link.csv','matching_tmc2gmns.csv','measurement_tmc.csv')