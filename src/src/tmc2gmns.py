# @author       Chris Yang Song (songyang0714@gmail.com)
# @time         2021/5/31 9:40
# @desc         [script description]

""" tmc2gmns
This script aims to transform tmc file into gmns format. Then after map-matching program
of MapMatching4GMNS, the link file generated from tmc file is matched to the underlying network (here is gmns format of osm map). 
In the end, the link performance file of underlying network is generated with tmc file and the corresponding reading file.
"""
#!/usr/bin/python
# coding:utf-8

import os
import datetime
import numpy as np
import pandas as pd

import os.path
import MapMatching4GMNS

'''step 1 Convert TMC Data into GMNS Format
Convert TMC Data into GMNS Format
'''

def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)


def Convert_TMC(tmc_path):
    '''build node_tmc.csv'''

    print('reading tmc data...')

    files= os.listdir(tmc_path) 
    for file in files: 
        if file[:18] == 'TMC_Identification':
            tmc = pd.read_csv(tmc_path + os.sep + file)
            break
    
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

    node_tmc.to_csv(tmc_path + os.sep + '/node_tmc.csv')
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


    link_tmc.to_csv(tmc_path + os.sep + '/link_tmc.csv')
    print('link_tmc.csv generated!')


    '''build link_performance_tmc.csv''' 

    reading = pd.read_csv(tmc_path + os.sep + 'Reading_VA.csv')
    # reading = reading[pd.to_datetime(reading['measurement_tstamp'], format='%Y-%m-%d %H:%M:%S')<datetime.datetime.strptime('2015-04-01 02:00:00', '%Y-%m-%d %H:%M:%S')]
    reading = reading.loc[0:2000]

    link_performance_tmc = pd.DataFrame()
    link_performance_tmc['name'] = None
    link_performance_tmc['corridor_id'] = None
    link_performance_tmc['corridor_link_order'] = None
    link_performance_tmc['from_node_id'] = None
    link_performance_tmc['to_node_id'] = None
    link_performance_tmc['timestamp'] = None
    link_performance_tmc['volume'] = None
    link_performance_tmc['travel_time'] = None
    link_performance_tmc['speed'] = None
    link_performance_tmc['reference_speed'] = None
    link_performance_tmc['density'] = None
    link_performance_tmc['queue'] = None
    link_performance_tmc['notes'] = None


    gp = reading.groupby('measurement_tstamp')
    for key, form in gp:
        # print(key)
        for i in link_tmc.index:
            form_selected = form[form['_vatmc_code']==link_tmc['name'][i]]
            if len(form_selected)>0:
            # break
                link_performance_tmc = link_performance_tmc.append({'name': link_tmc['name'][i],\
                                                'corridor_id': link_tmc['corridor_id'][i],\
                                                'corridor_link_order' : link_tmc['corridor_link_order'][i],\
                                                'from_node_id': link_tmc.loc[i,'from_node_id'], \
                                                'to_node_id': link_tmc.loc[i,'to_node_id'], \
                                                'timestamp': form_selected['measurement_tstamp'].values[0][0:10]+'T'+form_selected['measurement_tstamp'].values[0][11:13]+':'+form_selected['measurement_tstamp'].values[0][14:16],\
                                                'volume': None,\
                                                'travel_time': link_tmc['length'][i]/form_selected['speed'].values[0],\
                                                'speed': form_selected['speed'].values[0],\
                                                'reference_speed': form_selected['reference_speed'].values[0],\
                                                'density': None,\
                                                'queue': None,\
                                                'notes': None }, ignore_index=True)
            else:
                link_performance_tmc = link_performance_tmc.append({'name': link_tmc['name'][i],\
                                                'corridor_id': link_tmc['corridor_id'][i],\
                                                'corridor_link_order' : link_tmc['corridor_link_order'][i],\
                                                'from_node_id': link_tmc.loc[i,'from_node_id'], \
                                                'to_node_id': link_tmc.loc[i,'to_node_id'], \
                                                'timestamp': None,\
                                                'volume': None,\
                                                'travel_time': None,\
                                                'speed': None,\
                                                'reference_speed': None,\
                                                'density': None,\
                                                'queue': None,\
                                                'notes': None }, ignore_index=True)

    link_performance_tmc.to_csv(tmc_path + os.sep +'/link_performance_tmc.csv',index = False)
    print('link_performance_tmc.csv generated!')




    '''build trace.csv'''
    '''trace_id is numeric'''
    trace = pd.DataFrame()
    trace['corridor_id'] = None
    trace['agent_id'] = None
    trace['date'] = None
    trace['tmc'] = None
    trace['trace_id'] = None
    trace['hh'] = None
    trace['mm'] = None
    trace['ss'] = None
    trace['y_coord'] = None
    trace['x_coord'] = None

    agent_id = 1
    trace_id = 0

    for i in range(0,len(tmc)-1):
        if tmc.loc[i+1,'road_order'] > tmc.loc[i,'road_order']:
            trace = trace.append({'corridor_id': tmc.loc[i,'road'] + '_' + tmc.loc[i,'direction'],\
                                        'agent_id': agent_id,\
                                        'date': None, \
                                        'tmc': tmc.loc[i,'tmc'],\
                                        'trace_id': trace_id,\
                                        'hh': None,\
                                        'mm': 'None',\
                                        'ss': None,\
                                        'y_coord': tmc.loc[i,'start_latitude'],\
                                        'x_coord': tmc.loc[i,'start_longitude']}, ignore_index=True)
            trace_id +=1
        else:
            trace = trace.append({'corridor_id': tmc.loc[i,'road'] + '_' + tmc.loc[i,'direction'],\
                                        'agent_id': agent_id,\
                                        'date': None, \
                                        'tmc': tmc.loc[i,'tmc'],\
                                        'trace_id': trace_id,\
                                        'hh': None,\
                                        'mm': 'None',\
                                        'ss': None,\
                                        'y_coord': tmc.loc[i,'start_latitude'],\
                                        'x_coord': tmc.loc[i,'start_longitude']}, ignore_index=True)
            trace_id += 1

            trace = trace.append({'corridor_id': tmc.loc[i,'road'] + '_' + tmc.loc[i,'direction'],\
                                        'agent_id': agent_id,\
                                        'date': None, \
                                        'tmc': tmc.loc[i,'tmc'],\
                                        'trace_id': trace_id,\
                                        'hh': None,\
                                        'mm': 'None',\
                                        'ss': None,\
                                        'y_coord': tmc.loc[i,'end_latitude'],\
                                        'x_coord': tmc.loc[i,'end_longitude']}, ignore_index=True)
            agent_id += 1
            trace_id = 0
            
    trace = trace.append({'corridor_id': tmc.loc[i+1,'road'] + '_' + tmc.loc[i+1,'direction'],\
                                        'agent_id': agent_id,\
                                        'date': None, \
                                        'tmc': tmc.loc[i+1,'tmc'],\
                                        'trace_id': trace_id,\
                                        'hh': None,\
                                        'mm': 'None',\
                                        'ss': None,\
                                        'y_coord': tmc.loc[i+1,'start_latitude'],\
                                        'x_coord': tmc.loc[i+1,'start_longitude']}, ignore_index=True)
    trace_id +=1
    trace = trace.append({'corridor_id': tmc.loc[i+1,'road'] + '_' + tmc.loc[i+1,'direction'],\
                                        'agent_id': agent_id,\
                                        'date': None, \
                                        'tmc': tmc.loc[i+1,'tmc'],\
                                        'trace_id': trace_id,\
                                        'hh': None,\
                                        'mm': 'None',\
                                        'ss': None,\
                                        'y_coord': tmc.loc[i+1,'end_latitude'],\
                                        'x_coord': tmc.loc[i+1,'end_longitude']}, ignore_index=True)

    trace.to_csv(tmc_path + os.sep +'/trace.csv')
    print('trace.csv generated!')



'''step 2 map matching
'''
def mapmatch(tmc_path,osm_path):
    "Get the OSM Network"
    import osm2gmns as og
    files= os.listdir(osm_path) 
    for file in files: 
        if file[-3:] == 'osm':
            net = og.getNetFromOSMFile(osm_path + os.sep + file,network_type=('auto'), default_lanes=True, default_speed=True)
            og.consolidateComplexIntersections(net)
            og.outputNetToCSV(net, output_folder=osm_path)
        break
    
    create_folder(os.path.join(os.path.dirname(os.path.realpath('__file__')), 'data/testdata'))

    import shutil
    source_link_osm = os.path.join(os.path.join(os.path.dirname(os.path.realpath('__file__')),osm_path),'link.csv')
    source_node_osm = os.path.join(os.path.join(os.path.dirname(os.path.realpath('__file__')),osm_path),'node.csv')
    source_trace = os.path.join(os.path.join(os.path.dirname(os.path.realpath('__file__')),tmc_path),'trace.csv')
    source_list = [source_link_osm,source_node_osm,source_trace]

    destination_1 = os.path.join(os.path.dirname(os.path.realpath('__file__')), 'data/testdata')
    destination_2 = os.path.dirname(os.path.realpath('__file__'))
    destination_list = [destination_1,destination_2]
    
    for i in range(len(source_list)):
        for j in range(len(destination_list)):
            shutil.copy(source_list[i], destination_list[j])

    MapMatching4GMNS.map_match()
    print("Map matching is completed!")


'''step 3 link performance of osm
'''
def osm_link_performance(tmc_path,osm_path):
    link_road = pd.read_csv(osm_path + os.sep +'/link.csv', encoding='gbk',low_memory=False,index_col=None)
    trace = pd.read_csv(tmc_path + os.sep +'/trace.csv', encoding='gbk',low_memory=False,index_col=None)
    link_performance_tmc = pd.read_csv(tmc_path + os.sep +'/link_performance_tmc.csv', encoding='gbk',low_memory=False,index_col=None)
    mapping = pd.read_csv('link_performance.csv', encoding='gbk',low_memory=False,index_col=None)

    link_performance_osm = pd.DataFrame()
    link_performance_osm['link_id'] = None
    link_performance_osm['from_node_id'] = None
    link_performance_osm['to_node_id'] = None
    link_performance_osm['lanes'] = None
    link_performance_osm['length'] = None
    link_performance_osm['timestamp'] = None
    link_performance_osm['geometry'] = None
    link_performance_osm['volume'] = None
    link_performance_osm['speed'] = None


    link_performance_osm_sub = link_performance_osm

    gp = link_performance_tmc.groupby('timestamp')
    for key, form in gp:
        for i in mapping.index:
            link_road_selected = link_road[(link_road['from_node_id'] == mapping.loc[i,'from_node_id']) & (link_road['to_node_id'] == mapping.loc[i,'to_node_id'])]
            trace_selected = trace[(trace['agent_id'] == mapping.loc[i,'agent_id']) & (trace['trace_id'] == mapping.loc[i,'trace_id'])]
            form_selected = form[(form['name']==trace_selected['tmc'].values[0]) & (form['corridor_id']==trace_selected['corridor_id'].values[0])]
            if len(form_selected)>0:
                link_performance_osm_sub = link_performance_osm_sub.append({'link_id': link_road_selected['link_id'].values[0],\
                                                'from_node_id': link_road_selected['from_node_id'].values[0],\
                                                'to_node_id': link_road_selected['to_node_id'].values[0],\
                                                'lanes': link_road_selected['lanes'].values[0],\
                                                'length': link_road_selected['length'].values[0],\
                                                'timestamp': form_selected['timestamp'].values[0],\
                                                'geometry': link_road_selected['geometry'].values[0],\
                                                'volume': form_selected['volume'].values[0],\
                                                'speed': form_selected['speed'].values[0]}, ignore_index=True)
            else:
                link_performance_osm_sub = link_performance_osm_sub.append({'link_id': link_road_selected['link_id'].values[0],\
                                                'from_node_id': link_road_selected['from_node_id'].values[0],\
                                                'to_node_id': link_road_selected['to_node_id'].values[0],\
                                                'lanes': link_road_selected['lanes'].values[0],\
                                                'length': link_road_selected['length'].values[0],\
                                                'timestamp': None,\
                                                'geometry': link_road_selected['geometry'].values[0],\
                                                'volume': None,\
                                                'speed': None}, ignore_index=True)
        link_performance_osm_sub = link_performance_osm_sub.drop_duplicates(['from_node_id', 'from_node_id']) #avoid one osm link corresponding to mutiple tmc link

        link_performance_osm = link_performance_osm.append(link_performance_osm_sub)
        link_performance_osm_sub = pd.DataFrame()
    link_performance_osm.to_csv(tmc_path + os.sep +'/link_performance_osm.csv')
    print("link_performance_osm.csv is generated!")


# if __name__ == "__main__":
#     tmc_path = '/usr/local/home/ysx28/Desktop/GMNS/tmc2gmns/TMC2GMNS/src/tmc'
#     osm_path = '/usr/local/home/ysx28/Desktop/GMNS/tmc2gmns/TMC2GMNS/src/osm'

#     import tmc2gmns as tg
#     tg.Convert_TMC(tmc_path) 
#     tg.mapmatch(tmc_path,osm_path)
#     tg.osm_link_performance(tmc_path,osm_path)