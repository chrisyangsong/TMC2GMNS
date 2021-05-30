# @author       Kai (Frank) Zhang (xiaomo123zk@gmail.com)
# @time         2021/5/20 9:40
# @desc         [script description]

""" MapMatching4GMNS

Based on input network and given GPS trajectory data, the map-matching program
of Matching2Route aims to find most likely route in terms of node sequence in
the underlying network, with the following data flow chart.

The code is adopted and modified from
https://github.com/asu-trans-ai-lab/MapMatching4GMNS
"""
#!/usr/bin/python
# coding:utf-8

import threading
import time
import ctypes
import collections
import heapq
import os.path
from sys import platform
import shutil

if platform.startswith('win32'):
    _dll_file = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS.dll')
elif platform.startswith('linux'):
    _dll_file = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS.so')
elif platform.startswith('darwin'):
    _dll_file = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS.dylib')
else:
    raise Exception('Please build the shared library compatible to your OS\
                    using source files in engine_cpp!')


_cdll = ctypes.cdll.LoadLibrary(_dll_file)
# _cdll.MapMatching4GMNS.argtypes = [ctypes.c_int]

# print("To avoid complex data folder settings, please always first put the input data on the current directory.")


def map_match():
    '''
    input data: node.csv, link.csv and trace.csv
    output data: agent.csv and agent_performance.csv
    '''
    data_cleaning()
    print('call MapMatching4GMNS  dynamic library')
    print('\MapMatching4GMNS run starts')
    _cdll.MapMatching4GMNS()
    print('\MapMatching4GMNS run completes')


def data_cleaning():
    '''
    input_agent.csv
    agent.csv
    link_performance.csv
    '''
    # input_agent.csv
    '''
    There will be an extra comma(',') behind the geometry value in the output input_agent.csv 
    file due to different data specifications  by different programing languages. Therefore, we have to clean the data.
    '''
    if(os.path.exists("input_agent.csv")):
        f = open('input_agent.csv', encoding="utf-8")
        content = f.read()
        f.close()
        t = content.replace(",)", ")")  # delete the comma(',')
        with open("input_agent.CSV", "w", encoding='utf-8') as f1:
            f1.write(t)

    '''
    In the colab environment, if the output file already exists in the current directory, 
    an error will be reported when the newly generated output file is needed.
    '''
    # agent.csv
    if(os.path.exists("agent.csv")):
        os.remove("agent.csv")
    # link_performance.csv
    if(os.path.exists("link_performance.csv")):
        os.remove("link_performance.csv")


def _download_url(url, filename, loc_dir):
    try:
        import requests
    except ImportError:
        print('please print requests to preceed downloading!!')

    try:
        r = requests.get(url)
        r.raise_for_status()
        with open(loc_dir+filename, 'wb') as f:
            f.write(r.content)
    except requests.HTTPError:
        print('file not existing: '+url)
    except requests.ConnectionError:
        raise Exception(
            'check your connectcion!!! Then, if the input data exists, copy the data from local system')
    except Exception as e:
        raise e


def download_sample_data_sets_from_network():
    #              raw.githubusercontent.com/username/repo-name/branch-name/path
    url = 'https://raw.githubusercontent.com/xiaomo123zk/MapMatching4GMNS-0.2/master/data/'
    data_sets = [
        "testdata"
    ]

    files = [
        "node.csv",
        "link.csv",
        "trace.csv"
    ]

    print('downloading starts')

    # data folder under cdw
    loc_data_dir = 'data'
    if not os.path.isdir(loc_data_dir):
        os.mkdir(loc_data_dir)

    for ds in data_sets:
        web_dir = url + ds + '/'
        loc_sub_dir = os.path.join(loc_data_dir, ds) + '/'

        if not os.path.isdir(loc_sub_dir):
            os.mkdir(loc_sub_dir)

        # multi-threading
        threads = []
        for x in files:
            t = threading.Thread(
                target=_download_url,
                args=(web_dir+x, x, loc_sub_dir)
            )
            t.start()
            threads.append(t)

        for t in threads:
            t.join()

    print('downloading completes')

    print('check '+os.path.join(os.getcwd(), loc_data_dir) +
          ' for downloaded data sets')

    # then, copy the input data to current path
    #  os.path.join(os.getcwd(), 'bin/MapMatching4GMNS.dll')
    if(os.path.exists(os.getcwd()+"/data/testdata")):
        filePath = os.path.join(os.getcwd()+"/data/testdata")
        currentPath = os.path.join(os.getcwd())
        shutil.copy(filePath+'/node.csv', currentPath)
        shutil.copy(filePath+'/link.csv', currentPath)
        shutil.copy(filePath+'/trace.csv', currentPath)


if __name__ == "__main__":

    # first, check your operation system.
    # If you run the code on windows system without installed C++ environment,
    # some necessary dependency libraries need to be copied.

    # If an error occurs: Permission denied: 'C:/Windows/System32/*.dll',
    # you need to manually copy the dependency library from  https://github.com/xiaomo123zk/MapMatching4GMNS-0.2/tree/main/Dependent_libraries_missing_in_windows_system

    # First, download the input data of the test: node.csv, link.csv and trace.csv from Github.

    download_sample_data_sets_from_network()

    # If the online download fails, Please download manually the input data from https://github.com/asu-trans-ai-lab/osm_test_data_set/map_matching/.

    # 2, call the mapmatching4gmns library to calculate and output the result in the current directory.
    start = time.time()
    map_match()
    end = time.time()
    print('MapMatching4GMNS time cost: %.6f seconds' % (end - start))
    print("The output data is generated!")
