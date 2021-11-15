from rasterstats import zonal_stats
from pathlib import Path
from shapely.wkt import loads
import geopandas as gpd
import pandas as pd
import numpy as np
import datetime as dt
import re
import json
import os


def main():
    index_list = ['SPI-1', 'SPI-2',
                  'SPI-3', 'SPI-6',
                  'SPI-12', 'SPEI-1',
                  'SPEI-2', 'SPEI-3',
                  'SPEI-6', 'SPEI-12',
                  'VHI', 'VCI',
                  'SMA', 'CDI']

    # import shapefile
    nuts3 = gpd.read_file('/mnt/CEPH_PROJECTS/ADO/GIS/ERTS89_LAEA/alpinespace_eusalp_NUTS3_simple.shp')

    for i in index_list:
        print(i)
        compute_nuts3_stats(nuts3, i)

    merge_time_series(nuts3, index_list)


def merge_time_series(nuts3, index_list):
    # iterate through nuts zones
    for ifeat in nuts3.iterfeatures():
        iid = ifeat['properties']['NUTS_ID']

        ts_list = list()
        # iterate through indices
        for iind in index_list:
            ts_inpath = '../json/timeseries/NUTS3_' + iid + '_tmp' + iind + '.json'
            tmp = pd.read_json(ts_inpath)
            tmp = tmp.set_index('date')
            ts_list.append(tmp)
            os.remove(ts_inpath)

        # concatenate time series
        ts_out = pd.concat(ts_list)
        # reformat for output
        ts_outpath = '../json/timeseries/NUTS3_' + iid + '.json'
        ts_list = list()
        for i in ts_out.iterrows():
            ts_list.append({**{'date': i[0].strftime('%Y-%m-%d')}, **i[1].to_dict()})
        # write time series
        with open(ts_outpath, 'w') as outfile:
            json.dump(ts_list, outfile)


def compute_nuts3_stats(nuts3, index='SPI3'):
    today = dt.datetime.strptime('2018-08-24', '%Y-%m-%d') #TODO: modify to use actual date of the day the script is run
    drange = pd.date_range(today - dt.timedelta(days=365), today, freq='1D')

    basepath = get_basepath(index)
    path_list = list()
    for path in basepath.glob('*.tif'):
        path_list.append(path)

    # sort list of paths
    path_list = sorted(path_list, key=lambda i: parse_index_date(index, i))

    # create a working copy
    tmp_nuts3 = nuts3.copy()

    # add SPI property
    dictlist = [dict() for i in range(nuts3.shape[0])]
    tmp_nuts3[index] = dictlist

    # calculate statistics
    for i_date in drange:
        no_data = True
        for i in path_list:
            if compare_date(index, parse_index_date(index, i), i_date):
                i_path = i
                no_data = False
                break

        # date = parse_index_date(index, i_path)
        # if no data for time step, fill with None, else compute stats
        if no_data:
            for ifeat in tmp_nuts3.iterfeatures():
                ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = None
        else:
            r_stats = get_stats(i_path, index)
            print(i_date)
            for ifeat, istat in zip(tmp_nuts3.iterfeatures(), r_stats):
                if np.isinf(istat['nanmean']) or istat['count'] == 0:
                    ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = None
                else:
                    ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = round(istat['nanmean'].astype(float), 3)

    # date range to be deleted after time series export
    delrange = pd.date_range(today-dt.timedelta(days=365), today-dt.timedelta(days=335))

    # create time series per nuts region
    for ifeat in tmp_nuts3.iterfeatures():
        # check if file exists
        ts_outpath = '../json/timeseries/NUTS3_' + ifeat['properties']['NUTS_ID'] + '_tmp' + index + '.json'
        ts_list = list()
        for idate in ifeat['properties'][index].items():
            ts_list.append({'date': idate[0], index: idate[1]})
        # write time series
        with open(ts_outpath, 'w') as outfile:
            json.dump(ts_list, outfile)

        # remove dates from nuts3 maps
        for ddate in delrange:
            del ifeat['properties'][index][ddate.date().strftime('%Y-%m-%d')]

    tmp_nuts3_4326 = round_coordinates(tmp_nuts3.to_crs("EPSG:4326"))
    tmp_nuts3_4326.to_file('../json/' + index + '-latest.geojson', driver='GeoJSON', encoding='utf-8')


def compare_date(index, index_date, ts_date):
    if index == 'VCI' or index == 'VHI':
        doy4step = (ts_date.dayofyear // 4) * 4 + 4
        if doy4step > 364:
            doy4step = 364
        date4step = dt.datetime.strptime(str(ts_date.year) + '_' + f"{doy4step:03d}", '%Y_%j')
        return index_date == date4step
    else:
        return index_date == ts_date


def parse_index_date(index, path):
    if index == 'SPI-6' or index == 'SPI-12':
        index_date = dt.datetime.strptime(path.name[7:15], '%Y%m%d')
    elif index[0:3] == 'SPI':
        index_date = dt.datetime.strptime(path.name[6:14], '%Y%m%d')
    elif index == 'SPEI-6' or index == 'SPEI-12':
        index_date = dt.datetime.strptime(path.name[8:16], '%Y%m%d')
    elif index[0:3] == 'SPE':
        index_date = dt.datetime.strptime(path.name[7:15], '%Y%m%d')
    elif index[0:3] == 'VHI' or index[0:3] == 'VCI':
        index_date = dt.datetime.strptime(path.name[7:15], '%Y_%j')
    elif index[0:3] == 'SMA':
        index_date = dt.datetime.strptime(path.name[15:23], '%Y%m%d')
    elif index[0:3] == 'CDI':
        index_date = dt.datetime.strptime(path.name[4:12], '%Y%m%d')

    return index_date


def get_basepath(index):
    if index[0:3] == 'SPI' or index[0:3] == 'SPE':
        basepath = '/mnt/CEPH_PROJECTS/ADO/ARSO/ERA5_QM_NEW/' + index + '/'
    elif index[0:3] == 'VHI':
        basepath = '/mnt/CEPH_PROJECTS/ADO/VHI/03_results/vhi/'
    elif index[0:3] == 'VCI':
        basepath = '/mnt/CEPH_PROJECTS/ADO/VHI/03_results/vci/'
    elif index[0:3] == 'SMA':
        basepath = '/mnt/CEPH_PROJECTS/ADO/SM/ERA5_ERA5l_QM/anomalies_long_ref/'
    elif index[0:3] == 'CDI':
        basepath = '/mnt/CEPH_PROJECTS/ADO/CDI/'

    return Path(basepath)


def get_stats(r_path, index):
    if index == 'CDI':
        r_stats = zonal_stats("/mnt/CEPH_PROJECTS/ADO/GIS/ERTS89_LAEA/alpinespace_eusalp_NUTS3_simple.shp",
                              r_path,
                              add_stats={'nanmean': most_frequent}, nodata=np.nan)
    else:
        r_stats = zonal_stats("/mnt/CEPH_PROJECTS/ADO/GIS/ERTS89_LAEA/alpinespace_eusalp_NUTS3_simple.shp",
                              r_path,
                              add_stats={'nanmean': mymean})
    return r_stats


def mymean(x):
    return np.nanmean(x)


def most_frequent(x):
    return np.inf if x.compressed().size == 0 else np.argmax(np.bincount(x.compressed().astype(np.int64)))


def round_coordinates(df):
    simpledec = re.compile(r"\d*\.\d+")

    def mround(match):
        return "{:.5f}".format(float(match.group()))

    df.geometry = df.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

    return df


if __name__ == "__main__":
    main()