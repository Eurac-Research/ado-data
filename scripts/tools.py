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
import psycopg2


def generate_metadata(index, drange):
    geodf = dict()

    geodf['short_name'] = index
    # add index long name
    if index == 'CDI':
        geodf['long_name'] = 'Combined Drought Index'
    elif index == 'VHI':
        geodf['long_name'] = 'Vegetation Health Index'
    elif index == 'VCI':
        geodf['long_name'] = 'Vegetation Condition Index'
    elif index[0:3] == 'SPI':
        geodf['long_name'] = 'Standardized Precipitation Index - ' + index[4::]
    elif index[0:4] == 'SPEI':
        geodf['long_name'] = 'Standardaized Precipitation and Evapotranspiration Index - ' + index[5::]
    elif index == 'SMA':
        geodf['long_name'] = 'Soil Moisture Anomalies'
    elif index == 'SDI':
        geodf['long_name'] = 'Streamflow Drought Index'
    elif index == 'SSPI-10':
        geodf['long_name'] = 'Standardized Snowpack Index'

    # add start and end dates as well as the colormaps
    timerange = {"type": "date",
                 "properties": {
                     "firstDate": drange[0].strftime("%Y-%m-%d"),
                     "lastDate": drange[-1].strftime("%Y-%m-%d")
                 }}

    geodf['timerange'] = timerange

    if index[0:2] == 'SP':
        cmapname = 'SPI_SPEI'
    elif index[0] == 'V':
        cmapname = 'VCI_VHI'
    elif index == 'SDI':
        cmapname = 'SPI_SPEI'
    elif index == 'SSPI-10':
        cmapname = "SSPI"
    else:
        cmapname = index
    # get colormap
    cmap = pd.read_csv('../visualization/' + cmapname + '_colormap_hex.txt',
                       header=None,
                       index_col=False,
                       skiprows=2,
                       parse_dates=False)
    colormap = {"id": "data",
                "type": "fill",
                "paint": {
                    "fill-color": {
                        "property": "value",
                        "stops": [x[1].iloc[[0, 1]].to_list() for x in cmap.iterrows()]
                    }
                },
                "legend": {
                    "stops": [x[1].iloc[[0, 2, 1]].to_list() for x in cmap.iterrows()]
                }}
    geodf['colormap'] = colormap

    with open('../json/metadata/' + index + '.json', 'w') as json_file:
        json.dump(geodf, json_file, indent=4)


def update_metadata_drange(index, drange):
    json_path = '../json/hydro/metadata/' + index + '.json'
    with open(json_path) as json_file:
        data = json.load(json_file)

    data['timerange']['properties']['firstDate'] = (drange[-1] - dt.timedelta(days=30)).strftime("%Y-%m-%d")
    data['timerange']['properties']['lastDate'] = drange[-1].strftime("%Y-%m-%d")

    with open(json_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)


def compute_stats(drange, shapes, points=None, index='SPI3', aggunit='hydro', update=False):
    # find all index files
    basepath = get_basepath(index)

    if index != 'SDI':
        path_list = list()
        for path in basepath.glob('*.tif'):
            path_list.append(path)

        # sort list of paths
        path_list = sorted(path_list, key=lambda i: parse_index_date(index, i))

    if update:
        tmp_shapes = shapes
    else:
        # create a working copy
        tmp_shapes = shapes.copy()

        # add index property
        dictlist = [dict() for i in range(shapes.shape[0])]
        tmp_shapes[index] = dictlist

    if index == 'SDI':
        for ifeat in tmp_shapes.iterfeatures():
            # get dsc ts
            dsc_ts = get_dsc(ifeat['properties']['id_station'])
            # TODO: remove this work around to achieve time-series of equal length before operationalization
            dsc_ts = dsc_ts.fillna(method='bfill')
            dsc_ts = dsc_ts.fillna(method='ffill')
            # extend time series to the future
            if dsc_ts.index.max() < dt.datetime(year=2020, month=12, day=31):
                deltat = dt.datetime(year=2020, month=12, day=31) - dsc_ts.index.max()
                dsc_ts_ext = dsc_ts.reindex(index=pd.date_range(dsc_ts.index.min(), '2020-12-31', freq='1d'))
                dsc_ts_ext.loc[dsc_ts.index.max() + dt.timedelta(days=1):'2020-12-31', 'discharge'] = dsc_ts.loc[pd.date_range(dsc_ts.index.min(), dsc_ts.index.min() + deltat - dt.timedelta(days=1)), 'discharge'].values
                dsc_ts = dsc_ts_ext
            # extend time series to the past
            if dsc_ts.index.min() > dt.datetime(year=1981, month=1, day=1):
                deltat = dsc_ts.index.min() - dt.datetime(year=1981, month=1, day=1)
                dsc_ts_ext = dsc_ts.reindex(index=pd.date_range('1981-01-01', dsc_ts.index.max(), freq='1d'))
                dsc_ts_ext.loc['1981-01-01':dsc_ts.index.min() - dt.timedelta(days=1), 'discharge'] = dsc_ts.loc[pd.date_range(dsc_ts.index.max() - (deltat-dt.timedelta(days=1)), dsc_ts.index.max()), 'discharge'].values
                dsc_ts = dsc_ts_ext

            # compute anomalies
            SDI_ts = compute_sdi(dsc_ts)

            for i_date in drange:
                if SDI_ts is not None:
                    if i_date in SDI_ts.index:
                        ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = round(SDI_ts.loc[i_date].values[0], 2)
                    else:
                        ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = None
                else:
                    ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = None

    else:
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
                for ifeat in tmp_shapes.iterfeatures():
                    ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = None
            else:
                r_stats = get_stats(i_path, index, aggunit)
                print(i_date)
                for ifeat, istat in zip(tmp_shapes.iterfeatures(), r_stats):
                    if istat['mean'] is None or istat['count'] == 0:
                        ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = None
                    else:
                        if index == 'CDI':
                            ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = int(istat['mostf'])
                        else:
                            ifeat['properties'][index][i_date.strftime('%Y-%m-%d')] = round(istat['mean'], 3)

    if not update:
        # date range to be deleted after time series export
        delrange = drange[0:335]

    # create time series per nuts region
    if aggunit == 'hydro':
        ts_prefix = 'ID_STATION_'
        shpid = 'id_station'
    else:
        ts_prefix = 'NUTS3_'
        shpid = 'NUTS_ID'
    for ifeat in tmp_shapes.iterfeatures():
        # check if file exists
        ts_outpath = '../json/' + aggunit + '/timeseries/' + ts_prefix + ifeat['properties'][shpid] + '_tmp' + index + '.json'
        ts_list = list()
        for idate in ifeat['properties'][index].items():
            ts_list.append({'date': idate[0], index: idate[1]})
        # write time series
        with open(ts_outpath, 'w') as outfile:
            json.dump(ts_list, outfile)

        if update:
            for ikey in ifeat['properties'][index].keys():
                if dt.datetime.strptime(ikey, '%Y-%m-%d') < (drange[-1] - dt.timedelta(days=30)):
                    del ifeat['properties'][index][ikey]
        else:
            # remove dates from nuts3 maps
            for ddate in delrange:
                del ifeat['properties'][index][ddate.date().strftime('%Y-%m-%d')]

    if aggunit == 'hydro':
        # add point geometry
        if not update:
            tmp_shapes = tmp_shapes.merge(points.rename(columns={'geometry': 'point_geometry'}), on='id_station')
            tmp_shapes['point_geometry'] = tmp_shapes['point_geometry'].to_crs("EPSG:4326").astype(str)
        out_4326 = round_coordinates(tmp_shapes.to_crs("EPSG:4326"))
    else:
        out_4326 = round_coordinates(tmp_shapes.to_crs("EPSG:4326"))

    geo_outpath = '../json/' + aggunit + '/' + index + '-latest.geojson'
    out_4326.to_file(geo_outpath, driver='GeoJSON', encoding='utf-8')


def compare_date(index, index_date, ts_date):
    if index == 'VCI' or index == 'VHI':
        doy8step = (ts_date.dayofyear // 8) * 8 + 1
        if doy8step > 361:
            doy8step = 361
        date8step = dt.datetime.strptime(str(ts_date.year) + '_' + f"{doy8step:03d}", '%Y_%j')
        return index_date == date8step
    else:
        return index_date == ts_date


def parse_index_date(index, path):
    if index == 'SPI-6' or index == 'SPI-12' or index == 'SSPI-10':
        index_date = dt.datetime.strptime(path.name[7:15], '%Y%m%d')
    elif index[0:3] == 'SPI':
        index_date = dt.datetime.strptime(path.name[6:14], '%Y%m%d')
    elif index == 'SPEI-6' or index == 'SPEI-12':
        index_date = dt.datetime.strptime(path.name[8:16], '%Y%m%d')
    elif index[0:3] == 'SPE':
        index_date = dt.datetime.strptime(path.name[7:15], '%Y%m%d')
    elif index[0:3] == 'VHI' or index[0:3] == 'VCI':
        index_date = dt.datetime.strptime(path.name[10:17], '%Y%j')
    elif index[0:3] == 'SMA':
        index_date = dt.datetime.strptime(path.name[13:21], '%Y%m%d')
    elif index[0:3] == 'CDI':
        index_date = dt.datetime.strptime(path.name[4:12], '%Y%m%d')

    return index_date


def get_basepath(index):
    if index[0:3] == 'SPI' or index[0:3] == 'SPE':
        basepath = '/mnt/CEPH_PROJECTS/ADO/ARSO/ERA5_QM_NEW/' + index + '/'
    elif index[0:3] == 'VHI':
        basepath = '/mnt/CEPH_PROJECTS/ADO/VHI/03_results/eusalp_2001_2020/laea_vegmasked/vhi/'
    elif index[0:3] == 'VCI':
        basepath = '/mnt/CEPH_PROJECTS/ADO/VHI/03_results/eusalp_2001_2020/laea_vegmasked/vci/'
    elif index[0:3] == 'SMA':
        basepath = '/mnt/CEPH_PROJECTS/ADO/SM/ERA5/anomalies/'
    elif index[0:3] == 'CDI':
        basepath = '/mnt/CEPH_PROJECTS/ADO/CDI/'
    elif index == 'SSPI-10':
        basepath = '/mnt/CEPH_PROJECTS/ADO/ARSO/SNOWGRID/SSPI10/'
    elif index == 'SDI':
        basepath = ''

    return Path(basepath)


def get_stats(r_path, index, aggunit):
    if aggunit == 'hydro':
        shp_path = "/mnt/CEPH_PROJECTS/ADO/JSON/hydro/catchment_areas_LAEA_simple.geojson"
    else:
        shp_path = "/mnt/CEPH_PROJECTS/ADO/GIS/ERTS89_LAEA/alpinespace_eusalp_NUTS3_simple.shp"
    if index == 'CDI':
        r_stats = zonal_stats(shp_path,
                              r_path,
                              add_stats={'mostf': most_frequent}, nodata=np.nan)
    elif index == 'SMA':
        r_stats = zonal_stats(shp_path,
                              r_path, band=2)
                              #add_stats={'nanmean': mymean}, band=2)
    else:
        r_stats = zonal_stats(shp_path,
                              r_path)
                              #add_stats={'nanmean': mymean})
    return r_stats


def get_dsc(st_id):
    # establish connection do hydro db
    conn = psycopg2.connect(host='10.8.244.31',
                            database='climate_data',
                            user='ado_admin',
                            password='oda347hydro',
                            port=5432)
    cur = conn.cursor()

    # get discharge data
    query = "SELECT date, discharge FROM \"ML_discharge\".mod_disc WHERE id_station = '" + \
            st_id + \
            "' ORDER BY date;"

    tbl_dsc = pd.read_sql_query(query, conn,  index_col='date', parse_dates=['date'])

    # tbl_dsc.set_index('date', inplace=True)

    # TODO: use the real forecast once in operational mode
    # # get forecast
    # query = "SELECT date, 10_d_disch FROM \"ML_discharge\".pred_disch WHERE id_station = '" + \
    #         st_id + \
    #         "'"
    # tbl_forecast = pd.read_sql_query(query, conn, index_col='date', parse_dates=['date'])

    cur.close()
    conn.close()

    if tbl_dsc.size == 0:
        return None
    else:
        return tbl_dsc


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


def compute_sdi(ts):
    from scipy.stats import gamma, norm

    # calulate anomalies
    sdi = ts.copy()
    # for i in range(1, 367, 10):
    for i in np.unique(ts.index.dayofyear):
        # extract all values for day of the year
        doyts_in = ts.loc[ts.index.dayofyear == i]

        # determine the percentage of 0 precipitation values
        p_zero = doyts_in[doyts_in.discharge == 0].shape[0] / doyts_in.shape[0]
        doyts_in = doyts_in[doyts_in != 0]

        params = gamma.fit(doyts_in.dropna())
        gamma_dist = gamma(params[0], loc=params[1], scale=params[2])
        if not (p_zero is None):
            cdf = p_zero + (1 - p_zero) * gamma_dist.cdf(doyts_in)
        else:
            cdf = gamma_dist.cdf(doyts_in)

        norm_ppf = norm.ppf(cdf)
        norm_ppf[np.isinf(norm_ppf)] = np.nan

        sdi[sdi.index.dayofyear == i] = norm_ppf

    return sdi


def merge_time_series(points, index_list, aggunit='hydro', update=False):
    # iterate through nuts zones
    for ifeat in points.iterfeatures():

        ts_list = list()
        # iterate through indices
        if aggunit == 'hydro':
            ts_prefix = 'ID_STATION_'
            iid = ifeat['properties']['id_station']
        else:
            ts_prefix = 'NUTS3_'
            iid = ifeat['properties']['NUTS_ID']
        for iind in index_list:
            ts_inpath = '../json/'+ aggunit + '/timeseries/' + ts_prefix + iid + '_tmp' + iind + '.json'
            tmp = pd.read_json(ts_inpath)
            tmp = tmp.set_index('date')
            ts_list.append(tmp)
            os.remove(ts_inpath)

        # concatenate time series
        ts_out = pd.concat(ts_list, axis=1)
        # reformat for output
        ts_outpath = '../json/'+ aggunit + '/timeseries/' + ts_prefix + iid + '.json'
        ts_out.reset_index(level=0, inplace=True)
        ts_out['date'] = ts_out['date'].astype(str)
        if update:
            # read existing time series
            ts_out_old = pd.read_json(ts_outpath)
            ts_out_old['date'] = ts_out_old['date'].astype(str)
            N_new_rows = pd.to_datetime(ts_out.date).max() - pd.to_datetime(ts_out_old.date).max()
            ts_out = pd.concat([ts_out_old, ts_out.iloc[-N_new_rows.days:, :]], axis=0, ignore_index=True)
            # ts_out = ts_out[ts_out.index.drop_duplicates()]
            # drop rows at the beginning of the file
            ts_out = ts_out.iloc[N_new_rows.days:, :]
            # ts_out['date'] = ts_out['date'].date.astype(str)

        ts_out.to_json(ts_outpath, orient='records', double_precision=3)