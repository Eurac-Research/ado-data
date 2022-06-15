import geopandas as gpd
import pandas as pd
import datetime as dt
import tools
from joblib import Parallel, delayed


def main():
    index_list = ['SPI-1', 'SPI-2',
                  'SPI-3', 'SPI-6',
                  'SPI-12', 'SPEI-1',
                  'SPEI-2', 'SPEI-3',
                  'SPEI-6', 'SPEI-12',
                  'SSPI-10',
                  'SMA', 'VHI', 'VCI', 'CDI']

    # import shapefile
    nuts3 = gpd.read_file('/mnt/CEPH_PROJECTS/ADO/GIS/ERTS89_LAEA/alpinespace_eusalp_NUTS3_simple.shp',
                          encoding='utf-8')

    # set date range
    today = dt.datetime.strptime('2022-03-23',
                                 '%Y-%m-%d')  # TODO: modify to use actual date of the day the script is run
    drange = pd.date_range(today - dt.timedelta(days=364), today, freq='1D')

    today = dt.datetime.strptime('2020-08-01',
                                 '%Y-%m-%d')
    drange2 = pd.date_range(today - dt.timedelta(days=364), today, freq='1D')

    def call_compute_stats(index):
        usedrange = drange2 if index in ['VHI', 'VCI', 'CDI'] else drange
        tools.compute_stats(usedrange, nuts3, index=index, aggunit='nuts')
        # tools.update_metadata_drange(index, usedrange)

    Parallel(n_jobs=-1)(delayed(call_compute_stats)(i) for i in index_list)

    tools.merge_time_series(nuts3, index_list, aggunit='nuts')


if __name__ == "__main__":
    main()