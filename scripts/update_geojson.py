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

    def call_compute_stats(index):
        # determine the last date in the index geojson
        tmp_geojson = gpd.read_file('../json/nuts/' + index + '-latest.geojson')
        tmp_series = pd.DataFrame.from_dict(tmp_geojson[index][0], orient='index', columns=['droughtindex'])
        maxdate = tmp_series.index.max()
        drange = pd.date_range(maxdate + dt.timedelta(days=1), today, freq='1D')
        total_drange = pd.date_range(today - dt.timedelta(days=364), today, freq='1D')

        tools.compute_stats(drange, tmp_geojson, index=index, aggunit='nuts', update=True)
        tools.update_metadata_drange(index, total_drange)

    Parallel(n_jobs=-1)(delayed(call_compute_stats)(i) for i in index_list)

    tools.merge_time_series(nuts3, index_list, aggunit='nuts', update=True)


if __name__ == "__main__":
    main()