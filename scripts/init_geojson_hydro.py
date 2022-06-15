import geopandas as gpd
import pandas as pd
import datetime as dt
import tools
from joblib import Parallel, delayed


def main():
    index_list = ['SDI', 'SSPI-10',
                  'SPI-1', 'SPI-2',
                  'SPI-3', 'SPI-6',
                  'SPI-12', 'SPEI-1',
                  'SPEI-2', 'SPEI-3',
                  'SPEI-6', 'SPEI-12',
                  'VHI', 'VCI',
                  'SMA', 'CDI']

    # import shapefile
    shapes = gpd.read_file('/mnt/CEPH_PROJECTS/ADO/JSON/hydro/catchment_areas_LAEA_simple.geojson',
                           encoding='utf-8')
    points = gpd.read_file('/mnt/CEPH_PROJECTS/ADO/JSON/hydro/gauging_stations_LAEA.geojson',
                           encoding='utf-8')
    station_list = ['ADO_DSC_AT12_0280',
                    'ADO_DSC_AT31_0206',
                    'ADO_DSC_AT31_0254',
                    'ADO_DSC_CH03_0075',
                    'ADO_DSC_CH04_0011',
                    'ADO_DSC_CH05_0201',
                    'ADO_DSC_CH07_0006',
                    'ADO_DSC_CH07_0100',
                    'ADO_DSC_CH07_0147',
                    'ADO_DSC_FRK2_0041',
                    'ADO_DSC_FRK2_0042',
                    'ADO_DSC_ITC1_0020',
                    'ADO_DSC_ITC1_0037',
                    'ADO_DSC_ITC1_0072',
                    'ADO_DSC_ITH1_0012',
                    'ADO_DSC_ITH2_0035',
                    'ADO_DSC_ITH5_0006',
                    'ADO_DSC_SI03_0033',
                    'ADO_DSC_SI03_0148']
    shapes = shapes[shapes['id_station'].isin(station_list)]
    points = points[points['id_station'].isin(station_list)]

    # set date range
    today = dt.datetime.strptime('2020-03-23',
                                 '%Y-%m-%d')  # TODO: modify to use actual date of the day the script is run
    drange = pd.date_range(today - dt.timedelta(days=364), today, freq='1D')

    def call_compute_stats(index):
        tools.compute_stats(drange, shapes, points=points, index=index, aggunit='hydro')
        tools.update_metadata_drange(index, drange)

    Parallel(n_jobs=-1)(delayed(call_compute_stats)(i) for i in index_list)

    tools.merge_time_series(points, index_list, aggunit='hydro')


if __name__ == "__main__":
    main()