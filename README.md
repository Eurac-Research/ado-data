# portal_data

repo for data to be visualized on the ado portal

## Contents

* ./json: 
  * one geojson file per index with values for each NUTS3 region for the past 30 days
  * ./json/timeseries
    * one json file for each NUTS3 region containing the time-series of all indices for the last 365 days
* ./scripts
  * init_geoson.py: initializes all json and geojson files based on the raster data store in the ADO database
* ./visualization
  * a text file for each index (or group of indices) containing visualization parameters