# portal_data

Repository holding the data that is to be visualized on the public portal https://ado.eurac.edu (frontend repo https://github.com/Eurac-Research/ADO)

## Data Usage and Citation
If you are planning to use any data published in this repository, please read [About the Data](https://ado.eurac.edu/md/about-the-data) of the ADO webportal and cite them accordingly. 

## Contents

- ./factsheets/
  - holding the factsheets (.pdfs) describing the drought indices.
- ./html/ 
  - hydrological station reports (.htmls)
- ./json/
  - hydro/
    - metadata/: metadata for each index (abstract, colorscheme, ...)
    - timeseries/: timeseries data dating back to start of each index
    - hydrological basin maps with the last 6 months
  - impacts/
    - static impact maps (.json)
  - nuts/
    - metadata/: metadata for each index (abstract, colorscheme, ...)
    - timeseries/: timeseries data dating back to start of each index
    - nuts2/3 maps with the last 6 months
  - vulnerabilities/
    - static vulnerability maps (.json)
- ./markdown/navigation-items/
  - description about the data, about the project, imprint, ...
- ./visualization/
  - a text file for each index (or group of indices) containing visualization parameters

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
