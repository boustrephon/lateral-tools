
# pySeismic
Seismic tools for Python (under development)

The general scope of this repository include:

* Extraction of seismic ground motion data from files (ONLY for cases where other options do not exist, see below)
* Code calculations
* Seismic plots
* Sample notebooks

Note that the library `obspy` contains readers for many file formats. The documentation also lists many software tools that are based on obspy and are worth looking into before developing anything for yourself (see links below and [the list of formats that `obspy` can handle](https://docs.obspy.org/packages/index.html)).

Note that the ASC file format for the 2023 Turkey files in the resources folder do not appear to be covered by the `obspy` read function.

## References

The following libraries may be of use.

* [`eqsig` on PyPI](https://pypi.org/project/eqsig/) ([documentation](https://eqsig.readthedocs.io/en/latest/eqsig.html))
* [`openseespy` on PyPI](https://pypi.org/project/openseespy/) ([documentation](https://openseespydoc.readthedocs.io/en/latest/index.html))
* [`obspy` on PyPI](https://pypi.org/project/obspy/) ([documentation](https://docs.obspy.org/), [wiki](https://github.com/obspy/obspy/wiki/)) - a Python framework for processing seismological data. It provides parsers for common file formats, clients to access data centers and seismological signal processing routines which allow the manipulation of seismological time series. Installation using conda is recommended - `conda install -c conda-forge obspy` - install into a `conda-forge`-based environment.
* [`pyweed` on PyPI](https://pypi.org/project/pyweed/) ([documentation](https://iris-edu.github.io/pyweed/)) - a GUI-based application for retrieving event-based seismic data based on ObsPy and PyQt.
* [`hvsrpy` on PyPI](https://pypi.org/project/hvsrpy/) ([documentation](https://hvsrpy.readthedocs.io/en/latest/?badge=latest)) - horizontal to vertical spectral ratio processing
* [`openseespy` on PyPI]() ([documentation]()) 



