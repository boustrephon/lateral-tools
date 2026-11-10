# lateral-tools

Utilities for seismic and lateral load engineering.

> **Note:** the Git repository directory is named `pySeismic`, but the
> installable Python package (and its import name) is `lateral_tools`
> (PyPI distribution `lateral-tools`).

## Overview

`lateral-tools` is a Python library of utilities for structural engineers
working with earthquake and wind loading. It provides code-based design
calculations, seismic hazard and storey analysis, ground-motion data
importing, and plotting helpers.

## Package contents

### Core modules

- `math_utils` — interpolation (`interp`, `list_interp`), type conversion
  (`as_float`, `as_numeric`), and unit scaling (`scale_to_SI_units`).
- `risk_utils` — seismic hazard recurrence (`exceedance2return`,
  `return2exceedance`), PGA to Ss/S1 (`PGA2S475`, `S_calc`), and
  remaining-life PGA estimates.
- `analysis` — storey elevations, force distribution, shear/moment
  envelopes, and member shear interpolation.
- `print_utils` — table formatting, storey data display, CSV import, and
  rebar detailing.

### Design codes (`lateral_tools.codes`)

- **ASCE 7** — `ASCE7_05`, `ASCE7_10`, `ASCE7_16`, `ASCE7_22`
- **ASCE 41** — `ASCE41_06`, `ASCE41_13`, `ASCE41_17`, `tier1_checks`,
  `asce41_C1_assessment`
- **Eurocodes** — `EC0`, `EC1_1_4`, `EC8`
- **Chinese codes** — `GB50011`, `GB50023_2009`
- **Wind codes** — `HKWind93`, `HKWind2004`, `HKWind2019`, `MS1553_02`,
  `RSAEEP_2008`
- **Other** — `ACI318_08M`, `Thai`

### Plotting (`lateral_tools.plotting`)

- `spectra` — code response-spectrum plots
- `tripartite` — tripartite (log-log) spectral plots
- `storey` — storey shear and moment diagrams
- `eq` — earthquake time-history and response-spectrum plots

### Importing (`lateral_tools.importing`)

- `importing` — reads ASC-format ground-motion files (e.g. the 2023 Turkey
  earthquake records in `resources/`).

## Installation

Requires Python >= 3.9.

```bash
pip install -e .         # install the package (editable)
pip install -e ".[dev]"  # include development tooling (pytest, ruff, black)
```

## Quick start

```python
from lateral_tools.codes.ASCE7_16 import design_response_spectrum
from lateral_tools.plotting.spectra import plot_seismic_spectrum
```

## Testing

```bash
pytest
```

## Notes on external libraries

The library [`obspy`](https://pypi.org/project/obspy/) contains readers for
many file formats, and its documentation lists many software tools built on
it that are worth reviewing before developing your own (see
[the list of formats that `obspy` can handle](https://docs.obspy.org/packages/index.html)).

Note that the ASC file format used by the 2023 Turkey files in the
`resources/` folder does not appear to be covered by the `obspy` read
function.

## References

The following libraries may be of use.

- [`eqsig` on PyPI](https://pypi.org/project/eqsig/) ([documentation](https://eqsig.readthedocs.io/en/latest/eqsig.html))
- [`openseespy` on PyPI](https://pypi.org/project/openseespy/) ([documentation](https://openseespydoc.readthedocs.io/en/latest/index.html))
- [`obspy` on PyPI](https://pypi.org/project/obspy/) ([documentation](https://docs.obspy.org/), [wiki](https://github.com/obspy/obspy/wiki/)) - a Python framework for processing seismological data. It provides parsers for common file formats, clients to access data centers and seismological signal processing routines which allow the manipulation of seismological time series. Installation using conda is recommended - `conda install -c conda-forge obspy` - install into a `conda-forge`-based environment.
- [`pyweed` on PyPI](https://pypi.org/project/pyweed/) ([documentation](https://iris-edu.github.io/pyweed/)) - a GUI-based application for retrieving event-based seismic data based on ObsPy and PyQt.
- [`hvsrpy` on PyPI](https://pypi.org/project/hvsrpy/) ([documentation](https://hvsrpy.readthedocs.io/en/latest/?badge=latest)) - horizontal to vertical spectral ratio processing

## TODO

- ASCE41-17 Tier 1 calculations are not complete and need work
- Tests need to be added
- Calculations need QC review before use on projects
