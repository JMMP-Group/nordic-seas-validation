<a id="__init__"></a>

# \_\_init\_\_

<a id="utils"></a>

# utils

<a id="standardizer"></a>

# standardizer

<a id="standardizer.Standardizer"></a>

## Standardizer Objects

```python
@dataclass
class Standardizer()
```

Standardize raw data

<a id="standardizer.Standardizer.kogur"></a>

#### kogur

```python
@property
def kogur() -> Dataset
```

Standardized Kögur dataset

<a id="standardizer.Standardizer.latrabjarg_climatology"></a>

#### latrabjarg\_climatology

```python
@property
def latrabjarg_climatology() -> Dataset
```

Standardized Látrabjarg climatology dataset

<a id="standardizer.Standardizer.latrabjarg_survey"></a>

#### latrabjarg\_survey

```python
@property
def latrabjarg_survey() -> Dataset
```

Standardized Látrabjarg survey dataset

<a id="standardizer.Standardizer.fim_1m"></a>

#### fim\_1m

```python
@property
def fim_1m() -> Dataset
```

Standardized FIM in 1m depth bins

<a id="standardizer.Standardizer.fim_25m"></a>

#### fim\_25m

```python
@property
def fim_25m() -> Dataset
```

Standardized FIM in 25m depth bins

<a id="standardizer.Standardizer.osnap"></a>

#### osnap

```python
@property
def osnap() -> Dataset
```

Standardized OSNAP dataset

<a id="standardizer.Standardizer.ovide"></a>

#### ovide

```python
@property
def ovide() -> Dataset
```

Standardized OVIDE dataset

<a id="standardizer.Standardizer.eel"></a>

#### eel

```python
@property
def eel()
```

Standardized EEL dataset

<a id="section_finder"></a>

# section\_finder

<a id="section_finder.SectionFinder"></a>

## SectionFinder Objects

```python
@dataclass
class SectionFinder()
```

Parameters
----------
ds_domain: Dataset
    domain_cfg dataset

<a id="section_finder.SectionFinder.grids"></a>

#### grids

```python
@property
def grids()
```

Dictionary mapping each grid to a dataset with its coordinates

<a id="section_finder.SectionFinder.nearest_neighbor"></a>

#### nearest\_neighbor

```python
def nearest_neighbor(lons, lats, grid)
```

Given the coordinates defining a section, find the nearest points
on a model grid.

**Arguments**:

- `lons` _1D array-like_ - Longitudes defining a section
- `lats` _1D array-like_ - Latitudes defining a section
- `grid` _string_ - Model grid `{"u", "v", "t", "f"}`
  

**Returns**:

  Dataset

<a id="section_finder.SectionFinder.zigzag_section"></a>

#### zigzag\_section

```python
def zigzag_section(lons, lats, grid)
```

Given the coordinates defining a section, find the correspoinding zigzag section
on a model grid.

**Arguments**:

- `lons` _1D array-like_ - Longitudes defining a section
- `lats` _1D array-like_ - Latitudes defining a section
- `grid` _string_ - Model grid `{"u", "v", "t", "f"}`
  

**Returns**:

  Dataset

<a id="section_finder.SectionFinder.velocity_points_along_zigzag_section"></a>

#### velocity\_points\_along\_zigzag\_section

```python
def velocity_points_along_zigzag_section(lons, lats)
```

Given the coordinates defining a section, find the corrisponding velocity points
along a zigzag section (f-grid). Useful to compute accurate volume fluxes.

**Arguments**:

- `lons` _1D array-like_ - Longitudes defining a section
- `lats` _1D array-like_ - Latitudes defining a section
  

**Returns**:

- `dict` - Dictionary mapping u/v grids to their coordinates and indexes.

