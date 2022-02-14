<a id="nsv"></a>

# nsv

<a id="nsv.utils"></a>

# nsv.utils

<a id="nsv.standardizer"></a>

# nsv.standardizer

<a id="nsv.standardizer.Standardizer"></a>

## Standardizer Objects

```python
@dataclass
class Standardizer()
```

Standardize raw data

<a id="nsv.standardizer.Standardizer.kogur"></a>

#### kogur

```python
@property
@final_cleanup_before_returning
def kogur() -> Dataset
```

Standardized Kögur dataset

<a id="nsv.standardizer.Standardizer.latrabjarg_climatology"></a>

#### latrabjarg\_climatology

```python
@property
@final_cleanup_before_returning
def latrabjarg_climatology() -> Dataset
```

Standardized Látrabjarg climatology dataset

<a id="nsv.standardizer.Standardizer.latrabjarg_survey"></a>

#### latrabjarg\_survey

```python
@property
@final_cleanup_before_returning
def latrabjarg_survey() -> Dataset
```

Standardized Látrabjarg survey dataset

<a id="nsv.standardizer.Standardizer.fim"></a>

#### fim

```python
@final_cleanup_before_returning
def fim(resolution: int) -> Dataset
```

Standardized FIM dataset

**Arguments**:

- `sec` _int_ - resolution {1, 25}
  

**Returns**:

- `Dataset` - Standardized dataset

<a id="nsv.standardizer.Standardizer.osnap"></a>

#### osnap

```python
@property
@final_cleanup_before_returning
def osnap() -> Dataset
```

Standardized OSNAP dataset

<a id="nsv.standardizer.Standardizer.ovide"></a>

#### ovide

```python
@property
@final_cleanup_before_returning
def ovide() -> Dataset
```

Standardized OVIDE dataset

<a id="nsv.standardizer.Standardizer.ho2000"></a>

#### ho2000

```python
@property
@final_cleanup_before_returning
def ho2000() -> Dataset
```

Standardized Hansen & Osterhus 2000 dataset

<a id="nsv.standardizer.Standardizer.m82_1"></a>

#### m82\_1

```python
@final_cleanup_before_returning
def m82_1(sec_id) -> Dataset
```

Standardized Meteor cruise M82/1

**Arguments**:

- `sec` _int_ - Section ID {1, 2, 3, 4, 5, 6, 7, 8, 9}
  

**Returns**:

- `Dataset` - Standardized dataset

<a id="nsv.standardizer.Standardizer.eel"></a>

#### eel

```python
@property
@final_cleanup_before_returning
def eel() -> Dataset
```

Standardized EEL dataset

<a id="nsv.standardizer.Standardizer.kn203_2"></a>

#### kn203\_2

```python
@final_cleanup_before_returning
def kn203_2(sec_id: str) -> Dataset
```

Standardized Knorr cruise KN203-2

**Arguments**:

- `sec_id` _str_ - Section ID {"A", "B", "C", "D", "E"}
  

**Returns**:

- `Dataset` - Standardized dataset

<a id="nsv.utils"></a>

# nsv.utils

<a id="nsv._version"></a>

# nsv.\_version

<a id="nsv.section_finder"></a>

# nsv.section\_finder

<a id="nsv.section_finder.SectionFinder"></a>

## SectionFinder Objects

```python
@dataclass
class SectionFinder()
```

**Arguments**:

- `ds_domain` _Dataset_ - domain_cfg dataset

<a id="nsv.section_finder.SectionFinder.grids"></a>

#### grids

```python
@property
def grids() -> dict
```

Dictionary mapping each grid to a dataset with its coordinates

<a id="nsv.section_finder.SectionFinder.nearest_neighbor"></a>

#### nearest\_neighbor

```python
def nearest_neighbor(lons, lats, grid: str) -> Dataset
```

Given the coordinates defining a section, find the nearest points
on a model grid.

**Arguments**:

- `lons` _1D array-like_ - Longitudes defining a section
- `lats` _1D array-like_ - Latitudes defining a section
- `grid` _string_ - Model grid `{"u", "v", "t", "f"}`
  

**Returns**:

- `Dataset` - Dataset with model coordinates and indexes

<a id="nsv.section_finder.SectionFinder.zigzag_section"></a>

#### zigzag\_section

```python
def zigzag_section(lons, lats, grid: str) -> Dataset
```

Given the coordinates defining a section, find the correspoinding zigzag section
on a model grid.

**Arguments**:

- `lons` _1D array-like_ - Longitudes defining a section
- `lats` _1D array-like_ - Latitudes defining a section
- `grid` _string_ - Model grid `{"u", "v", "t", "f"}`
  

**Returns**:

- `Dataset` - Dataset with model coordinates and indexes

<a id="nsv.section_finder.SectionFinder.velocity_points_along_zigzag_section"></a>

#### velocity\_points\_along\_zigzag\_section

```python
def velocity_points_along_zigzag_section(lons, lats) -> dict
```

Given the coordinates defining a section, find the corrisponding velocity points
along a zigzag section (f-grid). Useful to compute accurate volume fluxes.

**Arguments**:

- `lons` _1D array-like_ - Longitudes defining a section
- `lats` _1D array-like_ - Latitudes defining a section
  

**Returns**:

- `dict` - Dictionary mapping u/v grids to their coordinates and indexes.

