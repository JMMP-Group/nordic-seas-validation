# nordic-seas-validation

## Quick-start

```shell
git clone git@github.com:JMMP-Group/nordic-seas-validation.git
cd nordic-seas-validation
conda env create -f environment.yml
conda activate nordic-seas-validation
pip install -e .
python -c "import nsv"
```

## Raw data
[JMMP GWS](https://gws-access.jasmin.ac.uk/public/jmmp/NORVAL/)

| NAME | REFERENCE |
|---|---|
| EEL | [Loïc Houpert GitHub](https://github.com/lhoupert/analysis_eel_data)
| FIM | [10.7489/2036-1](https://doi.org/10.7489/2036-1) |
| Kögur | [kogur.whoi.edu](http://kogur.whoi.edu/php/index.php#gridded) |
| Látrabjarg (hydrography) | [Mastropole et al., 2017](https://doi.org/10.1002/2016JC012007) |
| Látrabjarg (velocity) | [Våge et al., 2011](https://doi.org/10.1038/ngeo1234) |
| NOAA Arctic Climatology | [10.7289/v5qc01j0](https://doi.org/10.7289/v5qc01j0) |
| OSNAP | [10.7924/r4z60gf0f](https://doi.org/10.7924/r4z60gf0f) |
| OVIDE | [10.17882/46446](https://doi.org/10.17882/46446)


## Standardize raw data

```python
from nsv import Standardizer
ds = Standardizer(raw_data_path=None).kogur
```
