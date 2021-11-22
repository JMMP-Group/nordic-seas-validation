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

| NAME | REFERENCE | FUNCTION
|---|---|---|
| EEL | [Loïc Houpert GitHub](https://github.com/lhoupert/analysis_eel_data) | |
| FIM | [10.7489/2036-1](https://doi.org/10.7489/2036-1) | `nsv.Standardizer().fim_1m`<br> `nsv.Standardizer().fim_25m` |
| Kögur | [kogur.whoi.edu](http://kogur.whoi.edu/php/index.php#gridded) | `nsv.Standardizer().kogur` |
| Látrabjarg (climatology) | [Mastropole et al., 2017](https://doi.org/10.1002/2016JC012007) | `nsv.Standardizer().latrabjarg_climatology` |
| Látrabjarg (survey) | [Våge et al., 2011](https://doi.org/10.1038/ngeo1234) | `nsv.Standardizer().latrabjarg_survey` |
| NOAA Arctic Climatology | [10.7289/v5qc01j0](https://doi.org/10.7289/v5qc01j0) | |
| OSNAP | [10.7924/r4z60gf0f](https://doi.org/10.7924/r4z60gf0f) | `nsv.Standardizer().osnap` |
| OVIDE | [10.17882/46446](https://doi.org/10.17882/46446) | `nsv.Standardizer().ovide` |
