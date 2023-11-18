# nordic-seas-validation
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10149505.svg)](https://doi.org/10.5281/zenodo.10149505)
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
| EEL | [Loïc Houpert GitHub](https://github.com/lhoupert/analysis_eel_data) | `nsv.Standardizer().eel` |
| FIM | [10.7489/2036-1](https://doi.org/10.7489/2036-1) | `nsv.Standardizer().fim(1)`<br> `nsv.Standardizer().fim(25)` |
| Kögur | [kogur.whoi.edu](http://kogur.whoi.edu/php/index.php#gridded) | `nsv.Standardizer().kogur` |
| Látrabjarg (climatology) | [Mastropole et al., 2017](https://doi.org/10.1002/2016JC012007) | `nsv.Standardizer().latrabjarg_climatology` |
| Látrabjarg (survey) | [Våge et al., 2011](https://doi.org/10.1038/ngeo1234) | `nsv.Standardizer().latrabjarg_survey` |
| NOAA Arctic Climatology | [10.7289/v5qc01j0](https://doi.org/10.7289/v5qc01j0) | |
| OSNAP | [10.7924/r4z60gf0f](https://doi.org/10.7924/r4z60gf0f) | `nsv.Standardizer().osnap` |
| OVIDE | [10.17882/46446](https://doi.org/10.17882/46446) | `nsv.Standardizer().ovide` |
| HO2000 | [Hansen & Osterhus 2000](https://doi.org/10.1016/S0079-6611(99)00052-X) | `nsv.Standardizer().ho2000` |
| M82/1 | [10.1594/PANGAEA.890362](https://doi.pangaea.de/10.1594/PANGAEA.890362) | `nsv.Standardizer().m82_1(1)` <br> ...  <br> `nsv.Standardizer().m82_1(9)`|
| KN203-2 | [10.1594/PANGAEA.919251-2](https://doi.pangaea.de/10.1594/PANGAEA.919251) | `nsv.Standardizer().kn203_2("A")` <br> ...  <br> `nsv.Standardizer().kn203_2("E")` |
| POS503 | [10.1594/PANGAEA.890699](https://doi.org/10.1594/PANGAEA.890699) | `nsv.Standardizer().pos503(1)` <br> ...  <br> `nsv.Standardizer().pos503(5)` |
