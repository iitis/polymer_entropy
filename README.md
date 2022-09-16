[![Pylint](https://github.com/iitis/polymer_entropy/actions/workflows/pylint.yml/badge.svg)](https://github.com/iitis/polymer_entropy/actions/workflows/pylint.yml)

# polymer_entropy
computation of conformational entropy from polymer simulations

## Using software

1. create virtual environment
```
$ python3 -m venv venv
```

2. use virtual environment
```
$ source venv/bin/activate
```

3. install project dependencies
```
$ pip3 install -r requirements.txt
```

4. create directory for plots
```
$ mkdir pics
```

5. Example call:
```
$ python3 src/plot.py --datafolder data/HA/ --complex 'Albumin+HA'
```

tested with
```
$ python3 -V
Python 3.8.10
$ uname -osrmp
Linux 5.4.0-90-generic x86_64 x86_64 GNU/Linux
```

## Citing this work

While using this code please cite: P. Sionkowski, P. Bełdowski, N. Kruszewska, P. Weber, B. Marciniak and
K. Domino, 'Effect of ion and binding site on the conformation of chosen glycosaminoglycans at the albumin surface',	Entropy 2022, 24, 811.  https://www.mdpi.com/1099-4300/24/6/811/pdf  DOI: 	10.3390/e24060811 

Input data were computed at the Academic Computer Centre in Gdańsk.
