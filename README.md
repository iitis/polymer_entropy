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
