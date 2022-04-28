# polymer_entropy
computation of conformational entropy from polymer simulations

Running:

1. create virtual environment
$ python3 -m venv venv

2. use virtual environment
$ source venv/bin/activate

3. install project dependencies
$ pip3 install -r requirements.txt

4. create directory for plots
$ mkdir pics

5. Example call:
$ python3 src/plot.py --datafolder data/HA/ --complex 'Albumin+HA'

tested with
```
$ python3 -V
Python 3.8.10
```
