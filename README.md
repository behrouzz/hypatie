**Author:** [Behrouz Safari](https://behrouzz.github.io/)<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)<br/>

# hypatie
*A python package for querying NASA's Horizons API*


## Installation

You can install the latest version of *hypatie* from [PyPI](https://pypi.org/project/hypatie/):

    pip install hypatie

The only requirement is *numpy*.


## How to use

Let's get the positions of the sun between two times:

```python
import hypatie as hp

t1 = '2021-03-20 08:00:00'
t2 = '2021-03-20 10:00:00'
```

If you want the apparent RA and DEC of the Sun with respect to Earth (geocentric):

```python
dates, positions = hp.get_position('obs', 'sun', t1, t2, step=100)
```

If you want the cartesian positions (x,y,z) of the Sun with respect to the barycenter of Solar System:

```python
dates, positions = hp.get_position('vec', 'sun', t1, t2, step=100)
```
