**Author:** [Behrouz Safari](https://behrouzz.github.io/)<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)<br/>

# hypatie
*A python package for querying NASA's JPL HORIZONS API*


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
obs = hp.Observer('sun', t1, t2, step=5)
```

Now you can access the time intervals with *.time* attribute:

```python
print(obs.time)

[datetime.datetime(2021, 3, 20, 8, 0)
 datetime.datetime(2021, 3, 20, 8, 24)
 datetime.datetime(2021, 3, 20, 8, 48)
 datetime.datetime(2021, 3, 20, 9, 12)
 datetime.datetime(2021, 3, 20, 9, 36)
 datetime.datetime(2021, 3, 20, 10, 0)]
```

To acces the position you can use *obs.pos*, *obs.ra*, or *obs.dec*:

```python
print(obs.pos)

[[ 3.59938235e+02 -2.66803120e-02]
 [ 3.59953431e+02 -2.00920520e-02]
 [ 3.59968627e+02 -1.35038600e-02]
 [ 3.59983823e+02 -6.91573600e-03]
 [ 3.59999018e+02 -3.27680000e-04]
 [ 1.42132560e-02  6.26030600e-03]]
```

The first column in the above array is RA and the second column is DEC.

You can request the cartesian positions (x,y,z) of a target with *Vector* class.

```python
vec = hp.Vector('sun', t1, t2, step=5)
```

As with the *Observer* class, there are two attributes *.time* and *.pos* for *Vector* class.
Note the when creating a Vector class, you have *.x*, *.y* and *.z* attributes instead of *ra* and *des*.