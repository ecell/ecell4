# ecell4

This provides pure Python libraries (Miscellaneous utility functions) for [E-Cell System version 4](https://github.com/ecell/ecell4-base).

Try online
----------

You can try this package online with Google Colaboratory.
Please refer to the https://github.com/ecell/ecell4_docs

Quick Start
-----------

Here are 3 extremely simple examples, See https://github.com/ecell/ecell4_docs or http://ecell4.readthedocs.org for more details.

```
(base) root@5e5b25d2363d:~# python
Python 3.7.3 (default, Mar 27 2019, 22:11:17)
[GCC 7.3.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from ecell4_base.core import *
>>> sp = Species("A.B.C")
>>> sp.serial()
'A.B.C'
>>>
```

### Binding and unbinding reactions

Run this with Jupyter Notebook

```python
%matplotlib inline
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

run_simulation(10, {'A': 60, 'B': 60})
```

![png](./samples/output_7_0.png)

### Diffusion on a spherical surface

Run this with Jupyter Notebook

```python
%matplotlib inline
from ecell4_base.core import *
from ecell4 import *

with species_attributes():
    M | {'dimension': 2}
    A | {'D': 1.0, 'location': 'M'}

surface = Sphere(ones() * 0.5, 0.49).surface()
obs = FixedIntervalTrajectoryObserver(1e-4)
run_simulation(
    0.4, y0={'A': 10}, structures={'M': surface},
    solver='spatiocyte', observers=obs, return_type=None)

viz.plot_trajectory(obs, interactive=False)
```

![png](./samples/hairball.png)

Installation
------------

Please refer to https://github.com/ecell/ecell4_base#installation
