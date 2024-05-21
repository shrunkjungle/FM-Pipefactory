# PipeFactory

<p align="center">
  <img src="/resources/img/PipeIcon.png" alt="PipeFactory logo - digiLab" title="PipeFactory logo - digiLab" width="200" height="200">
</p>

> This repo contains the code for the parameteric mesh generation tool for Full Matrix Ltd

## Quickstart

Install the package with default feature flags:

```shell
poetry env use python@3.10
poetry install
```

OR install with parallelism features enabled:

```shell
poetry env use python@3.10
poetry install -E parallel
```

## Testing

Run the tests:

```shell
poetry run pytest
```

With code coverage:

```shell
poetry run pytest --cov=pipefactory --cov-report html
open htmlcov/index.html
```

## Library

Install the PipeFactory library:

```shell
poetry install
```

And then import it into your project:

```python
import numpy as np
import pipefactory as pf

Straight0 = {
    'Name': 'Straight0',
    'length': 350.0,
    'type': 'Straight',
}

Bend1 = {
    'Name': 'Bend1',
    'type': 'Bend',
    'param': {
        'axis' : "left_right",
        'radius': 200.0,
        'angle': 90.0
    }
}

Bend2 = {
    'Name': 'Bend2',
    'type': 'Bend',
    'param': {
        'axis' : "up_down",
        'radius': 150.0,
        'angle': 90.0
    }
}

section_list = [Straight0, Bend2, Bend1]

mesh = pf.Pipe(70.0, 3.0, section_list, ("hex", False), 6., 3)

mesh.export('foo.vtk')
```
