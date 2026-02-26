# calcART: calculations of Anisotopic Radiative Transport

## Installation
Add calcART to your Python environment
```bash
export PYTHONPATH="/path/to/calcART:$PYTHONPATH"
```

## Usage
```python
import ART as art
# calculate emissivity at scattering albedo of 0.6
epsilon = art.calc_emissivity_bezier(omega=0.6)
``` 
More detailed usage examples will be provided in the future. For now, please refer to the source code for available functions and their parameters.
