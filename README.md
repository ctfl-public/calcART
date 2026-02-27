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


## References
1. Yassin, A. H., and Poovathingal, S. J., “Characterization of Directional and Anisotropic Scattering Dependency of Emissivity for Fibrous Heat Shields under Non-Isothermal Conditions,” International Journal of Heat and Mass Transfer, Vol. 220, 2024, p. 124859. https://doi.org/10.1016/j.ijheatmasstransfer.2023.124859
2. Yassin, A. H., and Poovathingal, S. J., “Development of a physics-based radiative model for anisotropic scattering medium,” 11th International Symposium on Radiative Transfer, Kusadasi (Turkey), 2025. https://doi.org/10.1615/RAD-25.430
3. Banerjee, A., Yassin, A. H., Davuluri, R. S. C., Martin, A., and Poovathingal, S. J., “Radiative Coefficients and Their Influence on In-Depth Heating of Porous Ablators,” AIAA Journal, Vol. 60, No. 12, 2022, pp. 6520–6535. https://doi.org/10.2514/1.J061953

