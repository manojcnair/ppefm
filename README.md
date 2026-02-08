# PPEFM (Prompt Penetration Electric Field Model)

![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey)
![MATLAB](https://img.shields.io/badge/MATLAB-R2018b%2B-blue)

This repository provides MATLAB implementations of the Prompt Penetration Electric Field Model (PPEFM) and an example script that downloads OMNI 5‑minute solar wind data and computes equatorial electric fields in the time domain.

**Real‑time implementation (web tool):**
```
https://geomag.colorado.edu/index.php/online-calculators/real-time-model-ionospheric-electric-fields
```

**Contents**
- `ppefm_time_domain.m` — time‑domain PPEFM implementation using 2008 transfer‑function coefficients and LT response.
- `ppefm_omni_example.m` — minimal, clean example script that downloads OMNI 5‑minute GSM By/Bz + Vsw and runs PPEFM.
- `vdrift-model_jgr0499.f` — original Scherliess & Fejer (1999) vertical drift model (FORTRAN).

## Quick Start (MATLAB)
```matlab
% Make sure ppefm_time_domain.m is on your path
addpath(pwd);

% Run the example script
ppefm_omni_example
```

## Requirements
- MATLAB R2018b or newer (uses `datetime`, `webread`, `fillmissing`)
- No additional toolboxes required

## Outputs & Units
The model outputs the equatorial eastward electric field (EEF) and intermediate IEF components in **mV/m**. An example transfer‑function response is shown below.

![Transfer‑function response example](images/tf_response_v2.png)
*Transfer‑function response to synthetic IEF Ey input. Left: input signal. Right: zonal electric‑field output. Adapted from Manoj et al. (2008).*

## Climatology Model (Scherliess & Fejer 1999)
The climatological (quiet‑time) component is based on the Scherliess & Fejer (1999) model of equatorial vertical plasma drift. This repository does **not** implement the climatology in MATLAB; the real‑time modeling tool includes it, and the original FORTRAN code is provided here for reference. Users can run the FORTRAN model with local time, longitude, day‑of‑year, and F10.7, then add the resulting quiet‑time EEF to the PPEFM output if they want **quiet + prompt‑penetration** fields.

For the climatological model:
- It evaluates vertical drift as a function of local time, longitude, day‑of‑year, and F10.7.
- Quiet‑time equatorial ionospheric eastward electric field (climatology of EEF) may be obtained by multiplying the **vertical drift** output from the FORTRAN code by the local magnetic‑field strength (e.g., from the **International Geomagnetic Reference Field, IGRF** or the **World Magnetic Model, WMM**): **E (mV/m) = V (m/s) × B (nT) / 10⁶**.

The original FORTRAN is included as `vdrift-model_jgr0499.f` for reference.

## Citations
Please cite the following works when using this model:

**PPEFM / Prompt Penetration**
- Manoj, C., S. Maus, and P. Alken (2013), *Long-period prompt-penetration electric fields derived from CHAMP satellite magnetic measurements*, J. Geophys. Res. Space Physics, 118, 5919–5930, doi:10.1002/jgra.50511.
- Manoj, C., and S. Maus (2012), *A real-time forecast service for the ionospheric equatorial zonal electric field*, Space Weather, 10, S09002, doi:10.1029/2012SW000825.
- Manoj, C., S. Maus, H. Lühr, and P. Alken (2008), *Penetration characteristics of the interplanetary electric field to the daytime equatorial ionosphere*, J. Geophys. Res., 113, A12310, doi:10.1029/2008JA013381.

**Climatology (Quiet‑time) Model**
- Scherliess, L., and B. G. Fejer (1999), *Radar and satellite global equatorial F region vertical drift model*, J. Geophys. Res., 104, 6829–6842.

## Notes
- OMNI data are already propagated to the bow‑shock nose. The PPEFM implementation applies the additional ~17‑minute delay from bow shock to equatorial ionosphere in the LT response phase.
- Input cadence is expected to be 5 minutes (OMNI 5‑min resolution).
