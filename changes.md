# Applied Changes Overview
The first RK3+2nd order IBM implementation wasn't working correctly I made various changes.

## RK3
- RK3 coefficient wasnt put into the calculation for the divergence and correction schemes. Thats has been changed.
- Pressure gradient was inside the RHS calculation this was a problematic that has been removed and pressure gradient added at the intermediate velocity calculation.
- RK3 coefficients calculation moved inside RK3 steps. Since the size of the dt changes. 

## IBM 
- IBM arrays were shifted that has been corrected.