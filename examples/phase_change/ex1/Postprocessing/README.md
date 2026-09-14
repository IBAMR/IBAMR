# Exporting the liquid-solid interface from VisIt

1. Open the simulation output in VisIt and add a Contour plot of
   `liquid_fraction` with a single contour value of `0.5`.
2. Select **File > Export database**, then choose the output directory and file
   name.
3. Add `liquid_fraction` to the exported variables and select the **XYZ** format.
4. Select **Export all time states** to export the full time series. Otherwise,
   VisIt exports only the current time state.

The XYZ files contain two header lines followed by rows with a label, the x, y,
and z coordinates, and the liquid fraction. For this single-contour export,
remove the headers and label before reading the numeric data in MATLAB:

```sh
awk 'NR > 2 { print $2, $3, $4, $5 }' interface.xyz > interface.dat
```

VisIt interpolates the liquid fraction to locate the `0.5` contour. If a time
state has no such contour, record the missing state and omit its entry from a
copy of `dumps.visit` used for the export.

The accompanying `stefan_analytical_solution.m` computes the analytical interface
position. Its material parameters and time interval must correspond to the
simulation being compared.
