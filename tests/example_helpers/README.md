The four cases here cover the main behavior of the compiled example helpers.
They are built explicitly with `tests-example_helpers`, separately from the
ordinary library test target.

- One combined helper case per dimension checks geometry initialization and
  reset, Heaviside synchronization, and liquid-fraction force masking. Shifted
  coordinates and spatially varying fields expose indexing and averaging
  errors without separate input files for each helper.
- One Enthalpy and one Allen-Cahn application case run the real shared
  `CoupledApplication` lifecycle, including its timestep loop, material
  callbacks, mass diagnostic, and restart output. A periodic Fourier mode has
  a known discrete Crank-Nicolson decay; the other direction uses physical
  zero-flux boundaries. Each case runs four steps, restarts at step two, and
  checks the resumed field and diagnostic results.

Inputs contain the numerical configuration. Tests do not construct solver or
hierarchy databases in C++. The application cases demonstrate representative
shared execution; they are not a parameter sweep or a replacement for physical
validation of every example.

From a Debug build, use `attest -R '^example_helpers/'`. The discovery list
`expected-tests.txt` is checked by CI so missing executables or fixture links
cannot silently omit this coverage.
