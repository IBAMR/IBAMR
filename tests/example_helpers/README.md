The four cases here cover the main behavior of the compiled example helpers.
They are built explicitly with `tests-example_helpers`, separately from the
ordinary library test target.

- One helper case per dimension checks the multiphase ex3 sphere and velocity
  initializers, level-set resets, Heaviside synchronization, and liquid-fraction
  force masking. Shifted coordinates and spatially varying fields expose
  indexing and averaging errors.
- One Enthalpy and one Allen-Cahn application case run the shared
  `CoupledApplication` lifecycle, including its timestep loop, material
  callbacks, mass diagnostic, and restart output. A periodic Fourier mode has
  a known discrete Crank-Nicolson decay; the other direction uses physical
  zero-flux boundaries. Each case runs four steps, restarts at step two, and
  checks the resumed field and diagnostic results.

Input files specify the numerical configuration.

From a Debug build, use `attest -R '^example_helpers/'`. The discovery list
`expected-tests.txt` is checked by CI so missing executables or fixture links
cannot silently omit this coverage.
