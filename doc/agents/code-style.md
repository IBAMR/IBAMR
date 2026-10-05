# Code style and organization

- `include/ibamr/` and `src/` contain IBAMR interfaces and implementations;
  `ibtk/include/ibtk/` and `ibtk/src/` contain the reusable toolkit. Inline
  implementations belong in the corresponding `private/` headers. Cartesian
  Fortran kernels often require `m4` preprocessing; edit their sources, not
  generated files.
- Consult [doc/cmake.md](../../doc/cmake.md), the current CMake files, and
  [.github/workflows/push-pull.yml](../../.github/workflows/push-pull.yml) for build
  configuration. Reuse compatible dependency installations and build environments
  when available, or configure them for the current environment. Do not assume
  another developer's paths or an old copied command are current.
- Follow established naming in the same class, struct, or namespace. That local
  convention takes precedence over the defaults: static member functions and free
  functions use `snake_case()`, including private static helpers and function
  templates; non-static member functions use `camelCase()`. For example,
  `PETScMatUtilities` uses `camelCase()` for static operator-construction functions,
  so `constructSCInterpOpAxis` follows its convention. Preserve required external
  interface names and existing APIs; this rule does not authorize broad renaming.
  Member data use `d_`, and static member data use `s_`.
  Enum enumerators use uppercase names such as `SKEW_SYMMETRIC`; follow
  the existing `enum_to_string` and `string_to_enum` conventions when adding
  conversions. Follow the surrounding include style.
- Use `NDIM` for spatial dimension. Introduce a separate dimension template
  parameter only for an object whose dimension can actually differ from `NDIM`,
  not for hypothetical reuse.
- Use `auto` for lambdas, for range-based `for` loop variables, and when the exact
  type is evident from the right-hand side, such as
  `auto x = std::make_shared<T>(...);` or `auto x = T(...);`. Otherwise spell out
  the type: avoid opaque declarations such as `auto x = y;` or `auto x = y.z();`
  whose type the reader cannot see, and use the explicit index type for
  `box.lower()`. Deduction for genuinely dependent template types is also
  acceptable; do not add elaborate machinery or type erasure just to avoid `auto`.
- Use braces for all `if`, `else`, `for`, `while`, and `do` bodies, including
  single-statement and empty bodies. Ordinary `else if` chains are permitted.
- Prefer `using` aliases over `typedef`, and use `nullptr` and C++ named casts.
  Initialize variables when declared and use `const` for values that should not
  change. Choose initialization syntax for its meaning; braces are not mandatory
  for initialization, since container constructions such as `{n}` and `(n)` differ.
- Use `explicit` for constructors callable with one argument and conversion
  operators unless implicit conversion is an intentional part of the interface.
  Require `override` on overrides; use `final` only when preventing further
  inheritance or overriding is intentional.
- Prefer `enum class` for new enumerations, preserving established API and
  dependency requirements. Do not convert unrelated existing enums.
- Prefer the rule of zero: let members manage their resources and let the compiler
  supply special member functions when their behavior is appropriate. Use
  `= default` or `= delete` when explicitly specifying or disabling those operations;
  these declarations may appear with the interface in the declaration header.
- Prefer range-based loops when an index is unnecessary. Otherwise choose the
  index type appropriate to the container, numerical index, or called interface;
  no single signed or unsigned type fits every loop.
- Keep headers self-contained so they compile independently. Include the
  corresponding `ibamr/config.h` or `ibtk/config.h`. Use `namespaces.h` and
  `app_namespaces.h` only in source files, not library headers.
- Put class and struct definitions in declaration headers; an implementation-only
  nested type may remain in its owner's private section. Ordinary `.cpp`-local
  helper types stay in their owning implementation. Public headers declare
  interfaces, not function bodies. Put ordinary function definitions in `.cpp`
  files; put visible template, `constexpr`, and inline definitions in matching
  private inline headers included by the declaration header. A one-line body is
  still an implementation.
- Keep algorithm-specific helper types in their owning implementation, not the
  public API. A source-private shared header is appropriate for two genuine
  implementation consumers. Do not install it merely to make inclusion easier.
  Reusable kernel evaluation, for example, need not belong to PETSc matrix
  utilities.
- Pass only the information a helper needs; use an input filename rather than
  `argc` and `argv` when that is sufficient. A helper used only in one `.cpp` and
  needing no class access normally belongs in its unnamed namespace.
- Avoid `friend`. Prefer operations on the owning class or ordinary interfaces
  that belong to its contract. Do not add friendship or public accessors solely
  for tests. A required C++ customization or external interface may justify a
  narrow, documented exception.
- Order declarations `public`, `protected`, then `private`; keep definitions in
  declaration order. Avoid public member data. In declarations, omit top-level
  `const` on by-value parameters; definitions may use it for unchanged values.

Assign directly to members rather than adding aliases such as `auto& dx = d_dx;`
just to rename them. Spelling an explicit type for the same alias does not help:

```cpp
for (unsigned int d = 0; d < NDIM; ++d)
{
    d_dx[d] = dx0[d] / static_cast<double>(ratio(d));
}
```

## Name meaningful constants, not every literal

Use a constant when the name explains domain meaning or keeps related uses in
sync. Prefer `UPPER_CASE_SNAKE_CASE` at the narrowest useful scope, subject to
established local conventions. Do not replace constants with macros or expose
them publicly just for a test.

```cpp
constexpr int MAX_RETRIES = 3;
for (int attempt = 0; attempt < MAX_RETRIES; ++attempt)
{
    // ...
}
```

Keep ordinary `0`, `1`, and clear mathematical literals: `values.size() + 1`
for a terminator does not need `constexpr int ONE = 1`. When a policy needs
explanation, explain why it exists instead of restating its bound or current
numeric value. Stable mathematical facts, such as a three-dimensional formula,
may of course contain numbers. Do not rename every literal in a file while
fixing one function.

The same judgment applies beyond constants: use a small helper for a genuinely
repeated operation, not a wrapper that merely renames a single comparison.

## Fortran and m4 kernels

- Edit the Fortran/m4 sources, not generated files. Inspect generated fixed-form
  `.f` files after macro changes: keep portable fixed-form statements within
  72 columns and follow the compiler's continuation rules.
- Preserve CMake dependency tracking when sharing m4 code. The include scanner
  in `CMakeLists.txt` recognizes include lines using directory macros such as
  `TOP_SRCDIR`, `CURRENT_SRCDIR`, and `SAMRAI_FORTDIR`; follow that pattern.
