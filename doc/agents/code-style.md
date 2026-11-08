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
- New code should look like similar code in the same header or implementation
  file or, for a new file, in similar files: naming, structure, member order, the
  Doxygen commands used (`\brief`, `\param`, `\return`, `\see`), and the level of
  detail of comments. A reader should not be able to tell from style which parts
  of a file are new. This does not justify a comment longer than
  [Comments and class documentation](documentation.md) allows, and it does not
  authorize restyling existing code.
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
- End a line of streamed output with `"\n"`, not `std::endl`, which also flushes
  the stream. The flush slows output in loops and does nothing in `TBOX_ERROR()`
  and `TBOX_WARNING()`, which build their message in a string stream. Write
  `std::flush` only where output must appear at that point (C++ Core Guidelines
  SL.io.50).

Assign directly to members rather than adding aliases such as `auto& dx = d_dx;`
just to rename them. Spelling an explicit type for the same alias does not help:

```cpp
for (unsigned int d = 0; d < NDIM; ++d)
{
    d_dx[d] = dx0[d] / static_cast<double>(ratio(d));
}
```

## Choose names

A name states what the thing is or does in the words of its class and siblings.
Names follow the [Language](documentation.md#language) rules on word choice and
the [Terminology](documentation.md#terminology) rules. Like comments, they carry
no working vocabulary from the development of a change. These rules apply to new
and changed names.

- Name types and data other than flags and enumerators with a noun phrase, and
  functions with a verb. An adjective or participle alone is not enough: write
  `shifted_geometry`, not `shifted`, and `fine_data`, not `fine`.
- A short, precise name is better than a long, vague one. Short index names are
  fine. So are short fields of a small struct if the struct's name states their
  context, as in `NormalVelocityStencil::idx` and `weight`.
- Use an abbreviation only if it is established, or if it is more precise than
  a longer alternative and is explained where it is defined: `tt_corners`
  (traction-traction corners) is more precise than `traction_corners`.
- Name an enumerator for its case, not for its condition: `INTERIOR`, not
  `RATIO_4_OR_MORE`.
- Name a flag for its use or for the property that it indicates:
  `use_source_ghost_layer`, `restrict_to_active_patches`, or
  `has_lagrangian_data`.
- A name must not state an action that the code does not perform: write
  `evaluateBcCoefficients()`, not `readBcCoefficients()`. Do not use a verb that
  Language rules 2 and 3 exclude for code, such as decide, hold, keep, or need.
- If a class or function could refer to more than one object of a kind, say
  which one: `d_shifted_geometry` and `d_original_geometry`, or
  `getShiftedGeometry()`.
- Name a copy for its result, a builder for its product, and a container for its
  contents: `d_shifted_patch`, `buildShiftedPatch()`, and `d_flagged_boxes`.
- Build paired and sibling names the same way. Give the same kind of thing the
  same name and a different kind a different name: `componentUsesCubicScratch()`
  and `levelUsesCubicScratch()`, not `needsCubicScratch()` and
  `usesCubicScratch()`.
- Use established IBAMR and SAMRAI names and the notation of the method:
  `periodic_shift` (not `period`) for a periodic shift, `*_idx` for patch data
  indices, `ln` for the level number in a loop and `finest_ln` for the number of
  the finest level, `p_` for a pointer obtained by a cast, `a`, `b`, and `g` for
  the Robin coefficients, and "limit" for a limiting operation.
- A generic word such as `Scope` or `Info` is not enough as a type name, alone or
  as a suffix: write `ScopedGeometryReplacement`, not
  `ShiftedPatchGeometry::Scope`, and `LevelSchedules`, not `LevelInfo`. Name a
  struct for one kind of data, such as schedules, only if it contains nothing
  but that data and closely related data, as `LevelSchedules` does.
- Choose plain words that suggest no unintended image or technique: write
  `DepartingMarker`, not `Leaver`. Do not borrow the vocabulary of another
  technique: `Tape`, `slot`, and `replay` together suggest automatic
  differentiation.

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
