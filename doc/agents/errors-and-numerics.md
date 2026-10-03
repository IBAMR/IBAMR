# Errors, ownership, and numerical behavior

- Prefer validating construction-time arguments and invariants during construction;
  omit repeated checks of unchanged validated state. Continue checking new
  arguments and state or lifecycle requirements that can change. Mutations must
  preserve invariants.
- Use `static_assert` at the earliest meaningful compile-time scope for compile-time
  requirements. Retain `static_assert` checks of function-template arguments only
  known at use; they have no runtime overhead.
- Use `TBOX_ERROR()` for fatal IBAMR runtime errors, with the surrounding
  stream-style formatting. Avoid C++ exceptions and catch-and-rethrow/translate
  scaffolding unless an existing external interface requires that boundary.
  Retain `TBOX_ERROR()` for invalid external input that must fail in Release.
- Use `noexcept` when required by an interface or when it has a concrete
  performance benefit, such as enabling standard-library containers to move
  elements instead of copying them. Verify that all operations performed support
  the nonthrowing guarantee. Do not annotate ordinary functions solely because
  IBAMR avoids exceptions, or add error-handling machinery just to make an
  operation `noexcept`.
- Call `TBOX_ERROR()` on the rank that detects the error. Do not broadcast an
  error string or add a collective just so all ranks report it. Communication
  needed to establish successful shared state or a global predicate is different:
  use the required reduction with consistent participation before evaluating a
  genuinely global condition, such as subdomain coverage.
- Handle PETSc failures through `IBTK_CHKERRQ()` or the established local macro.
  Make owned versus borrowed objects explicit, preserve lifetime requirements,
  and release owned resources through the normal teardown path.
- Pass opaque PETSc handles (`Mat`, `Vec`, `KSP`, `PC`, `SNES`, `IS`, `AO`, etc.)
  by value when the caller's handle is not replaced, even when modifying the
  object's entries or state. Modifying a `Vec`'s entries does not replace its
  handle. Borrowed getters normally return handles by value too: for example,
  `void setOperatorMat(Mat A);` and `Mat getOperatorMat() const;`.
  In ordinary new IBAMR C++ interfaces, use `Handle&` for a required output/in-out
  handle that is created, replaced, or reset, as in `void assembleOperator(Mat& A);`.
  Use `Handle*` when the output slot itself is optional or a required external
  interface dictates it; optional outputs are not a preferred API everywhere.
- Copying a handle does not copy its object. Passing by value or reference, or
  applying `const` to the pointer typedef, establishes neither object immutability
  nor ownership. Document borrowing, retention, transfer, and lifetime where the
  contract needs them; do not change reference counts to satisfy a signature rule.
  Handle arrays and containers retain their buffer/container interfaces, such as
  `const std::vector<Mat>&` for a read-only container. Use ordinary value-type rules
  for `PetscInt`, `PetscScalar`, and structs.
- Apply the handle rules consistently to new or reworked interfaces and matching
  declarations, definitions, and callers. Preserve required PETSc/C callback and
  override signatures and unchanged inherited APIs. Recognize necessary interface
  exceptions without using them to justify new inconsistencies; this is not a
  repository-wide ABI migration or a rewrite of native PETSc signatures.
- Prefer RAII and `std::unique_ptr` for new exclusive ownership. Use
  `std::shared_ptr` only for actual shared ownership. Borrowed pointers and
  references remain appropriate; preserve established SAMRAI and PETSc ownership
  conventions rather than imposing standard smart pointers on their interfaces.
- Preserve numerical contracts, signs, scaling, and boundary conditions unless
  intentionally correcting them. Explain and test an approved correction rather
  than preserving a known defect or silently changing an expected result.
- Register variables determining the maximum required ghost width before SAMRAI
  geometry computes that width; a later increase can be rejected.
- Composite-vector norms can mask coarse cells covered by finer levels through
  control-volume weights. For a whole-level comparison, select the intended
  levels and masking explicitly; an appropriate operation with volume index `-1`
  can include cells that a composite norm would omit.
- Choose exact, absolute-tolerance, or relative-tolerance floating-point comparisons
  according to the mathematical requirement. Reuse suitable existing helpers;
  do not prescribe one comparison helper or tolerance for every quantity.
- Use an extension's declared contract, such as an operator's stencil width,
  rather than inferring capabilities from its name. Correct a misreported contract
  at its source. Validate supported-mode restrictions at the earliest valid
  lifecycle point.
- Do not select a code path by testing an object's dynamic type with `typeid` or
  `dynamic_cast`, for example to take a faster path only for the base class or to
  keep work that a derived class might use. Remove work that nothing uses; when
  implementations differ in what they support, express that in the interface.
  Before adding a virtual function or a capability query so that callers can skip
  work, check whether every existing implementation already supports the cheaper
  path; if so, use it directly.
- Review frequently called code affected by the change for repeated conversions,
  lookups, allocations, and indirect calls. Move invariant work outside loops
  when practical, and use profiling when an uncertain cost could matter to the
  change. This does not require a benchmark for every PR or prohibit indirect
  calls where they are appropriate. Include repeated setup and rebuild costs;
  choose data structures for their access pattern and allocation cost, especially
  for per-degree-of-freedom storage. Debug timings are not performance evidence.
  Weigh the work an optimization saves, such as the number of ghost fills or
  copies per time step it removes, against the interface and code it adds.

```cpp
// Prefer a direct fatal error at the point of detection.
if (width <= 0)
{
    TBOX_ERROR("Interpolator::initialize():\n" << "  kernel width must be positive\n");
}
```

Avoid throwing `std::runtime_error`, catching it in the caller, broadcasting
`what()`, and finally calling `TBOX_ERROR()` for the same invalid width.

## Configuration and small utility functions

- Check existing input keys and their meanings before adding a new mechanism.
  Keep related settings with the component they configure and preserve established
  precedence, including command-line PETSc overrides.
- Preserve and document each component's precedence between restart data and
  input settings. An intentional change needs justification and tests; do not
  impose a universal restart/input ordering on existing components.
- Validate the rules IBAMR introduces. Let SAMRAI handle database syntax and PETSc
  handle PETSc option syntax and values instead of implementing a second parser
  or a more restrictive set of rules in IBAMR.
- Prefer direct processing when no intermediate representation is needed. Use a
  standard-library conversion or algorithm when it fits, rather than adding a
  custom text helper or Boost usage. Check that conversion preserves the intended
  value; do not add normalization or locale handling without a concrete need.
- Prefer `std::filesystem::path` consistently for filesystem operations, converting
  to strings where an interface requires them. Use nonthrowing overloads when
  handling filesystem errors through `TBOX_ERROR()`, and check the reported error.

| Prefer | Avoid |
| --- | --- |
| Walk a database directly when finding and applying settings is the whole operation. | Collect, sort, normalize, and replay settings without an ordering requirement. |
| Add a distinct input key for a new mechanism. | Reinterpret an existing key such as `petsc_options_prefix`. |
| Use `std::to_chars` with appropriate precision to format a floating-point setting. | Use `std::to_string` for a tiny tolerance without checking whether it becomes `"0.000000"`; write a custom numeric formatter unnecessarily. |
| Use an enum for a genuinely closed set of choices. | Force an extensible application-supplied catalog into a closed enum. |
