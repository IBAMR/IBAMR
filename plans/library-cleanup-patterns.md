# Library Cleanup Patterns

This note collects small implementation patterns that came up repeatedly while
cleaning solver and preconditioner code. The goal is to make these choices
deliberate and reusable across later cleanup passes.

## Prefer explicit ownership

- Prefer `std::unique_ptr` over bare owning pointers.
- Use raw pointers only for non-owning references or when an external API
  requires them.
- When an API needs a raw pointer, prefer calling `.get()` on a smart pointer
  instead of storing the object as a raw pointer.

## Prefer direct PETSc destroy calls

- For PETSc objects initialized to `nullptr`, prefer direct calls to
  `XXXDestroy(&obj)` without redundant `nullptr` guards.
- After destroy, reset cached handles explicitly when the surrounding code
  benefits from showing that the object is no longer owned.

## Match teardown order to declaration order

- In `deallocate*()` methods, prefer destroying members in the same order they
  appear in the header.
- Use comments that describe the actual resources being released, not comments
  that describe bookkeeping strategy.

## Prefer descriptive names for cached data

- Accessors that return cached state should say so in the name, e.g.
  `getCachedPressureDOFs()`.
- Reserve shorter `get*()` names for lightweight accessors whose values are
  always already available or trivial to compute.

## Prefer simple shared helpers over ad hoc lambdas

- When two code paths differ only by type, prefer a small function template to
  duplicated code or a stateful lambda.
- Helper arguments should expose meaningful policy choices directly, e.g.
  pass ghost width explicitly instead of hardwiring it inside the helper.

## Prefer `clear()` to `resize(0)`

- For vectors and similar containers, use `clear()` when the intent is to
  remove all entries.
- Reserve `resize()` for cases where the target size itself matters.

## Keep comments concrete

- Prefer comments that name the actual resources being managed.
- Avoid comments that restate obvious ownership rules when the code already
  makes them clear.
