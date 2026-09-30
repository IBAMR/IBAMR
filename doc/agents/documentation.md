# Comments and class documentation

- Write comments for information a caller or maintainer needs. Use direct
  mathematical and library terminology; explain unfamiliar terms when needed.
  State concrete behavior, requirements, or reasons for non-obvious choices.
  Omit development history, hypothetical extensions, and prose that merely
  restates the declaration or advertises the design.
- A class header should explain what the class represents and how to use or
  extend its interface. Include mathematical meaning, inputs, results, limitations,
  ownership, lifetime, invalidation, and calling requirements as needed. Describe
  allocation, caching, storage, or forwarding only when it establishes a requirement
  or guarantee the reader needs. Include a small, complete example when it makes
  usage clearer.
- Review API documentation by asking, in order: "Would removing this sentence leave
  someone unsure how to use or extend the class correctly?" If so, "Is that
  information already explained where the API is defined, usually in a base class?"
  These questions concern API documentation, not helpful implementation comments.
- Document shared API requirements where the API is defined. Derived classes
  explain particular behavior, additional requirements, and meaningful differences
  without repeating the base class. Explain requirements spanning several methods
  once at class or group level. Keep shared configuration requirements in the
  relevant configuration documentation rather than repeating them in constructors.
- If shared documentation is missing or unclear, consider improving the base class.
  This normally belongs in a separate focused PR unless the clarification is needed
  to document the current change correctly. Do not silently broaden the change or
  defer a required clarification to later cleanup.
- Compare documentation of sibling classes implementing the same API while reviewing
  the current change. Comparable behavior merits comparable detail, not identical
  wording or length; extra detail should explain genuine differences. If siblings
  need the same explanation, consider putting it in the base class. This does not
  require rewriting every sibling.
- Document member functions, including private ones. A concise Doxygen sentence
  is enough for a simple helper; documenting every member does not mean repeating
  class documentation. Use Doxygen links, `\see`, or `\copydoc` when useful, with
  valid targets; no automatic cross-reference is needed at every override. State
  additional requirements of an overload beside that overload.
- Explain implementation choices and non-obvious mathematical relationships beside
  the implementation; keep private-data descriptions brief and local. Omit obvious
  mechanics, commented-out code, and obsolete explanations.

For example, declare a simple query in the public header and define it in the
corresponding implementation file:

```cpp
/*! \brief Return the number of levels in the hierarchy. */
int getNumberOfLevels() const;
```

The short description above is sufficient. A comment explaining why a borrowed
matrix must outlive a solver is useful; `// Set the matrix` above
`setMatrix(matrix)` is not.

For helpers that assemble matrix rows or retain an evaluator, make the operation
explicit, for example:

| Before | After |
| --- | --- |
| Assemble one velocity component with compile-time direction choices. | Assemble interpolation-matrix rows for velocity component `Axis`. |
| Bind immutable evaluator storage to generic matrix construction. | Return a matrix builder that owns a const evaluator. |

Verify mathematical and behavioral claims against the implementation or tests.
Define specialized method terminology and cite the relevant method when needed.
Document necessary preconditions that the code cannot check. Apply the same
accuracy review to changelog entries and PR descriptions.
