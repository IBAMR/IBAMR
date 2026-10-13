# Comments and class documentation

- Write comments for information a caller or maintainer needs. Use direct
  mathematical and library terminology; explain unfamiliar terms when needed.
  State concrete behavior, requirements, or reasons for non-obvious choices.
  Omit development history, hypothetical extensions, and prose that merely
  restates the declaration or advertises the design.
- Keep comments short, and prefer code that needs none: a clear name or a simpler
  structure is better than an explanation. Test every comment and every sentence
  of documentation:
  - If it is not needed to understand how to use the API, it does not belong in
    the header.
  - If a header needs it to explain the implementation, the API is probably wrong:
    raise that instead of writing the comment.
  - If it is not needed to explain how the code in an implementation file works,
    leave it out.
  - If removing it would not reduce a reader's understanding of the interface or
    the implementation, remove it.
- A Doxygen comment is normally one sentence, plus what the caller must supply or
  may rely on. An implementation comment is normally one or two lines. Put
  reasoning that only a reviewer needs in the commit message or PR description.
  The same applies to comments in test sources and input files.
- Do not carry working vocabulary from the development of a change into
  comments, names, input keys, file names, or messages: for example `parity`
  for agreement between two code paths, `live`, `pin`, `probe`, or the labels
  of a PR stack.
- A class header should explain what the class represents and how to use or
  extend its interface. Include mathematical meaning, inputs, results, limitations,
  ownership, lifetime, invalidation, and calling requirements only where a user or
  implementer of the interface needs them. Describe allocation, caching, storage,
  or forwarding only when it establishes a requirement or guarantee the reader
  needs. Include a small, complete example when it makes usage clearer.
- Review API documentation by asking, in order: "Would removing this sentence leave
  someone unsure how to use or extend the class correctly?" If so, "Is that
  information already explained where the API is defined, usually in a base class?"
  Review an implementation comment the same way: would removing it leave a
  maintainer unsure how the code works?
- Document shared API requirements where the API is defined. Derived classes
  explain particular behavior, additional requirements, and meaningful differences
  without repeating the base class. Explain requirements spanning several methods
  once at class or group level. Keep shared configuration requirements in the
  relevant configuration documentation rather than repeating them in constructors.
- Keep each class's documentation about that class. A base class does not
  describe what particular derived classes or callers do, a derived class does
  not re-explain its base, and a class does not describe the internals of a
  sibling or collaborator. State a base contract in terms that every
  implementation satisfies. If a base class cannot be documented without
  describing a derived class, the interface is probably wrong: raise that instead
  of working around it in a comment.
- If shared documentation is missing or unclear, consider improving the base class.
  This normally belongs in a separate focused PR unless the clarification is needed
  to document the current change correctly. Do not silently broaden the change or
  defer a required clarification to later cleanup.
- Compare documentation of sibling classes implementing the same API while reviewing
  the current change. Comparable behavior merits comparable detail, not identical
  wording or length; extra detail should explain genuine differences. If siblings
  need the same explanation, consider putting it in the base class. This does not
  require rewriting every sibling.
- Match the documentation of neighboring members in the same file: the Doxygen
  commands they use, a one-sentence brief followed by details where the file does
  that, and its line width. When shortening a comment, keep what its neighbors
  consistently carry, such as their `\see` or `\param` lines.
- Document member functions. A concise Doxygen sentence is enough for most;
  documenting every member does not mean repeating class documentation. Private
  members may be documented in groups: every private member belongs to a
  documented group, and needs a comment of its own only when the group's comment
  does not cover it. Use Doxygen links, `\see`, or `\copydoc` when useful, with
  valid targets; no automatic cross-reference is needed at every override. State
  additional requirements of an overload beside that overload.
- Explain implementation choices and non-obvious mathematical relationships
  briefly, beside the implementation; keep private-data descriptions brief and
  local. Omit obvious mechanics, commented-out code, and obsolete explanations.

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
