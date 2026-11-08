# Comments and class documentation

- Write comments for information a caller or maintainer needs. Use direct
  mathematical and library terminology; explain unfamiliar terms when needed.
  State concrete behavior, requirements, or reasons for non-obvious choices.
  Omit development history, hypothetical extensions, and prose that merely
  restates the declaration or promotes the design.
- Keep comments short, and prefer code that is clear without one: a clear name or
  a simpler structure is better than an explanation. Test every comment and every
  sentence of documentation:
  - If it is not needed to understand how to use the API, it does not belong in
    the header.
  - If a header must explain the implementation, the API is probably wrong: raise
    that instead of writing the comment.
  - If it is not needed to explain how the code in an implementation file works,
    leave it out.
  - If removing it would not reduce a reader's understanding of the interface or
    the implementation, remove it.
- Follow the "Shape" and "Language" sections below.
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
  the current change. Document comparable behavior in comparable detail, not in
  identical wording or length, and add detail only for genuine differences. If
  siblings need the same explanation, consider putting it in the base class.
  Rewriting every sibling is not required.
- Match the documentation of neighboring members in the same file: the Doxygen
  commands they use, a one-sentence brief followed by details where the file does
  that, and its line width. When shortening a comment, do not remove lines that
  all its neighbors contain, such as `\see` or `\param` lines.
- Document member functions. A concise Doxygen sentence is enough for most;
  documenting every member does not mean repeating class documentation. Private
  members may be documented in groups: every private member belongs to a
  documented group, and is documented separately only when the group's comment
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

For helpers that assemble matrix rows or own an evaluator, make the operation
explicit, for example:

| Before | After |
| --- | --- |
| Assemble one velocity component with compile-time direction choices. | Assemble interpolation-matrix rows for velocity component `Axis`. |
| Bind immutable evaluator storage to generic matrix construction. | Return a matrix builder that owns a const evaluator. |

Verify mathematical and behavioral claims against the implementation or tests.
Define specialized method terminology and cite the relevant method when needed.
Document necessary preconditions that the code cannot check. Apply the same
accuracy review to changelog entries and PR descriptions.

## Shape

- A class comment starts with a `\brief` sentence on what the class computes or
  represents. Add, only as needed, one or two sentences on its inputs, a short
  list of requirements, one note on accuracy or limits, and a small example. A
  class comment may be longer than six lines.
- A member comment is one sentence, plus what the caller must supply or may
  rely on. State a member's behavior at that member, not in the class comment.
- An implementation comment is one or two lines.
- Define each symbol in the comment that first uses it, or in the class comment.
- Proofs and reasoning that only a reviewer needs go in the commit message or PR
  description, not in a comment.
- These limits also apply to comments in test sources and input files.
  `scripts/maintenance/check_pull_request.py` reports added comment blocks
  longer than six lines, counting delimiter lines. A shorter comment can still
  be too long.
- Short labels that mark sections of a function or the kind of operation that
  follows, such as allocation, ghost filling, computation, or clean-up, help a
  reader. Give sibling classes the same layout and the same labels in the same
  words. Label a step in every sibling if any of them performs it in a way that
  warrants a comment; otherwise label it in none. Remove a label that only
  repeats one obvious call, in a function short enough to read at a glance, from
  every sibling at once.
- Choose wording for its content, not for how it wraps: do not reword a sentence
  only to change its line breaks.

## Language

These rules apply to comments, error messages, commit messages, changelog
entries, PR titles and descriptions, GitHub comments, and review replies. Apply
the Language rules to each new or changed sentence, and the Shape rules to each
new or changed comment. `scripts/maintenance/check_pull_request.py` reports
common violations.

Revising and rewriting differ. When revising text, such as comments a person
wrote or a review of someone else's PR, change only real errors and violations of
these rules, and do not replace technically correct wording only because another
phrasing reads better. When rewriting text that falls short of this guide, such
as agent-drafted comments, apply every rule and the Shape limits in full. In
either mode, base each change on a rule.

1. Only people and code act: reading and deciding belong to the person or code
   that does them. With data, values, cells, faces, boxes, patches, levels,
   stencils, or markers as the subject, use a verb that states a property, a
   relation, or an effect: a box contains a cell, a stencil extends to a face,
   the fine boxes determine the stencil, a flag indicates unset data, and a
   value comes from the input file. When the subject does not act, name what
   does: write "the FAC preconditioner requires no scratch data on a single
   level", not "a single level uses no scratch data".
2. Code, which includes functions, classes, tests, programs, commits, and PRs,
   may act with any verb except decide, answer, say, ask, know, want, and
   expect. For those, write determines, returns, checks, or calls.
3. Do not use hold, keep, retain, need, or live when the subject is code or one
   of the things listed in rule 1. Do not use kept or retained as adjectives for
   them either: write "the recorded levels" or "the remaining terms". Use "take"
   only for the arguments of a function, the values of a mathematical function,
   or a count of iterations. For hold, write contain or store. For keep, write
   is unchanged or still uses. For need, write is used or is required. With
   code or a configuration as the subject, write requires, as in "the function
   requires two ghost cells" or "the KSP type requires a linear
   preconditioner". For take a value, write uses. For live in, write is defined in or
   exists only during.
4. Write that results change, not that they move: "the expected output
   changes", and "the flux changes from 0.1 to 0.2".
5. Rules 1 through 4 do not apply to standard usage in mathematics, science, and
   engineering:
   - an identity or estimate holds, and two approximations agree to round-off;
   - a function takes the value 0 at x, takes two arguments, vanishes, or maps
     one space to another, and a matrix takes a form;
   - a method approximates or interpolates a quantity, and a scheme is
     second-order accurate;
   - a sequence converges, a solver converges in two iterations, a V-cycle
     visits each level, and a NaN propagates;
   - an error or a convergence rate rises or falls, and a norm sums or
     accumulates contributions;
   - a setting takes precedence or takes effect;
   - a point lies in a box, a level covers a cell, and a patch touches a
     boundary;
   - an object owns another, a function takes ownership of an object, and a
     function requires a precondition;
   - a marker moves, a test fails, and a program stops with an error.
6. Give each sentence one main clause and at most one relative clause. Clauses
   with if, when, or because do not count. A colon may introduce a second main
   clause that explains or illustrates the first. A semicolon may not. Split
   longer sentences. Call each thing by one name throughout.
7. Use definite, specific, and concrete language without needless words. Prefer
   the positive form when it says the same thing: write "is unchanged", not
   "does not change". A plain negative such as "no error was emitted" is fine.
   Use the active voice when the actor matters. Otherwise use the passive voice,
   as in "the first time that a face beyond the corner is filled".
8. Use American spelling and the serial (Oxford) comma. Do not use the first
   person ("I", "we", "our", or "us") in any text these rules cover. In a
   comment, make the code or the method the subject: "the solver computes", not
   "we compute".
9. Describe behavior in the present tense: "returns", not "will return".
   Function briefs, commit subjects, PR titles, the opening sentences of a PR
   description, and changelog entries use the imperative.
10. Omit "Note that", "simply", "just", "basically", and "of course".
11. Do not begin a sentence with a bare "This" or "It" that refers to the
    previous sentence. Write "This fill" or name the thing.
12. A comment, message, or document may state, say, or describe.
13. Use "where" for place and "when" or "while" for time: "the cells where the
    fine level covers the coarse level", "when the solver converges". For a
    condition, write "if".
14. Use "problem" for a recognized class of mathematical problems. For a
    specific mathematical object, write "system", "operator", or "equations";
    for a simple numerical experiment, "test"; for a fuller experiment tied to
    published results, "benchmark"; for a formulation, "model".
15. Punctuate a displayed formula as part of the sentence that contains it.

| Avoid | Write |
| --- | --- |
| The scratch data hold the fine values. | The scratch data contain the fine values. |
| the patch that holds the marker | the patch whose box contains the marker |
| The boxes of the fine level decide which cells the stencil reads. | The stencil uses only cells that the fine level covers. |
| A patch without viscosity data takes mu = 1. | The operator uses mu = 1 on a patch without viscosity data. |
| The boundary objects answer true. | `dependsOnTargetPatchData()` returns true. |
| the kept levels | the recorded levels |
| the box that the ghost cells need | the box that contains the ghost cells |
| The stencil lives in the new file. | The stencil is defined in the new file. |
| The two sets agree. | The two sets are equal. |
| A single level uses no scratch data. | The FAC preconditioner requires no scratch data on a single level. |
| Two expected outputs move. | Two expected outputs change. |
| The test expects the error. | The test checks that the error occurs. |

## Terminology

- Use "Eulerian grid" for the fluid grid and "Lagrangian mesh", "markers", or
  "nodes" for the structure. Do not call the structure "the boundary" when it may
  be volumetric.
- Spreading maps Lagrangian force densities to Eulerian force densities with a
  regularized delta function. Interpolation maps Eulerian velocities to
  Lagrangian velocities. The continuum coupling through the Dirac delta
  distribution is not spreading. State whether a quantity is a force, a force
  density, or a nodal force.
- Define "explicit", "semi-implicit", and "fully implicit" by the terms evaluated
  at the start and end of the time step. Distinguish the time step from the
  time-step size, and write "time-stepping scheme".
- Use "row" only for a row of a matrix. Cells parallel to a boundary form layers,
  as in "the two layers of ghost cells" and `num_inward_layers`.
- Write "geometric multigrid". Name a Schur complement by its variables, as in
  "pressure Schur complement". A matrix acts on degrees of freedom, not on a
  physical domain.
- In PR descriptions and changelog entries, write "indicates" for an empirical
  result and "proves" only for an analytical one. Limit claims of accuracy,
  stability, or scalability to the cases tested.

## Dimensions

In Doxygen comments and documentation pages, write the number of dimensions in
words rather than as 1d, 2d, 3d, 1D, 2D, or 3D. For the dimension of space,
write "one spatial dimension", "two spatial dimensions", or "three spatial
dimensions". For any other dimension, such as that of a cell, a surface, or a
quadrature rule, write "one-dimensional", "two-dimensional", or
"three-dimensional". Do not call a dimension spatial unless it is the dimension
of space. Identifiers, input keys, and test names such as `NDIM`, `SETUP_2D`, and
`_2d` are unchanged.

## Formula notation

Formulas in Doxygen comments distinguish quantities by typeface:

| Quantity | Typeface | Letters | Macro | Example |
| --- | --- | --- | --- | --- |
| Scalar | italic serif | any | none | `p`, `\mu` |
| Vector | bold italic serif | any | `\vect{}` | `\vect{u}`, `\vect{\chi}` |
| Second-order tensor | bold italic sans-serif | any | `\tens<letter>` | `\tensF`, `\tenssigma` |
| Fourth-order tensor | blackboard bold | uppercase Latin | `\ftens{}` | `\ftens{C}` |
| Matrix | bold upright sans-serif | uppercase Latin | `\mat{}` | `\mat{A}` |
| Discrete vector | bold upright sans-serif | lowercase Latin | `\dvec{}` | `\dvec{x}` |

1. The rules apply to every quantity, whether it is a function, a field, or a
   constant: a fixed normal `\vect{n}` is a vector and the identity `\tensI` is a
   tensor, while the identity matrix is `\mat{I}`.
2. Write formulas in Doxygen comments with the macros, so that the typefaces can
   be changed in one place. `doc/doxygen/scripts/make_formula_macros.py`
   generates their definitions in `doc/doxygen/formula-macros.tex`; change the
   script, not the generated file.
3. Markdown files that Doxygen does not read, such as `README.md` and the files
   in `doc/agents`, cannot use the macros. Write plain TeX there:
   `\boldsymbol{}` for vectors, `\boldsymbol{\mathsf{}}` for matrices and discrete
   vectors, and the Unicode character for a second-order tensor.
4. Components of vectors and tensors are scalars: write `F_{iJ}` and
   `\sigma_{ij}`.
5. Labels in subscripts and superscripts are upright text, written with
   `\mathrm{}`, as in `\mu_\mathrm{f}` and `\vect{u}_\mathrm{b}`. Indices and
   variables stay italic, as in `u_i` and `\tensF_t`.
6. An overbar or other accent goes on the base symbol only: write `\bar{I}_1`,
   not `\overline{I_1}` or `\bar{I_1}`.
7. The transpose is an upright sans-serif T: write `\tensF^\tran` and
   `\tensF^\invtran`, not `^T` or `^{-T}`.

Apply these rules to new and changed documentation; do not reword untouched
documentation only to follow them.
