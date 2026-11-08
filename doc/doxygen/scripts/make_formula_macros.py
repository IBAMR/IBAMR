#!/usr/bin/env python3
# ---------------------------------------------------------------------
#
# Copyright (c) 2026 by the IBAMR developers
# All rights reserved.
#
# This file is part of IBAMR.
#
# IBAMR is free software and is distributed under the 3-clause BSD
# license. The full text of the license can be found in the file
# COPYRIGHT at the top level directory of IBAMR.
#
# ---------------------------------------------------------------------

"""Write the formula macros used in IBAMR's documentation.

Usage: python3 doc/doxygen/scripts/make_formula_macros.py > doc/doxygen/formula-macros.tex

Do not edit formula-macros.tex by hand: change this script and rerun it.

Doxygen reads the output through FORMULA_MACROFILE. See the "Formula notation"
section of doc/agents/documentation.md for the convention.

Second-order tensors use bold italic sans-serif letters. No TeX command gives
those for Greek letters in MathJax, so each \\tens<letter> macro names the
Unicode character directly. The other macros take an argument.
"""

import unicodedata

GREEK_UPPER = """Alpha Beta Gamma Delta Epsilon Zeta Eta Theta Iota Kappa Lambda
Mu Nu Xi Omicron Pi Rho - Sigma Tau Upsilon Phi Chi Psi Omega""".split()
GREEK_LOWER = """alpha beta gamma delta epsilon zeta eta theta iota kappa lambda
mu nu xi omicron pi rho varsigma sigma tau upsilon phi chi psi omega""".split()
# Greek letters that look like Latin ones have no TeX command of their own.
TEX_GREEK_UPPER = {"Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Upsilon", "Phi", "Psi", "Omega"}
# The variant forms follow the lowercase letters in the Unicode block. TeX's
# \epsilon and \phi are the lunate and closed forms, and \varepsilon and
# \varphi are the curly and open forms, the reverse of the Unicode names.
GREEK_VARIANTS = {"epsilon": 1, "vartheta": 2, "varkappa": 3, "phi": 4, "varrho": 5, "varpi": 6}
UNICODE_TO_TEX_GREEK_LOWER = {"epsilon": "varepsilon", "phi": "varphi"}

LATIN_UPPER_START = 0x1D63C  # MATHEMATICAL SANS-SERIF BOLD ITALIC CAPITAL A
LATIN_LOWER_START = 0x1D656  # ... SMALL A
GREEK_UPPER_START = 0x1D790  # ... CAPITAL ALPHA
GREEK_LOWER_START = 0x1D7AA  # ... SMALL ALPHA
GREEK_VARIANT_START = 0x1D7C3  # ... PARTIAL DIFFERENTIAL


def tensor_macros():
    macros = []
    for i in range(26):
        macros.append((chr(ord("A") + i), LATIN_UPPER_START + i))
    for i in range(26):
        macros.append((chr(ord("a") + i), LATIN_LOWER_START + i))
    for i, name in enumerate(GREEK_UPPER):
        if name in TEX_GREEK_UPPER:
            macros.append((name, GREEK_UPPER_START + i))
    for i, name in enumerate(GREEK_LOWER):
        if name not in ("omicron", "varsigma"):
            macros.append((UNICODE_TO_TEX_GREEK_LOWER.get(name, name), GREEK_LOWER_START + i))
    for name, offset in GREEK_VARIANTS.items():
        macros.append((name, GREEK_VARIANT_START + offset))
    for _, code_point in macros:
        assert "SANS-SERIF BOLD ITALIC" in unicodedata.name(chr(code_point))
    return macros


def main():
    # Doxygen accepts only \newcommand and \renewcommand lines in this file:
    # no comments or blank lines.
    #
    # Vector quantities: bold italic serif, any letter.
    print(r"\newcommand{\vect}[1]{\boldsymbol{#1}}")
    # Matrices (uppercase) and discrete vectors (lowercase): bold upright
    # sans-serif, Latin letters only.
    print(r"\newcommand{\mat}[1]{\boldsymbol{\mathsf{#1}}}")
    print(r"\newcommand{\dvec}[1]{\boldsymbol{\mathsf{#1}}}")
    # Fourth-order tensors: blackboard bold, uppercase Latin letters only.
    print(r"\newcommand{\ftens}[1]{\mathbb{#1}}")
    # Transpose and inverse transpose exponents, an upright sans-serif T, used
    # as \tensF^\tran and \tensF^\invtran. Doxygen wraps each definition in
    # braces, so a definition that includes the ^ would put the exponent on an
    # empty group rather than on the symbol.
    print(r"\newcommand{\tran}{\mkern-1.5mu\mathsf{T}}")
    print(r"\newcommand{\invtran}{-\mathsf{T}}")
    # Second-order tensors: bold italic sans-serif, one macro per letter,
    # e.g., \tensF and \tenssigma.
    for name, code_point in tensor_macros():
        print(rf"\newcommand{{\tens{name}}}{{\unicode{{x{code_point:X}}}}}")


if __name__ == "__main__":
    main()
