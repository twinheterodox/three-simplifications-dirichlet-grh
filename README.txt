Three simplifications in the Alpöge–Claude argument for Dirichlet GRH
twinheterodox
21 November 2026

Read output/pdf/three-simplifications-dirichlet-grh.pdf. The editable source
is paper/three-simplifications-dirichlet-grh.tex.

The paper is for graduate readers familiar with Dirichlet L-functions and
finite character theory. It follows the centered prime sum of a fixed family
of characters from a lattice inequality to a zero-free conclusion. Theorem
3.1 compares all thick contributions at one argument. Corollary 4.1 shows
how much the weights can vary while preserving the prime-sum exponent.
Proposition 2.1 gives the sharp bound for the integrated rank loss.
Proposition 5.1 explains why an upper bound on arbitrarily late power
intervals suffices, even for a family with a possible real zero.

Theorem 1.1 assumes the arithmetic data and estimates that produce those
intervals. The manuscript draws on the Alpöge–Claude release, Jaswal's
exposition, and the primary sources cited at the steps where they enter.
The companion paper “From theta lifts to a prime-counting bound” gives more
detail on the geometric input.

Build from this directory with pdfLaTeX and the packages listed in the
TeX preamble. A TeX Live installation including Latin Modern suffices.

  mkdir -p build output/pdf
  pdflatex -halt-on-error -interaction=nonstopmode -output-directory=build paper/three-simplifications-dirichlet-grh.tex
  pdflatex -halt-on-error -interaction=nonstopmode -output-directory=build paper/three-simplifications-dirichlet-grh.tex
  cp build/three-simplifications-dirichlet-grh.pdf output/pdf/

The finite calculations use only the Python 3 standard library. Run from
this directory.

  python3 checks/finite_average.py
  python3 checks/check_projection.py

The first program tests the Frobenius-coset identity, the Fourier
projection and the direction of the ramp comparison. The second independently
checks the conjugation convention and normalization with unequal weights
and nonreal characters. The general inequalities are proved in the text.
The scope of the mathematical and exposition reviews is in review/reading.txt.
