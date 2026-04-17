#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Plot sparsity pattern of a Matrix Market sparse matrix (e.g. calibration Jacobian).

Dependencies: scipy, matplotlib

Examples:
  python scripts/plot_calibration_jacobian_spy.py /tmp/calib_J.mtx
  python scripts/plot_calibration_jacobian_spy.py /tmp/calib_J.mtx --mathtext -o /tmp/spy.png

Very large matrices may be slow or memory-heavy for matplotlib's spy().

By default, text is rendered with your system LaTeX (matplotlib text.usetex;
requires pdflatex and related tools on PATH). Pass --mathtext to use only
matplotlib's built-in mathtext (no TeX install).

Fonts: Computer Modern (serif) for math; serif text (LaTeX defaults when using
TeX; mathtext ``cm`` + matplotlib serif family when using --mathtext).
"""

from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Spy plot of a Matrix Market (.mtx) sparse matrix."
    )
    parser.add_argument(
        "input_mtx",
        help="Path to Matrix Market file (e.g. from export_calibration_jacobian).",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output PNG path (default: same basename as INPUT with .png).",
    )
    parser.add_argument(
        "--mathtext",
        action="store_true",
        help="Use matplotlib mathtext only (no system LaTeX). Default is LaTeX.",
    )
    args = parser.parse_args()

    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy import sparse as sp_sparse
        from scipy.io import mmread
    except ImportError as e:
        print("Requires scipy and matplotlib:", e, file=sys.stderr)
        return 1

    use_latex = not args.mathtext

    # Serif + Computer Modern (not sans-serif)
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["axes.unicode_minus"] = False
    if use_latex:
        mpl.rcParams["text.usetex"] = True
        mpl.rcParams["text.latex.preamble"] = (
            r"\usepackage{amssymb}"
            r"\usepackage{amsmath}"
            r"\renewcommand{\familydefault}{\rmdefault}"
        )
    else:
        mpl.rcParams["mathtext.fontset"] = "cm"

    if not os.path.isfile(args.input_mtx):
        print(f"Not a file: {args.input_mtx}", file=sys.stderr)
        return 1

    a = mmread(args.input_mtx)
    out = args.output
    if out is None:
        base, _ = os.path.splitext(args.input_mtx)
        out = base + ".pdf"

    n_rows, n_cols = a.shape
    if sp_sparse.issparse(a):
        nnz = a.nnz
    else:
        nnz = int(np.count_nonzero(a))

    basename = os.path.basename(args.input_mtx)
    basename_tex = basename.replace("_", r"\_") if use_latex else basename

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.spy(a, markersize=0.5, aspect="auto")
    title_line1 = r"$\mathbf{J}\ \mathrm{sparsity}$"
    title_line2 = (
        rf"$\mathbf{{J}} \in \mathbb{{R}}^{{{n_rows}\times{n_cols}}}$, "
        rf"$|\mathrm{{nnz}}| = {nnz}$"
    )
    if use_latex:
        title_line3 = rf"\textrm{{{basename_tex}}}"
        title = title_line1 + "\n" + title_line2 + "\n" + title_line3
    else:
        title = title_line1 + "\n" + title_line2 + "\n" + basename_tex
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(r"parameter index $j$ (column)", fontsize=12)
    ax.set_ylabel(r"residual index $i$ (row)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
