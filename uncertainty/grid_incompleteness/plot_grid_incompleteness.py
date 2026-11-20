# Copyright 2026 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Plots the grid incompleteness errors from grid_incompleteness.py.

The incompleteness error of a grid is its result less the result of the
reference grid of the same molecule: the largest radial size, with the largest
angular size run with it (e.g., 75 x 5810). The reference grid itself is not
plotted, since its error is zero by construction. For E_xc
(--quantity exc) that is the signed difference (or its absolute value, with
--absolute). For V_xc (--quantity vxc) it is the largest absolute difference
among the elements of the AO matrix. There is one panel per molecule, with the
error as a function of the radial size and one line per angular size.
"""

import argparse
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np

MY_DIR = os.path.dirname(os.path.realpath(__file__))

#: Where the plots are written by default
PLOT_DIR = os.path.join(MY_DIR, "plots")

MOLECULE_NAMES = {"he": "He", "h2": "H$_2$", "water": "H$_2$O"}

#: Categorical palette, assigned in this fixed order (coarse to fine angular
#: sizes, over all molecules, so an angular size has the same color in every
#: panel)
COLORS = [
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#6250d6",
    "#e34948",
]

#: Distinct markers, so lines that coincide (e.g., He) remain identifiable
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

#: Angular sizes plotted by default: roughly doubling, so at most one per
#: color. Sizes that were not run are skipped.
DEFAULT_PLOTTED_ANGULAR = [14, 26, 50, 110, 302, 590, 974, 2030]

#: Errors below this are at the SCF convergence/rounding floor, so the y axis
#: is linear (rather than logarithmic) between +/- this value
LINEAR_THRESHOLD = 1.0e-10

#: With --absolute, errors below this are drawn at this value, so that they
#: remain on the log axis
ABSOLUTE_FLOOR = 1.0e-12


def exc_error(result, reference):
    """The signed error in E_xc."""
    return result["E_xc"] - reference["E_xc"]


def vxc_error(result, reference):
    """The largest absolute error among the elements of V_xc."""
    return float(np.max(np.abs(result["V_xc"] - reference["V_xc"])))


#: For each quantity: the error function, whether the error is signed, and the
#: (LaTeX) symbol of the error used in the y-axis label
QUANTITIES = {
    "exc": (
        exc_error,
        True,
        r"E_\mathrm{XC} - E_\mathrm{XC}^\mathrm{ref}",
    ),
    "vxc": (
        vxc_error,
        False,
        r"\max_{\mu\nu}|V^\mathrm{XC}_{\mu\nu} - "
        r"V^\mathrm{XC,ref}_{\mu\nu}|",
    ),
}


def basis_label(basis):
    """The conventional spelling of the ChemCache basis key @p basis."""
    if basis.startswith("aug-cc-pv") or basis.startswith("cc-pv"):
        prefix, zeta = basis.rsplit("pv", 1)
        return prefix + "pV" + zeta.upper()
    return basis.upper()


def plot(
    results, output, plotted, absolute=False, quantity="exc", basis="sto-3g"
):
    error_fxn, signed, symbol = QUANTITIES[quantity]
    absolute = absolute or not signed
    molecules = [m for m in MOLECULE_NAMES if any(k[0] == m for k in results)]
    fig, axes = plt.subplots(
        1, len(molecules), figsize=(4 * len(molecules), 4.2), sharey=True
    )
    all_angular = sorted({k[2] for k in results} & set(plotted))
    if len(all_angular) > len(COLORS):
        raise ValueError("More angular sizes than colors; use --angular")

    legend = {}
    for ax, molecule in zip(axes, molecules):
        keys = [k for k in results if k[0] == molecule]
        radial = sorted({k[1] for k in keys})

        # The reference is the largest radial size with the largest angular
        # size that has been run with it, plotted or not
        ref_key = (
            molecule,
            radial[-1],
            max(k[2] for k in keys if k[1] == radial[-1]),
        )
        reference = results[ref_key]

        lines = []  # One {radial size: error} per line
        for n_angular in sorted({k[2] for k in keys} & set(plotted)):
            # Skips grids that have not been run, and the reference, whose
            # error is zero by construction
            r_ok = [
                r
                for r in radial
                if (molecule, r, n_angular) in keys
                and (molecule, r, n_angular) != ref_key
            ]
            errors = [
                error_fxn(results[(molecule, r, n_angular)], reference)
                for r in r_ok
            ]
            lines.append(dict(zip(r_ok, errors)))
            if absolute:
                errors = [max(abs(e), ABSOLUTE_FLOOR) for e in errors]
            i = all_angular.index(n_angular)
            (legend[n_angular],) = ax.plot(
                r_ok,
                errors,
                color=COLORS[i],
                marker=MARKERS[i],
                markersize=6,
                linewidth=2,
                label=str(n_angular),
                zorder=3,
            )

        # Lines that coincide hide each other, so say so. Only the radial
        # sizes that every line has are compared.
        shared = set.intersection(*(set(line) for line in lines))
        spread = max(
            (
                max(line[r] for line in lines) - min(line[r] for line in lines)
                for r in shared
            ),
            default=float("inf"),
        )
        if len(lines) > 1 and spread < LINEAR_THRESHOLD:
            ax.text(
                0.97,
                0.95,
                "All angular sizes coincide",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9,
                color="#52514e",
            )

        if absolute:
            ax.axhspan(0.0, LINEAR_THRESHOLD, color="#ecebe8", zorder=1)
            ax.set_yscale("log")
            ax.set_ylim(ABSOLUTE_FLOOR / 3.0, 0.3)
        else:
            ax.axhspan(
                -LINEAR_THRESHOLD, LINEAR_THRESHOLD, color="#ecebe8", zorder=1
            )
            ax.axhline(0.0, color="#8a8984", linewidth=0.8, zorder=2)
            ax.set_yscale("symlog", linthresh=LINEAR_THRESHOLD, linscale=0.5)
        # Log scale, so small and large radial sizes are both readable
        ax.set_xscale("log")
        ax.set_xticks(radial)
        ax.set_xticklabels([str(r) for r in radial])
        ax.minorticks_off()
        ax.set_title(
            "%s (reference %d x %d)"
            % (MOLECULE_NAMES[molecule], ref_key[1], ref_key[2])
        )
        ax.set_xlabel("Radial points")
        ax.grid(True, which="major", color="#e1e0dc", linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    if absolute:
        label = symbol if not signed else "|" + symbol + "|"
        axes[0].set_ylabel(r"$" + label + r"$ (E$_\mathrm{h}$)")
        axes[0].set_yticks([1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
        note = (
            "Shaded band: |error| < 1e-10 Eh (SCF convergence floor). "
            "Errors < 1e-12 are drawn at 1e-12."
        )
    else:
        axes[0].set_ylabel(r"$" + symbol + r"$ (E$_\mathrm{h}$)")
        ticks = [-1e-1, -1e-4, -1e-7, 0.0, 1e-7, 1e-4, 1e-1]
        axes[0].set_yticks(ticks)
        axes[0].set_yticklabels(
            ["$-10^{-1}$", "$-10^{-4}$", "$-10^{-7}$", "0"]
            + ["$10^{-7}$", "$10^{-4}$", "$10^{-1}$"]
        )
        note = (
            "Shaded band: |error| < 1e-10 Eh (linear scale; SCF convergence "
            "floor)."
        )
    note += " Reference grid (radial x angular) is not plotted."

    labels = sorted(legend)
    fig.legend(
        [legend[a] for a in labels],
        [str(a) for a in labels],
        title=basis_label(basis) + ": angular points",
        loc="upper center",
        ncol=len(labels),
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
    )
    fig.text(
        0.5,
        0.01,
        note,
        ha="center",
        fontsize=8,
        color="#52514e",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.9))
    fig.savefig(output, dpi=200)
    print("Wrote", output)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--basis",
        default="sto-3g",
        help="Basis set (ChemCache key) whose results are plotted",
    )
    parser.add_argument(
        "--input", help="Defaults to grid_incompleteness_<basis>.pkl"
    )
    parser.add_argument(
        "--output", help="Defaults to the plots/ directory next to this script"
    )
    parser.add_argument(
        "--angular",
        nargs="+",
        type=int,
        default=DEFAULT_PLOTTED_ANGULAR,
        help="Angular sizes to plot (at most %d)" % len(COLORS),
    )
    parser.add_argument(
        "--absolute",
        action="store_true",
        help="Plot the absolute error on a log scale (always, for vxc)",
    )
    parser.add_argument("--quantity", choices=QUANTITIES, default="exc")
    args = parser.parse_args()
    if args.output is None:
        name = "grid_incompleteness_" + args.basis + "_" + args.quantity
        name += "_abs.png" if args.absolute and args.quantity == "exc" else ""
        name += ".png" if not name.endswith(".png") else ""
        args.output = os.path.join(PLOT_DIR, name)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    if args.input is None:
        name = "grid_incompleteness_" + args.basis + ".pkl"
        args.input = os.path.join(MY_DIR, name)
    with open(args.input, "rb") as f:
        results = pickle.load(f)
    plot(
        results,
        args.output,
        args.angular,
        args.absolute,
        args.quantity,
        args.basis,
    )


if __name__ == "__main__":
    main()
