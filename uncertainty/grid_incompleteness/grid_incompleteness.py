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

"""Grid incompleteness of SVWN3 for He, H2, and water.

For each molecule, an SVWN3 SCF (STO-3G by default; see --basis) is run on
every combination of the requested radial and angular grid sizes. For each
grid, the total energy, the XC energy, and the XC potential (AO basis) are
recorded. The XC energy and potential are those of the converged density.

The XC energy and potential are captured by ``Recorder`` modules, which wrap
the SVWN3 quadrature modules and store each result they compute. Every SCF
iteration overwrites the previous one, so once the SCF finishes the stored
values are those of the final (converged) density.

Results are kept in a dictionary keyed by ``(molecule, radial, angular)``,
which is pickled after each calculation (one file per basis set). When the
script starts, it reads the pickled dictionary (if any) and skips the
calculations it already holds.
"""

import argparse
import os
import pickle
import time

import chemcache
import chemist
import integrals
import numpy as np
import nux
import parallelzone as pz
import pluginplay as pp
import scf
import simde

#: Radial (Mura-Knowles) grid sizes, coarse to fine
DEFAULT_RADIAL = [10, 20, 35, 50, 75]

#: Angular (Lebedev-Laikov) grid sizes, coarse to fine
DEFAULT_ANGULAR = [14, 26, 50, 110, 302]

MOLECULES = ["he", "h2", "water"]

#: Basis set, by its ChemCache module key
DEFAULT_BASIS = "sto-3g"

MY_DIR = os.path.dirname(os.path.realpath(__file__))


def default_output(basis):
    """The pickle the results for @p basis are kept in by default."""
    return os.path.join(MY_DIR, "grid_incompleteness_" + basis + ".pkl")


class Recorder(pp.ModuleBase):
    """Runs the wrapped module and records the result it computes.

    The inputs are forwarded to the wrapped module without being unwrapped,
    since not every BraKet the quadratures take is exported to Python.
    """

    def __init__(self, pt, record, key):
        pp.ModuleBase.__init__(self)
        self.satisfies_property_type(pt)
        self.description("Records the result of the wrapped module")
        self.add_submodule(pt, "Wrapped")
        self.pt = pt
        self.record = record
        self.key = key

    def run_(self, inputs, submods):
        rv = submods["Wrapped"].value().run(inputs)
        (tensor,) = self.pt.unwrap_results(rv)
        self.record[self.key] = np.array(tensor)
        return rv


def make_module_manager(runtime, record, n_radial, n_angular, tolerances):
    """Makes a ModuleManager set up for an SVWN3 SCF on the requested grid.

    A new ModuleManager is made for each calculation so that no memoized
    result can carry over from one grid to the next.
    """
    mm = pp.ModuleManager(runtime)
    mm.load_plugin(scf.SCF())
    mm.load_plugin(integrals.Integrals())
    mm.load_plugin(chemcache.ChemCache())
    mm.load_plugin(nux.NUX())
    mm.wire_defaults()

    # SVWN3. SCF dispatches SVWN3 to the native quadratures by default.
    mm.change_input("SCF Driver", "DFT Functional", "SVWN3")

    # Record the XC energy and potential
    recorders = [
        ("SVWN3 Energy", scf.xc_energy_quadrature(), "E_xc"),
        ("SVWN3 Potential", scf.xc_potential_quadrature(), "V_xc"),
    ]
    for xc_key, pt, record_key in recorders:
        quad_key = xc_key + " Quadrature"
        recorder_key = "Recording " + quad_key
        mm.add_module(recorder_key, Recorder(pt, record, record_key))
        mm.change_submod(recorder_key, "Wrapped", quad_key)
        mm.change_submod(xc_key, "Quadrature", recorder_key)

    # These modules get new inputs (a new density) every SCF iteration, so
    # their caches never hit. Looking inputs up in them is expensive, though:
    # each lookup compares the grid point by point against every cached
    # entry. Turning memoization off does not change the results and makes
    # the SCF about 2.6x faster (water, 20 x 110 grid). The grid and the AOs
    # on the grid do not change between iterations, so they stay memoized.
    per_iteration = [
        "SVWN3 Energy",
        "SVWN3 Potential",
        "SVWN3 Energy Quadrature",
        "SVWN3 Potential Quadrature",
        "Recording SVWN3 Energy Quadrature",
        "Recording SVWN3 Potential Quadrature",
        "Density2Grid",
    ]
    for key in per_iteration:
        mm.at(key).turn_off_memoization()

    # The grid
    mm.change_input("Grid From IntegratorXX", "Radial Size", n_radial)
    mm.change_input("Grid From IntegratorXX", "Angular Size", n_angular)

    # Converge tightly enough that the density error does not mask the grid
    e_tol, d_tol, g_tol, max_iter = tolerances
    mm.change_input("Loop", "energy tolerance", e_tol)
    mm.change_input("Loop", "density tolerance", d_tol)
    mm.change_input("Loop", "gradient tolerance", g_tol)
    mm.change_input("Loop", "max iterations", max_iter)
    return mm


def make_molecule(mm, name):
    """Makes the molecule named @p name (coordinates in bohr)."""
    if name == "water":
        return mm.at("NWX Molecules").run_as(simde.MoleculeFromString(), name)

    molecule = chemist.Molecule()
    if name == "he":
        molecule.push_back(chemist.Atom("He", 2, 7296.2971, 0.0, 0.0, 0.0))
    elif name == "h2":
        # Same geometry as the SCF tests' H2
        molecule.push_back(chemist.Atom("H", 1, 1837.4260, 0.0, 0.0, 0.0))
        molecule.push_back(chemist.Atom("H", 1, 1837.4260, 0.0, 0.0, 1.3984))
    else:
        raise ValueError("Unknown molecule: " + name)
    return molecule


def run_calculation(runtime, name, basis, n_radial, n_angular, tolerances):
    """Runs one SVWN3 SCF and returns what was recorded."""
    record = {}
    mm = make_module_manager(runtime, record, n_radial, n_angular, tolerances)
    molecule = make_molecule(mm, name)
    aos = mm.at(basis).run_as(simde.MolecularBasisSet(), molecule)
    sys = chemist.ChemicalSystem(molecule)

    start = time.perf_counter()
    e_total = mm.run_as(simde.AOEnergy(), "SCF Driver", aos, sys)
    wall_time = time.perf_counter() - start

    return {
        "n_points": molecule.size() * n_radial * n_angular,
        "E_total": float(np.array(e_total)),
        "E_xc": float(record["E_xc"]),
        "V_xc": record["V_xc"],
        "wall time": wall_time,
    }


def load_results(path):
    """Returns the pickled results at @p path, or an empty dictionary."""
    if not os.path.exists(path):
        return {}
    with open(path, "rb") as f:
        return pickle.load(f)


def save_results(results, path):
    """Pickles @p results to @p path, without clobbering it if interrupted."""
    tmp_path = path + ".tmp"
    with open(tmp_path, "wb") as f:
        pickle.dump(results, f)
    os.replace(tmp_path, path)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--molecules", nargs="+", choices=MOLECULES, default=MOLECULES
    )
    parser.add_argument(
        "--radial", nargs="+", type=int, default=DEFAULT_RADIAL
    )
    parser.add_argument(
        "--angular", nargs="+", type=int, default=DEFAULT_ANGULAR
    )
    parser.add_argument(
        "--basis",
        default=DEFAULT_BASIS,
        help="ChemCache key of the basis set (e.g., aug-cc-pvdz)",
    )
    parser.add_argument(
        "--output",
        help="Defaults to grid_incompleteness_<basis>.pkl next to this script",
    )
    # Tight so that the density error does not mask the grid error. The change
    # in the density stalls at a few 1E-10 on some grids (e.g., water 75x14),
    # so the density tolerance can not go lower. N.B. "gradient tolerance"
    # applies to the SQUARED norm of the orbital gradient.
    parser.add_argument("--energy-tol", type=float, default=1.0e-12)
    parser.add_argument("--density-tol", type=float, default=1.0e-9)
    parser.add_argument("--gradient-tol", type=float, default=1.0e-16)
    parser.add_argument("--max-iterations", type=int, default=50)
    return parser.parse_args()


def main():
    args = parse_args()
    tolerances = (
        args.energy_tol,
        args.density_tol,
        args.gradient_tol,
        args.max_iterations,
    )
    if args.output is None:
        args.output = default_output(args.basis)
    results = load_results(args.output)
    runtime = pz.runtime.RuntimeView()

    header = "{:>6} {:>7} {:>8} {:>9} {:>20} {:>20} {:>9}"
    row = "{:>6} {:>7} {:>8} {:>9} {:>20.12f} {:>20.12f} {:>9.1f}"
    print(
        header.format(
            "mol", "radial", "angular", "points", "E_total", "E_xc", "time (s)"
        )
    )
    for name in args.molecules:
        for n_radial in args.radial:
            for n_angular in args.angular:
                key = (name, n_radial, n_angular)
                if key in results:
                    status = "(cached)"
                else:
                    results[key] = run_calculation(
                        runtime,
                        name,
                        args.basis,
                        n_radial,
                        n_angular,
                        tolerances,
                    )
                    save_results(results, args.output)
                    status = ""
                r = results[key]
                print(
                    row.format(
                        name,
                        n_radial,
                        n_angular,
                        r["n_points"],
                        r["E_total"],
                        r["E_xc"],
                        r["wall time"],
                    ),
                    status,
                    flush=True,
                )


if __name__ == "__main__":
    main()
