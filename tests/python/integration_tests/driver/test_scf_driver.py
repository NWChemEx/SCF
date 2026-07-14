# Copyright 2025 NWChemEx-Project
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

import unittest

import chemist
import numpy as np
import nux
import nwchemex as nwx
import parallelzone as pz
import pluginplay as pp
import simde


class TestSCFDriver(unittest.TestCase):
    def test_scf_driver(self):
        self.mm.change_input("Loop", "DIIS", False)
        egy = self.mm.run_as(self.ao_energy, "SCF Driver", self.aos, self.sys)
        self.assertAlmostEqual(np.array(egy), -74.94208027122616, places=6)

    def test_scf_driver_diis(self):
        egy = self.mm.run_as(self.ao_energy, "SCF Driver", self.aos, self.sys)
        self.assertAlmostEqual(np.array(egy), -74.94208027122616, places=6)

    def test_dft_driver(self):
        func = chemist.qm_operator.xc_functional.PBE
        RKS_op = "Restricted Kohn-Sham Op"
        rks_op = "Restricted One-Electron Kohn-Sham Op"
        self.mm.change_input(RKS_op, "XC Potential", func)
        self.mm.change_input(rks_op, "XC Potential", func)
        self.mm.change_submod("Loop", "One-electron Fock operator", rks_op)
        self.mm.change_submod("Loop", "Fock operator", RKS_op)
        self.mm.change_submod("Core guess", "Build Fock Operator", rks_op)
        egy = self.mm.run_as(self.ao_energy, "SCF Driver", self.aos, self.sys)
        self.assertAlmostEqual(np.array(egy), -75.22989870343218, places=6)

    def setUp(self):
        # Setup Module Manager
        self.mm = pp.ModuleManager(pz.runtime.RuntimeView())
        nux.load_modules(self.mm)
        nwx.load_modules(self.mm)

        # integrals.set_defaults(mm) isn't exposed to Python (EXPORT_PLUGIN
        # only binds load_modules), so replicate its submodule wiring here --
        # see integrals::libint::set_defaults, integrals::ao_integrals::
        # set_defaults, and integrals::utils::set_defaults in
        # cxx/src/integrals/{libint,ao_integrals,utils}/*.{cpp,hpp}.
        self.mm.change_submod(
            "CauchySchwarz Estimator", "Decontract Basis Set",
            "Decontract Basis Set"
        )
        self.mm.copy_module("ERI4", "Benchmark ERI4")
        self.mm.change_input("Benchmark ERI4", "Threshold", 1.0e-16)
        self.mm.change_submod("CauchySchwarz Estimator", "ERI4", "Benchmark ERI4")
        self.mm.change_submod("Analytic Error", "ERI4s", "Benchmark ERI4")
        self.mm.change_submod(
            "Raw Primitive ERI4", "Decontract Basis Set", "Decontract Basis Set"
        )
        self.mm.change_submod(
            "Primitive Contractor ERI4", "Raw Primitive ERI4",
            "Raw Primitive ERI4"
        )
        self.mm.change_submod(
            "Primitive Contractor ERI4", "Primitive Normalization",
            "Primitive Normalization"
        )
        self.mm.change_submod(
            "AO integral driver", "Coulomb matrix", "Four center J builder"
        )
        self.mm.change_submod(
            "AO integral driver", "Exchange matrix", "Four center K builder"
        )
        # These four calls live directly in integrals::set_defaults (not in
        # one of the libint/ao_integrals/utils sub-namespace helpers above).
        self.mm.change_submod("AO integral driver", "Kinetic", "Kinetic")
        self.mm.change_submod(
            "AO integral driver", "Electron-Nuclear attraction", "Nuclear"
        )
        self.mm.change_submod("Four center J builder", "Four-center ERI", "ERI4")
        self.mm.change_submod("Four center K builder", "Four-center ERI", "ERI4")
        self.mm.change_submod(
            "Density Fitting Integral", "Three-center ERI", "ERI3"
        )
        self.mm.change_submod("Coulomb Metric", "Two-center ERI", "ERI2")
        self.mm.change_submod(
            "Density Fitted J builder", "DF ERI", "Density Fitting Integral"
        )
        self.mm.change_submod(
            "Density Fitted K builder", "DF ERI", "Density Fitting Integral"
        )
        self.mm.change_submod(
            "Density Fitting Integral", "Coulomb Metric", "Coulomb Metric"
        )
        self.mm.change_submod("UQ Driver", "ERIs", "ERI4")
        self.mm.change_submod("UQ Driver", "ERI Error", "Primitive Error Model")
        self.mm.change_submod("UQ Atom Symm Blocked Driver", "ERIs", "ERI4")
        self.mm.change_submod(
            "UQ Atom Symm Blocked Driver", "ERI Error", "Primitive Error Model"
        )
        self.mm.change_submod(
            "Screen Primitive Pairs", "Primitive Pair Estimator",
            "Black Box Primitive Pair Estimator"
        )

        # Set Submods
        self.mm.change_submod(
            "SCF Driver", "Hamiltonian", "Born-Oppenheimer Approximation"
        )
        self.mm.change_submod(
            "SCF integral driver", "Fundamental matrices", "AO integral driver"
        )
        self.mm.change_submod(
            "Diagonalization Fock update", "Overlap matrix builder", "Overlap"
        )
        self.mm.change_submod("Loop", "Overlap matrix builder", "Overlap")
        self.mm.change_submod("SAD guess", "SAD Density", "sto-3g SAD density")

        # Property Types
        self.mol = simde.MoleculeFromString()
        self.basis_set = simde.MolecularBasisSet()
        self.ao_energy = simde.AOEnergy()

        # Inputs
        self.water = self.mm.at("NWX Molecules").run_as(self.mol, "water")
        self.aos = self.mm.at("STO-3G").run_as(self.basis_set, self.water)
        self.sys = chemist.ChemicalSystem(self.water)
