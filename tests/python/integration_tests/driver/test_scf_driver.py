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

import nux
import simde
import parallelzone as pz
import chemist
import pluginplay as pp
import nwchemex as nwx
import numpy as np
import unittest


class TestSCFDriver(unittest.TestCase):

    def test_scf_driver(self):
        egy = self.mm.run_as(self.ao_energy, "SCF Driver", self.aos, self.sys)
        self.assertAlmostEqual(np.array(egy), -74.94208027122616, places=6)

    def test_scf_driver_diis(self):
        self.mm.change_input("Loop", "DIIS", True)
        egy = self.mm.run_as(self.ao_energy, "SCF Driver", self.aos, self.sys)
        self.assertAlmostEqual(np.array(egy), -74.94208027122616, places=6)

    def setUp(self):
        # Setup Module Manager
        self.mm = pp.ModuleManager(pz.runtime.RuntimeView())
        nux.load_modules(self.mm)
        nwx.load_modules(self.mm)

        # Set Submods
        self.mm.change_submod("SCF Driver", "Hamiltonian",
                              "Born-Oppenheimer Approximation")
        self.mm.change_submod("SCF integral driver", "Fundamental matrices",
                              "AO integral driver")
        self.mm.change_submod("Diagonalization Fock update",
                              "Overlap matrix builder", "Overlap")
        self.mm.change_submod("Loop", "Overlap matrix builder", "Overlap")

        # Property Types
        self.mol = simde.MoleculeFromString()
        self.basis_set = simde.MolecularBasisSet()
        self.ao_energy = simde.AOEnergy()

        # Inputs
        self.water = self.mm.at("NWX Molecules").run_as(self.mol, "water")
        self.aos = self.mm.at('STO-3G').run_as(self.basis_set, self.water)
        self.sys = chemist.ChemicalSystem(self.water)
