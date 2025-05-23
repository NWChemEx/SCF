# Copyright 2023 NWChemEx-Project
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

from pluginplay import ModuleManager
import scf
from simde import AOEnergy, MoleculeFromString, MolecularBasisSet
import nwchemex
from chemist import ChemicalSystem
import unittest


class TestTAMMSCF(unittest.TestCase):

    def test_4cHF(self):
        if self.has_TAMM:
            self.mm.change_input(self.key, 'molecule_name', 'water')
            egy = self.mm.run_as(AOEnergy(), self.key, self.aos, self.cs)
            self.assertAlmostEqual(egy, -74.3670617803483, places=6)

    def test_dft(self):
        if self.has_TAMM:
            self.mm.change_input(self.key, 'xc_type', ["pbe0"])
            self.mm.change_input(self.key, 'molecule_name', 'water')
            egy = self.mm.run_as(AOEnergy(), self.key, self.aos, self.cs)
            self.assertAlmostEqual(egy, -74.81168986385825, places=6)

    def setUp(self):
        self.mm = ModuleManager()
        scf.load_modules(self.mm)
        self.key = 'SCF Energy via TAMM'
        self.has_TAMM = self.mm.count(self.key) != 0
        if self.has_TAMM:
            nwchemex.load_modules(self.mm)
            mol = self.mm.run_as(MoleculeFromString(), "NWX Molecules",
                                 "water")
            self.aos = self.mm.run_as(MolecularBasisSet(), 'sto-3g', mol)
            self.cs = ChemicalSystem(mol)
