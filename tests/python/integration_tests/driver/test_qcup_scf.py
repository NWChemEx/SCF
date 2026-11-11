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

import unittest

import chemcache
import chemist
import integrals
import numpy as np
import nux
import parallelzone as pz
import pluginplay as pp
import scf
import simde


@unittest.skipUnless(scf.sigma_enabled, "QCUP SCF requires Sigma")
class TestQCUPSCF(unittest.TestCase):
    def test_qcup_scf(self):
        self.mm.change_input("QCUP SCF", "UQ Type", "thresholded affine")
        self.mm.change_submod(
            "QCUP SCF", "ERI Error Model", "Primitive Error Model"
        )
        self.mm.change_input("QCUP SCF", "ERI Error Estimate", "Coarse")
        egy = self.mm.run_as(self.ao_energy, "QCUP SCF", self.aos, self.sys)
        self.assertIsNotNone(egy)

    def test_inactive_settings(self):
        self.mm.change_input("QCUP SCF", "ERI Error Estimate", "Coarse")
        self.mm.run_as(self.ao_energy, "QCUP SCF", self.aos, self.sys)
        inactive = self.mm.at("QCUP SCF").inactive_settings()
        self.assertTrue(any("Primitive Error Model" in x for x in inactive))

    def test_ordinary_scf_is_unaffected(self):
        self.mm.run_as(self.ao_energy, "QCUP SCF", self.aos, self.sys)
        egy = self.mm.run_as(self.ao_energy, "SCF Driver", self.aos, self.sys)
        self.assertAlmostEqual(np.array(egy), -74.94208027122616, places=6)

    def setUp(self):
        # Setup Module Manager
        self.mm = pp.ModuleManager(pz.runtime.RuntimeView())
        self.mm.load_plugin(scf.SCF())
        self.mm.load_plugin(integrals.Integrals())
        self.mm.load_plugin(chemcache.ChemCache())
        self.mm.load_plugin(nux.NUX())
        self.mm.wire_defaults()

        # Property Types
        self.mol = simde.MoleculeFromString()
        self.basis_set = simde.MolecularBasisSet()
        self.ao_energy = simde.AOEnergy()

        # Inputs
        self.water = self.mm.at("NWX Molecules").run_as(self.mol, "water")
        self.aos = self.mm.at("STO-3G").run_as(self.basis_set, self.water)
        self.sys = chemist.ChemicalSystem(self.water)
