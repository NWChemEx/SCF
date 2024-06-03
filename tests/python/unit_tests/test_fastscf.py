#!/Users/jwaldrop/venvs/nwx/bin/python
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
from fastscf import load_modules
from simde import AOEnergy, MoleculeFromString, MolecularBasisSet
import chemcache as ccache
from molecules import make_h2
import parallelzone as pz
from chemist import ChemicalSystem
import unittest
import os
import sys


class TestFastSCF(unittest.TestCase):

    def test_scf(self):
        mol_name = "water"
        mol = self.mm.run_as(MoleculeFromString(), "NWX Molecules", mol_name)
        cs = ChemicalSystem(mol)

        basis_name = "sto-3g"
        aos = self.mm.run_as(MolecularBasisSet(), basis_name, mol)

        key = 'FastSCF Energy'
        egy = self.mm.run_as(AOEnergy(), key, aos, cs)
        self.assertAlmostEqual(egy, -1.094184522864, places=5)

    def setUp(self):
        self.mm = ModuleManager()
        ccache.load_modules(self.mm)
        load_modules(self.mm)


if __name__ == '__main__':
    # fastscf.tamm_initialize(argc, argv)
    rv = pz.runtime.RuntimeView()

    my_dir = os.path.dirname(os.path.realpath(__file__))

    loader = unittest.TestLoader()
    tests = loader.discover(my_dir)
    testrunner = unittest.runner.TextTestRunner()
    ret = not testrunner.run(tests).wasSuccessful()

    # fastscf.tamm_finalize()

    sys.exit(ret)
