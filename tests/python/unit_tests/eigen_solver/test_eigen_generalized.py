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
import pluginplay as pp
import parallelzone as pz
import scf
import simde
import tensorwrapper as tw
import numpy as np


class TestEigenGeneralized(unittest.TestCase):

    def test_eigen_generalized(self):
        # Property Type
        pt = simde.GeneralizedEigenSolve()

        # Inputs
        A = tw.Tensor(np.array([[1.0, 2.0], [2.0, 3.0]]))
        B = tw.Tensor(np.array([[1.0, 0.0], [0.0, 1.0]]))

        # Run module
        mod_name = "Generalized eigensolve via Eigen"
        values, vector = self.mm.run_as(pt, mod_name, A, B)
        test_values = np.array(values)
        self.assertAlmostEqual(test_values[0], -0.236068, places=6)
        self.assertAlmostEqual(test_values[1], 4.236068, places=6)

    def setUp(self):
        self.mm = pp.ModuleManager(pz.runtime.RuntimeView())
        scf.load_modules(self.mm)
