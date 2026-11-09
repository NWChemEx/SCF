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

import scf


class TestMOUncertainty(unittest.TestCase):
    def test_fields(self):
        pt = scf.MOUncertainty()
        inputs = [
            "Fock matrix",
            "Overlap matrix",
            "Orbital energies",
            "Orbital coefficients",
            "Density change",
        ]
        self.assertEqual(len(pt.inputs()), len(inputs))
        for label in inputs:
            self.assertIn(label, pt.inputs())
        self.assertEqual(len(pt.results()), 1)
        self.assertIn("Orbital coefficients", pt.results())


class TestConvergenceUncertainty(unittest.TestCase):
    def test_fields(self):
        pt = scf.ConvergenceUncertainty()
        self.assertEqual(len(pt.inputs()), 2)
        self.assertIn("Energy", pt.inputs())
        self.assertIn("Energy change", pt.inputs())
        self.assertEqual(len(pt.results()), 1)
        self.assertIn("Energy", pt.results())
