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


class TestQuadrature(unittest.TestCase):
    def check_fields(self, pt):
        self.assertEqual(len(pt.inputs()), 2)
        self.assertIn("BraKet", pt.inputs())
        self.assertIn("Grid", pt.inputs())
        self.assertEqual(len(pt.results()), 1)
        self.assertIn("tensor representation", pt.results())

    def test_xc_energy_quadrature(self):
        self.check_fields(scf.xc_energy_quadrature())

    def test_xc_potential_quadrature(self):
        self.check_fields(scf.xc_potential_quadrature())
