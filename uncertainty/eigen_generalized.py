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

import numpy as np
from pluginplay import ModuleManager
from nwchemex import load_modules
from simde import GeneralizedEigenSolve
from tensorwrapper import Tensor


class GeneralizedEigenSolverTester:

    def __init__(self):
        self.mm = ModuleManager()
        load_modules(self.mm)
        self.solver = self.mm.at("Generalized eigensolve via Eigen")

    def solve_gen_eigenproblem(self, A, B, verify=True):
        """Solver function that can use pre-initialized manager and solver"""
        A = np.asarray(A, dtype=np.float64)
        B = np.asarray(B, dtype=np.float64)
        A_tensor, B_tensor = Tensor(A), Tensor(B)

        # Solve
        λ, V = map(
            np.array,
            self.solver.run_as(GeneralizedEigenSolve(), A_tensor, B_tensor))

        if verify:
            max_residual = self._verify_results(A, B, λ, V)
            print(f"Maximum residual norm: {max_residual:.2e}")

        return λ, V

    def _verify_results(self, A, B, λ, V, rtol=1e-8):
        """Internal verification method"""
        max_residual = 0
        for i, ev in enumerate(λ):
            res = A @ V[:, i] - ev * (B @ V[:, i])
            max_residual = max(max_residual, np.linalg.norm(res))
        return max_residual

    def analyze_uncertainty(self,
                            num_trials=100,
                            matrix_size=5,
                            noise_level=1e-10):
        """
        Analyze numerical uncertainty in generalized eigensolve by testing pertubed matricies

        Args:
            num_trials = Number of random matricies to test
            matrix_size = Size of the random matricies 
            noise_level: Magnitude of pertubations to add
        """

        eigenvalue_errors = []
        eigenvector_errors = []

        for _ in range(num_trials):
            #Generate random symmetric matrices
            A = np.random.rand(matrix_size, matrix_size)
            A = (A + A.T) / 2
            B = np.eye(matrix_size)

            A_perturbed = A + noise_level * np.random.randn(*A.shape)
            A_perturbed = (A_perturbed + A_perturbed.T) / 2

            λ, V = self.solve_gen_eigenproblem(A, B, False)
            λ_p, V_p = self.solve_gen_eigenproblem(A_perturbed, B, False)

            eigenvalue_errors.append(np.abs(λ - λ_p))
            eigenvector_errors.append([
                np.linalg.norm(V[:, i] - V_p[:, i]) for i in range(matrix_size)
            ])

        #results
        print("\n=== Uncertainty Analysis ===")
        print(
            f"Tested {num_trials} random {matrix_size}*{matrix_size} matricies"
        )
        print(f"Pertubation level: {noise_level:.1e}")

        print("\n Eigenvalue Error Statistics:")
        print(f"Mean absolute error: {np.mean(eigenvalue_errors):.2e}")
        print(f"Max absolute error: {np.max(eigenvalue_errors):.2e}")

        print("\n Eigenvector Error Statistics:")
        print(
            f"Mean angular differences: {np.mean(eigenvector_errors):.2e} radians"
        )
        print(
            f"Max angular difference: {np.max(eigenvector_errors):.2e} radians"
        )

    def run_standard_example(self):
        """Run a fixed example"""
        A = np.array([[1.0, 2.0], [2.0, 3.0]])
        B = np.eye(2)
        λ, V = self.solve_gen_eigenproblem(A, B)
        print("Eigenvalues:\n", λ)


def main():
    tester = GeneralizedEigenSolverTester()

    tester.analyze_uncertainty(num_trials=100,
                               matrix_size=5,
                               noise_level=1e-10)

    print("\n=== Standard Example ===")
    tester.run_standard_example()


if __name__ == "__main__":
    main()
