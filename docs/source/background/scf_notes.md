# SCF Procedure

## SCF Driver

1. Obtain the AO parameters as well as the chemical system (inputs)
2. Assign the AO parameters as atomic orbitals
3. Create the Hamitonian from the Chemical System (Submod "Hamiltonian" used, currently only the BO Approximation from NUX)
    1. For the BO Approximation Hamiltonian:
        1. From the Chemical System, the Electrons and Nuclei are used to form the Hamiltonian
        2. Each electron gets a 1 electron kinetic energy term
        3. Each nuclei gets an Electron-Nuclei potential energy term
        4. If more than 1 electron, an Electron-Electron potential energy repulson term is added
        5. If more than 1 nuclei, a Nuclei-Nuclei potential energy repulsion term is added
        6. In order, the Electron kinetic energy, Electron-nucleus attraction energy, and Electron-electron repulsion energy are
           formed to create the Electronic Hamiltonian
        7. The Nuclear-nuclear repulsion energy is also created
        8. The BO Approximation returns the constructed Hamiltonian, H.
4. With the guess submod, (of type InitialGuess), the inputs are the Hamiltonian and the AOs, and the return value is the
   initial wavefunction
    1. The main module used is the "Core guess" module, which first builds the fock operator with 0 density
    2. I am not sure on what happens during lines 78-92:
        ```cpp
        simde::type::cmos cmos(tensor{}, aos, tensor{});
        NElectronCounter visitor;
        H.visit(visitor);
        auto n_electrons = visitor.n_electrons;
        if(n_electrons % 2 != 0)
            throw std::runtime_error("Assumed even number of electrons");

        typename rscf_wf::orbital_index_set_type occs;
        using value_type = typename rscf_wf::orbital_index_set_type::value_type;
        for(value_type i = 0; i < n_electrons / 2; ++i) occs.insert(i);

        rscf_wf zero_guess(occs, cmos);
        auto& update_mod = submods.at("Guess updater");
        const auto& Psi0 = update_mod.run_as<update_pt>(f, zero_guess);```

    3. From there, Psi0 is returned

### SCF Loop

5. The work horse of the module, the submod "Optimizer", uses `scf_loop.cpp` and the module `Loop` to optimize
   the wavefunction. It takes in a Bracket via <Psi0 | H | Psi0> = H00, and Psi0. The property type it satisfies is
    ```cpp
    using pt = simde::Optimize<egy_pt<WfType>, WfType>;
    using egy_pt = simde::eval_braket<WfType, hamiltonian, WfType>;
    using wf_type = simde::type::rscf_wf;
    satisfies_property_type<pt<wf_type>>;
    ```
    and the return result type is the `WfType`, or in this case the `simde::type::rscf_wf`.
    
   1. Step one is to get the nuclear-nuclear repulsion via lines 125-141
      1. This is a little messy, so I'll skip the explanation.
   2. Compute the overlap matrix, S
      1. This runs the S_mod, which is the "Overlap" module provided by the "Integrals" plugin, creating the overlap matrix.
   3. Build the density matrix
   4. Start the scf process
      1. The fock matrix is build, starting with the one-electron fock builder, then the Fock builder, which who knows what is going on there.
      2. The new wavefunction is build from the new fock matrix and the old wavefunction
      3. the new density is built from the new wavefunction
      4. The new electronic hamiltonian is built from the H_core and the new Fock matrix
      5. the H_00 BraKet is built, and then used to calculate the new electronic energy
      6. We then check for convergence:
         1. We calculate delta_e for the change in energy, the delta_p for the change in density, and then the
            gradient norm (FPS - SPF).
         2. check to see if they have hit convergence
      7. energy total is then calculated by adding the electronic energy to the nuclear energy calculated earlier. 

### Return results
6. From here, the energy that is of type Tensor and the optimized wavefunction is returned. 
