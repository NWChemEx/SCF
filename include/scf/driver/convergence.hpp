#pragma once
#include "simde/simde.hpp"
#include <pluginplay/pluginplay.hpp>

namespace scf {

/** @brief The property type for modules that return a molecule from a string
 *         input.
 *
 * The usage of this property type can be pretty varied. Modules that satisfy
 * this property type could potentially produce a molecule from a string with
 * XYZ inputs, a SMILES string, or a molecule name outright.
 *
 */

using simde::type::electronic_hamiltonian;
using simde::type::hamiltonian;
using simde::type::op_base_type;

template<typename WfType>
using egy_pt = simde::eval_braket<WfType, hamiltonian, WfType>;

template<typename WfType>
using elec_egy_pt = simde::eval_braket<WfType, electronic_hamiltonian, WfType>;

template<typename WfType>
using pt = simde::Optimize<egy_pt<WfType>, WfType>;

template<typename WfType>
using update_pt = simde::UpdateGuess<WfType>;

using density_t = simde::type::decomposable_e_density;

using fock_pt = simde::FockOperator<density_t>;

using density_pt = simde::aos_rho_e_aos<simde::type::cmos>;

using v_nn_pt = simde::charge_charge_interaction;

using fock_matrix_pt = simde::aos_f_e_aos;
using s_pt           = simde::aos_s_e_aos;
 
DECLARE_PROPERTY_TYPE(ConvergenceProp);

PROPERTY_TYPE_INPUTS(ConvergenceProp) {
    using wf_type = simde::type::rscf_wf;
    auto rv     = pluginplay::declare_input()
                    .add_field<BraKet::result_type>("New Energy")
                    .add_field<elec_egy_pt<wf_type>>("Old Energy")
                    .add_field<elec_egy_pt<wf_type>>("Fock Operator")
                    .add_field<elec_egy_pt<wf_type>>("Wave Function")
                    .add_field<double>("Energy Tolerance")
                    .add_field<double>("Density Tolerance")
                    .add_field<double>("Gradient Tolerance");

    rv.at("New Energy").set_description(
      "The string identifying the desired molecule");
    rv.at("Old Energy").set_description(
      "The string identifying the desired molecule");
    rv.at("Fock Operator").set_description(
      "The string identifying the desired molecule");
    rv.at("Wave Function").set_description(
      "The string identifying the desired molecule");
    rv.at("Energy Tolerance").set_description(
      "The string identifying the desired molecule");
    rv.at("Density Tolerance").set_description(
      "The string identifying the desired molecule");
    rv.at("Gradient Tolerance").set_description(
      "The string identifying the desired molecule");
    return rv;
}

PROPERTY_TYPE_RESULTS(ConvergenceProp) {
    auto rv     = pluginplay::declare_result().add_field<bool>("Convergence Status");
    rv.at("Convergence Status")
      .set_description("The molecule corresponding to the input string");
    return rv;
}

} // namespace simde

