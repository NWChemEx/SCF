#pragma once
#include <simde/simde.hpp>
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

template<typename EnergyProp, typename FockProp, typename WfProp>
DECLARE_TEMPLATED_PROPERTY_TYPE(ConvergenceProp, EnergyProp, FockProp, WfProp);

template<typename EnergyProp, typename FockProp, typename WfProp>
TEMPLATED_PROPERTY_TYPE_INPUTS(ConvergenceProp, EnergyProp, FockProp, WfProp) {
    using wf_type = simde::type::rscf_wf;
    auto rv     = pluginplay::declare_input()
                    .add_field<EnergyProp>("New Energy")
                    .template add_field<EnergyProp>("Old Energy")
                    .template add_field<FockProp>("Fock Operator")
                    .template add_field<WfProp>("Wave Function")
                    .template add_field<double>("Energy Tolerance")
                    .template add_field<double>("Density Tolerance")
                    .template add_field<double>("Gradient Tolerance");

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

template<typename EnergyProp, typename FockProp, typename WfProp>
TEMPLATED_PROPERTY_TYPE_RESULTS(ConvergenceProp, EnergyProp, FockProp, WfProp) {
    auto rv     = pluginplay::declare_result().add_field<bool>("Convergence Status");
    rv.at("Convergence Status")
      .set_description("The molecule corresponding to the input string");
    return rv;
}

} // namespace simde

