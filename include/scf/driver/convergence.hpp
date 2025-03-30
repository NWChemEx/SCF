#pragma once
#include <pluginplay/pluginplay.hpp>
#include <simde/types.hpp>


namespace scf {

DECLARE_PROPERTY_TYPE(ConvergenceProp);

PROPERTY_TYPE_INPUTS(ConvergenceProp) {
    auto rv     = pluginplay::declare_input()
                    .add_field<simde::type::tensor>("New Energy")
                    .template add_field<simde::type::tensor>("Old Energy")
                    .template add_field<chemist::DecomposableDensity<chemist::Electron>>("New Rho")
                    .template add_field<chemist::DecomposableDensity<chemist::Electron>>("Old Rho")
                    .template add_field<simde::type::tensor>("P_new")
                    .template add_field<simde::type::tensor>("S")
                    .template add_field<simde::type::fock>("Fock Operator")
                    .template add_field<chemist::wavefunction::AOs>("Wave Function")
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

PROPERTY_TYPE_RESULTS(ConvergenceProp) {
    auto rv     = pluginplay::declare_result().add_field<bool>("Convergence Status");
    rv.at("Convergence Status")
      .set_description("The molecule corresponding to the input string");
    return rv;
}

} // namespace simde

