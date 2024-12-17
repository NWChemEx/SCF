// #include "matrix_builders.hpp"

// namespace scf::matrix_builders {

// using pt             = simde::aos_f_e_aos;
// using op_base        = chemist::qm_operator::OperatorBase;
// using aos            = simde::type::aos;
// using aos_squared    = simde::type::aos_squared;
// using braket_type_2c = simde::type::braket<aos, op_base aos>;
// using braket_type_4c = simde::type::braket<aos_squared, op_base,
// aos_squared>; using pt_2c          = simde::EvaluateBraKet<braket_type_2c>;
// using pt_4c          = simde::EvaluateBraKet<braket_type_4c>;

// const auto desc = R"(
// )";

// MODULE_CTOR(Fock) {
//     description(desc);
//     satisfies_property_type<pt>();
//     add_submodule<pt_2c>("Two center evaluator");
//     add_submodule<pt_4c>("Four center evaluator");
// }

// } // namespace scf::matrix_builders