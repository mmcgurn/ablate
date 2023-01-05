#include "exactFunction.hpp"
#include "registrar.hpp"

ablate::domain::ExactFunction::ExactFunction(std::string fieldName, std::shared_ptr<ablate::mathFunctions::MathFunction> solutionField,
                                             std::shared_ptr<ablate::mathFunctions::MathFunction> timeDerivative)
    : FieldMathFunction(fieldName, solutionField), timeDerivative(timeDerivative) {}

REGISTER_DEFAULT(ablate::domain::ExactFunction, ablate::domain::ExactFunction, "a field description that can be used for initialization or exact solution ",
                 ARG(std::string, "fieldName", "the field name"), ARG(ablate::mathFunctions::MathFunction, "field", "the math function used to describe the field"),
                 OPT(ablate::mathFunctions::MathFunction, "timeDerivative", "the math function used to describe the field time derivative"));
