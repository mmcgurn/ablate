#include "fieldMathFunction.hpp"
#include "registrar.hpp"

ablate::domain::FieldMathFunction::FieldMathFunction(std::string fieldName, std::shared_ptr<mathFunctions::MathFunction> solutionField,
                                                    std::shared_ptr<ablate::domain::Region> region)
    :  fieldName(fieldName), solutionField(solutionField), region(region) {}

REGISTER_DERIVED(ablate::domain::FieldFunction, ablate::domain::FieldMathFunction);

REGISTER_DEFAULT(ablate::domain::FieldMathFunction, ablate::domain::FieldMathFunction, "a field description that can be used for initialization or exact solution ",
                 ARG(std::string, "fieldName", "the field name"), ARG(ablate::mathFunctions::MathFunction, "field", "the math function used to describe the field"),
                 OPT(ablate::domain::Region, "region", "A subset of the domain to apply the field function"));
