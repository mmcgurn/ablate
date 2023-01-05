#ifndef ABLATELIBRARY_EXACTFUNCTION_HPP
#define ABLATELIBRARY_EXACTFUNCTION_HPP
#include <memory>
#include <string>
#include "fieldMathFunction.hpp"
#include "mathFunctions/mathFunction.hpp"

namespace ablate::domain {
class ExactFunction : public FieldMathFunction {
   private:
    const std::shared_ptr<mathFunctions::MathFunction> timeDerivative;

   public:
    ExactFunction(std::string fieldName, std::shared_ptr<mathFunctions::MathFunction> solutionField, std::shared_ptr<mathFunctions::MathFunction> timeDerivative = {});
    virtual ~ExactFunction() = default;

    const std::string& GetName() const { return fieldName; }

    bool HasSolutionField() const { return solutionField != nullptr; }
    mathFunctions::MathFunction& GetSolutionField() { return *solutionField; }

    bool HasTimeDerivative() const { return timeDerivative != nullptr; }
    mathFunctions::MathFunction& GetTimeDerivative() { return *timeDerivative; }

    std::shared_ptr<mathFunctions::MathFunction> GetFieldFunction() const { return solutionField; }
    std::shared_ptr<mathFunctions::MathFunction> GetTimeDerivativeFunction() const { return timeDerivative; }
};
}  // namespace ablate::domain
#endif  // ABLATELIBRARY_EXACTFUNCTION_HPP
