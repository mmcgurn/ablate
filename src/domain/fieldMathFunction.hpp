#ifndef ABLATELIBRARY_FIELDMATHFUNCTION_HPP
#define ABLATELIBRARY_FIELDMATHFUNCTION_HPP
#include <memory>
#include <string>
#include "domain/region.hpp"
#include "fieldFunction.hpp"
#include "mathFunctions/mathFunction.hpp"

namespace ablate::domain {
class FieldMathFunction : public FieldFunction {
   protected:
    //! The field name to which to apply the vector
    const std::string fieldName;

    //! The math
    const std::shared_ptr<mathFunctions::MathFunction> solutionField;
    const std::shared_ptr<ablate::domain::Region> region;

   public:
    FieldMathFunction(std::string fieldName, std::shared_ptr<mathFunctions::MathFunction> solutionField, std::shared_ptr<ablate::domain::Region> region = nullptr);
    virtual ~FieldMathFunction() = default;

    /**
     * The name of the field
     * @return
     */
    const std::string& GetName() const { return fieldName; }

    /**
     * Project the math function onto this global vector
     * @param globVec
     * @param time
     */
    void ProjectField(Vec globVec, PetscReal time = 0.0) override {}

    bool HasSolutionField() const { return solutionField != nullptr; }
    mathFunctions::MathFunction& GetSolutionField() { return *solutionField; }

    std::shared_ptr<mathFunctions::MathFunction> GetFieldFunction() const { return solutionField; }
};
}  // namespace ablate::domain
#endif  // ABLATELIBRARY_FIELDMATHFUNCTION_HPP
