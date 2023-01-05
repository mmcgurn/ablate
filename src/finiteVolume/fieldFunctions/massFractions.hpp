#ifndef ABLATELIBRARY_FIELD_FUNCTION_MASSFRACTIONS_HPP
#define ABLATELIBRARY_FIELD_FUNCTION_MASSFRACTIONS_HPP

#include <eos/eos.hpp>
#include "compressibleFlowState.hpp"
#include "domain/fieldMathFunction.hpp"

namespace ablate::finiteVolume::fieldFunctions {

/**
 * Class that allows you to under specify the yi and assume zero for all others
 */
class MassFractions : public ablate::domain::FieldMathFunction {
   private:
    std::vector<std::shared_ptr<mathFunctions::MathFunction>> massFractionFunctions;
    std::vector<std::shared_ptr<ablate::domain::FieldMathFunction>> massFractionFieldFunctions;

    static PetscErrorCode ComputeYiFunction(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar* u, void* ctx);

   public:
    explicit MassFractions(std::shared_ptr<ablate::eos::EOS> eos, std::vector<std::shared_ptr<domain::FieldMathFunction>> massFractionFunctions);
};

}  // namespace ablate::finiteVolume::fieldFunctions
#endif  // ABLATELIBRARY_FIELD_SOLUTION_EULER_HPP
