#ifndef ABLATELIBRARY_FIELDFUNCTION_HPP
#define ABLATELIBRARY_FIELDFUNCTION_HPP
#include <memory>
#include <string>
#include "domain/region.hpp"
#include "mathFunctions/mathFunction.hpp"

namespace ablate::domain {
/**
 * Interface to describe all projection of fields to vectors
 */
class FieldFunction {

   public:

    virtual ~FieldFunction() = default;

    /**
     * Projects the field onto this supplied globVec according to the dm in the global Vedc
     * @param globVec
     * @param time
     */
    virtual void ProjectField(Vec globVec, PetscReal time) = 0;
};
}  // namespace ablate::mathFunctions
#endif  // ABLATELIBRARY_FIELDMATHFUNCTION_HPP
