#include "densityProgressVariables.hpp"
#include "finiteVolume/compressibleFlowFields.hpp"

ablate::finiteVolume::fieldFunctions::DensityProgressVariables::DensityProgressVariables(std::shared_ptr<ablate::finiteVolume::fieldFunctions::CompressibleFlowState> flowStateIn,
                                                                                         std::shared_ptr<ablate::domain::Region> region)
    : ablate::domain::FieldMathFunction(CompressibleFlowFields::DENSITY_PROGRESS_FIELD, flowStateIn->GetFieldFunction(CompressibleFlowFields::DENSITY_PROGRESS_FIELD), region) {}

#include "registrar.hpp"
REGISTER(ablate::domain::FieldMathFunction, ablate::finiteVolume::fieldFunctions::DensityProgressVariables,
         "initializes the density Progress Variable conserved field variables based upon a CompressibleFlowState",
         ARG(ablate::finiteVolume::fieldFunctions::CompressibleFlowState, "state", "The CompressibleFlowState used to initialize"),
         OPT(ablate::domain::Region, "region", "A subset of the domain to apply the field function"));