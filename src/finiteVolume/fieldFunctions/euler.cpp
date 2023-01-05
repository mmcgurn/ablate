#include "euler.hpp"

ablate::finiteVolume::fieldFunctions::Euler::Euler(std::shared_ptr<ablate::finiteVolume::fieldFunctions::CompressibleFlowState> flowStateIn, std::shared_ptr<ablate::domain::Region> region)
    : ablate::domain::FieldMathFunction("euler", flowStateIn->GetFieldFunction("euler"), region) {}

#include "registrar.hpp"
REGISTER(ablate::domain::FieldMathFunction, ablate::finiteVolume::fieldFunctions::Euler, "initializes the euler conserved field variables based upon a CompressibleFlowState",
         ARG(ablate::finiteVolume::fieldFunctions::CompressibleFlowState, "state", "The CompressibleFlowState used to initialize"),
         OPT(ablate::domain::Region, "region", "A subset of the domain to apply the field function"));
