#ifndef ABLATELIBRARY_ZERORK_SOURCECALCULATOR_HPP
#define ABLATELIBRARY_ZERORK_SOURCECALCULATOR_HPP

#include "eos/chemistryModel.hpp"
#include "zerork_cfd_plugin.h"
#include "zerork/mechanism.h"
#include "zerork/utilities.h"

#include "zerork_cfd_plugin.h"

namespace ablate::eos {
class zerorkEOS;
}

namespace ablate::eos::zerorkeos {

/**
 * public class to to compute the source for each specified node
 */
class SourceCalculator : public ChemistryModel::SourceCalculator, private utilities::Loggable<SourceCalculator> {
   public:
    /**
     * Allow the user to set the reactor type
     */
    enum class ReactorType { ConstantPressure, ConstantVolume };

    //! hold a struct that can be used for chemistry constraints
    struct ChemistryConstraints {
        double relTolerance = 1.0E-4;
        double absTolerance = 1.0E-8;

        // zerork output state
        bool verbose = false;

        // For large mechanisms(~Nspec>100) it is recommended to use a sparse math jacobian
        bool sparseJacobian = false;

        // store the reactor type in the chemistry constrains
        ReactorType reactorType = ReactorType::ConstantVolume;

        // store an optional threshold temperature.  Only compute the reactions if the temperature is above thresholdTemperature
        double thresholdTemperature = 0.0;

        void Set(const std::shared_ptr<ablate::parameters::Parameters>&);
    };
    /**
     * create a batch source for this size specified in cellRange
     * @param zerorkEos
     * @param constraints
     * @param cellRange
     */
    SourceCalculator(const std::vector<domain::Field>& fields, std::shared_ptr<zerorkEOS> zerorkEos, ablate::eos::zerorkeos::SourceCalculator::ChemistryConstraints constraints, const ablate::domain::Range& cellRange);

    /**
     * The compute source can be used as a prestep allowing the add source to be used at each stage without reevaluating
     */
    void ComputeSource(const ablate::domain::Range& cellRange, PetscReal time, PetscReal dt, Vec globalSolution) override;

    /**
     * Adds the source that was computed in the ComputeSource to the supplied vector
     */
    void AddSource(const ablate::domain::Range& cellRange, Vec localXVec, Vec localFVec) override;

   private:
    std::vector<double> sourceZeroRKAtI;
    zerork_handle zrm_handle;
    //! copy of constraints
    ablate::eos::zerorkeos::SourceCalculator::ChemistryConstraints chemistryConstraints;
    /**
     * Hold access to the zerorkEOS eos needed to create eos
     */
    std::shared_ptr<eos::zerorkEOS> eos;

    const size_t numberSpecies;

    //! the id for the required euler field
    PetscInt eulerId;

    //! the id for the required densityYi field
    PetscInt densityYiId;



};

/**
 * Support function
 * @param os
 * @param v
 * @return
 */
std::ostream& operator<<(std::ostream& os, const SourceCalculator::ReactorType& v);

/**
 * Support function
 * @param os
 * @param v
 * @return
 */
std::istream& operator>>(std::istream& is, SourceCalculator::ReactorType& v);

}  // namespace

#endif  // ABLATELIBRARY_BATCHSOURCE_HPP
