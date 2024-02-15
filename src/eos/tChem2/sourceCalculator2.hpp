#ifndef ABLATELIBRARY_TCHEM2_SOURCECALCULATOR_HPP
#define ABLATELIBRARY_TCHEM2_SOURCECALCULATOR_HPP

#include <TChem_KineticModelGasConstData.hpp>
#include "eos/chemistryModel.hpp"
#include "eos/tChem/sourceCalculator.hpp"
#include "zerork_cfd_plugin.h"
#include "zerork/mechanism.h"

namespace tChemLib = TChem;

namespace ablate::eos {
class TChem2;
}

namespace ablate::eos::tChem2 {

/**
 * public class to to compute the source for each specified node
 */
class SourceCalculator2 : public ChemistryModel::SourceCalculator, private utilities::Loggable<SourceCalculator2> {
   public:
    /**
     * Allow the user of TChem to set the reactor type
     */
    enum class ReactorType { ConstantPressure, ConstantVolume };

    //! hold a struct that can be used for chemistry constraints

    /**
     * create a batch source for this size specified in cellRange
     * @param tChemEos
     * @param constraints
     * @param cellRange
     */
    SourceCalculator2(const std::vector<domain::Field>& fields, std::shared_ptr<TChem2> tChemEos, ablate::eos::tChem::SourceCalculator::ChemistryConstraints constraints, const ablate::domain::Range& cellRange);
//    SourceCalculator(const std::vector<domain::Field>& fields, std::shared_ptr<TChem> tChemEos, ChemistryConstraints constraints, const ablate::domain::Range& cellRange);

    /**
     * The compute source can be used as a prestep allowing the add source to be used at each stage without reevaluating
     */
    void ComputeSource(const ablate::domain::Range& cellRange, PetscReal time, PetscReal dt, Vec globalSolution) override;

    /**
     * The compute source can be used as a prestep allowing the add source to be used at each stage without reevaluating
     */
    //    static void ComputeSource(SourceCalculator& sourceCalculator, const ablate::domain::Range& cellRange, PetscReal time, PetscReal dt, Vec globalSolution);

    /**
     * Adds the source that was computed in the ComputeSource to the supplied vector
     */
    void AddSource(const ablate::domain::Range& cellRange, Vec localXVec, Vec localFVec) override;

   private:

    zerork_handle zrm_handle;
    //! copy of constraints
//    ChemistryConstraints chemistryConstraints;
    ablate::eos::tChem::SourceCalculator::ChemistryConstraints chemistryConstraints;
    /**
     * Hold access to the tchem eos needed to create eos
     */
    std::shared_ptr<eos::TChem2> eos;

    const size_t numberSpecies;

    //! the id for the required euler field
    PetscInt eulerId;

    //! the id for the required densityYi field
    PetscInt densityYiId;

    // tchem memory storage on host/device.  These will be sized for the number of active nodes in the domain
    real_type_2d_view stateDevice;
    real_type_2d_view_host stateHost;

    // store the end state for the device/host
    real_type_2d_view endStateDevice;

    // the time advance information
    time_advance_type_1d_view timeAdvanceDevice;
    time_advance_type timeAdvanceDefault{};

    // store host/device memory for computing state
    real_type_1d_view internalEnergyRefDevice;
    real_type_1d_view_host internalEnergyRefHost;
    real_type_2d_view perSpeciesScratchDevice;

    // store the source terms (density* energy + density*species)
    real_type_2d_view_host sourceTermsHost;
    real_type_2d_view sourceTermsDevice;

    // tolerance constraints
    real_type_2d_view tolTimeDevice;
    real_type_1d_view tolNewtonDevice;
    real_type_2d_view facDevice;

    // store the time and delta for the ode solver
    real_type_1d_view timeViewDevice;
    real_type_1d_view dtViewDevice;

    // Hard code some values needed for the constant volume reactor
    static inline constexpr bool solveTla = false;   // do not calculate tangent linear approximation (TLA) for the const volume reactions
    static inline constexpr real_type thetaTla = 0;  // this is not used when solveTla is false

    // store device specific kineticModelGasConstants
    tChemLib::KineticModelConstData<typename Tines::UseThisDevice<exec_space>::type> kineticModelGasConstDataDevice;
    kmd_type_1d_view_host kineticModelDataClone;
    Kokkos::View<KineticModelGasConstData<typename Tines::UseThisDevice<exec_space>::type>*, typename Tines::UseThisDevice<exec_space>::type> kineticModelGasConstDataDevices;
};

/**
 * Support function for the TChemBase::ReactorType Enum
 * @param os
 * @param v
 * @return
 */
std::ostream& operator<<(std::ostream& os, const SourceCalculator2::ReactorType& v);

/**
 * Support function for the TChemBase::ReactorType Enum
 * @param os
 * @param v
 * @return
 */
std::istream& operator>>(std::istream& is, SourceCalculator2::ReactorType& v);

}  // namespace ablate::eos::tChem

#endif  // ABLATELIBRARY_BATCHSOURCE_HPP