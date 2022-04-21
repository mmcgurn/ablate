#include "particles.hpp"
#include <petscviewerhdf5.h>
#include <vector>
#include "utilities/mpiError.hpp"
#include "utilities/petscError.hpp"
#include "utilities/petscOptions.hpp"

ablate::particles::Particles::Particles(std::string solverId, std::shared_ptr<domain::Region> region, std::shared_ptr<parameters::Parameters> options, int ndims, std::vector<ParticleField> fields,
                                        std::shared_ptr<particles::initializers::Initializer> initializer, std::vector<std::shared_ptr<mathFunctions::FieldFunction>> fieldInitialization,
                                        std::shared_ptr<mathFunctions::MathFunction> exactSolution)
    : Solver(solverId, region, options),
      swarmDm(nullptr),
      ndims(ndims),
      timeInitial(0.0),
      timeFinal(0.0),
      exactSolution(exactSolution),
      dmChanged(false),
      initializer(initializer),
      fieldInitialization(fieldInitialization) {
    // create and associate the dm
    DMCreate(PETSC_COMM_WORLD, &swarmDm) >> checkError;
    DMSetType(swarmDm, DMSWARM) >> checkError;
    DMSetDimension(swarmDm, ndims) >> checkError;
    DMSwarmSetType(swarmDm, DMSWARM_PIC) >> checkError;

    // Record the default fields
    std::vector<std::string> coordComponents;
    switch (ndims) {
        case 1:
            coordComponents = {"X"};
            break;
        case 2:
            coordComponents = {"X", "Y"};
            break;
        case 3:
            coordComponents = {"X", "Y", "Z"};
            break;
        default:
            throw std::invalid_argument("Particles ndims must be 1, 2, or 3. " + std::to_string(ndims) + " is not valid.");
    }
    auto coordField = ParticleField{.name = DMSwarmPICField_coor, .components = coordComponents, .type = domain::FieldLocation::SOL, .dataType = PETSC_REAL};
    particleFieldDescriptors.push_back(coordField);
    particleSolutionFieldDescriptors.push_back(coordField);
    particleFieldDescriptors.emplace_back(ParticleField{.name = DMSwarmField_pid, .type = domain::FieldLocation::AUX, .dataType = PETSC_INT64});

    // if the exact solution was provided, register the initial particle location in the field
    if (exactSolution) {
        fields.push_back(ParticleField{.name = ParticleInitialLocation, .components = coordComponents, .type = domain::FieldLocation::AUX, .dataType = PETSC_REAL});
    }

    // register each field
    for (auto &field : fields) {
        RegisterParticleField(field);
    }
}

void ablate::particles::Particles::Setup() {}

void ablate::particles::Particles::Initialize() {
    // if more than one solution field is provided, create a new field to hold them packed together
    if (particleSolutionFieldDescriptors.size() > 1) {
        auto packedSolutionComponentSize = 0;

        for (const auto &solution : particleSolutionFieldDescriptors) {
            packedSolutionComponentSize += solution.components.size();
        }

        // Compute the size of the exact solution (each component added up)
        RegisterParticleField(
            ParticleField{.name = PackedSolution, .components = std::vector<std::string>(packedSolutionComponentSize, "_"), .type = domain::FieldLocation::AUX, .dataType = PETSC_REAL});
    }

    // before setting up the flow finalize the fields
    DMSwarmFinalizeFieldRegister(swarmDm) >> checkError;

    // associate the swarm with the cell dm
    DMSwarmSetCellDM(swarmDm, subDomain->GetDM()) >> checkError;

    // name the particle domain
    PetscObjectSetOptions((PetscObject)swarmDm, petscOptions) >> checkError;
    PetscObjectSetName((PetscObject)swarmDm, GetId().c_str()) >> checkError;
    DMSetFromOptions(swarmDm) >> checkError;

    // initialize the particles
    initializer->Initialize(*subDomain, swarmDm);

    // Setup particle position integrator
    TSCreate(subDomain->GetComm(), &particleTs) >> checkError;
    PetscObjectSetOptions((PetscObject)particleTs, petscOptions) >> checkError;
    TSSetApplicationContext(particleTs, this) >> checkError;

    // Link thw dm
    TSSetDM(particleTs, swarmDm);
    TSSetProblemType(particleTs, TS_NONLINEAR) >> checkError;
    TSSetExactFinalTime(particleTs, TS_EXACTFINALTIME_MATCHSTEP) >> checkError;
    TSSetMaxSteps(particleTs, 100000000) >> checkError;  // set the max ts to a very large number. This can be over written using ts_max_steps options

    // finish ts setup
    TSSetFromOptions(particleTs) >> checkError;

    // set the functions to compute error is provided
    if (exactSolution) {
        StoreInitialParticleLocations();
        TSSetComputeExactError(particleTs, ComputeParticleError) >> checkError;
    }

    // project the initialization field onto each local particle
    for (auto &field : fieldInitialization) {
        this->ProjectFunction(field->GetName(), field->GetSolutionField());
    }
}

ablate::particles::Particles::~Particles() {
    if (swarmDm) {
        DMDestroy(&swarmDm) >> checkError;
    }
    if (particleTs) {
        TSDestroy(&particleTs) >> checkError;
    }
    if (petscOptions) {
        ablate::utilities::PetscOptionsDestroyAndCheck(GetId(), &petscOptions);
    }
}

void ablate::particles::Particles::RegisterParticleField(const ParticleField &field) {
    // add the value to the field
    DMSwarmRegisterPetscDatatypeField(swarmDm, field.name.c_str(), field.components.size(), field.dataType) >> checkError;

    // store the field
    particleFieldDescriptors.push_back(field);
    if (field.type == domain::FieldLocation::SOL) {
        particleSolutionFieldDescriptors.push_back(field);
    }
}

void ablate::particles::Particles::StoreInitialParticleLocations() {
    // copy over the initial location
    PetscReal *coord;
    PetscReal *initialLocation;
    PetscInt numberParticles;
    DMSwarmGetLocalSize(swarmDm, &numberParticles) >> checkError;
    DMSwarmGetField(swarmDm, DMSwarmPICField_coor, NULL, NULL, (void **)&coord) >> checkError;
    DMSwarmGetField(swarmDm, ParticleInitialLocation, NULL, NULL, (void **)&initialLocation) >> checkError;

    // copy the raw data
    for (int i = 0; i < numberParticles * ndims; ++i) {
        initialLocation[i] = coord[i];
    }
    DMSwarmRestoreField(swarmDm, DMSwarmPICField_coor, NULL, NULL, (void **)&coord) >> checkError;
    DMSwarmRestoreField(swarmDm, ParticleInitialLocation, NULL, NULL, (void **)&initialLocation) >> checkError;
}

PetscErrorCode ablate::particles::Particles::ComputeParticleExactSolution(TS particleTS, Vec exactSolution) {
    PetscFunctionBeginUser;

    // get a pointer to this particle class
    ablate::particles::Particles *particles;
    TSGetApplicationContext(particleTS, (void **)&particles) >> checkError;

    // get the abs time for the particle evaluation, this is the ts relative time plus the time at the start of the particle ts solve
    PetscReal time;
    TSGetTime(particleTS, &time) >> checkError;
    time += particles->timeInitial;

    // Create a vector of the current solution
    PetscScalar *exactSolutionArray;
    VecGetArrayWrite(exactSolution, &exactSolutionArray) >> checkError;

    // exact the exact solution from the initial location
    PetscInt np;
    DMSwarmGetLocalSize(particles->swarmDm, &np) >> checkError;
    const PetscInt dim = particles->ndims;

    // Calculate the size of solution field
    PetscInt solutionFieldSize = 0;
    for (const auto &field : particles->particleSolutionFieldDescriptors) {
        solutionFieldSize += field.components.size();
    }

    // get the initial location array
    const PetscScalar *initialParticleLocationArray;
    DMSwarmGetField(particles->swarmDm, ParticleInitialLocation, NULL, NULL, (void **)&initialParticleLocationArray) >> checkError;

    // extract the petsc function for fast update
    void *functionContext = particles->exactSolution->GetContext();
    ablate::mathFunctions::PetscFunction functionPointer = particles->exactSolution->GetPetscFunction();

    // for each local particle, get the exact location and other variables
    for (PetscInt p = 0; p < np; ++p) {
        // compute the array offset
        const PetscInt initialPositionOffset = p * dim;
        const PetscInt fieldOffset = p * solutionFieldSize;

        // Call the update function
        functionPointer(dim, time, initialParticleLocationArray + initialPositionOffset, solutionFieldSize, exactSolutionArray + fieldOffset, functionContext) >> checkError;
    }
    VecRestoreArrayWrite(exactSolution, &exactSolutionArray) >> checkError;

    // cleanup
    DMSwarmRestoreField(particles->swarmDm, ParticleInitialLocation, NULL, NULL, (void **)&initialParticleLocationArray) >> checkError;

    PetscFunctionReturn(0);
}

PetscErrorCode ablate::particles::Particles::ComputeParticleError(TS particleTS, Vec u, Vec errorVec) {
    PetscFunctionBeginUser;

    // get a pointer to this particle class
    ablate::particles::Particles *particles;
    TSGetApplicationContext(particleTS, (void **)&particles) >> checkError;

    // get the abs time for the particle evaluation, this is the ts relative time plus the time at the start of the particle ts solve
    PetscReal time;
    TSGetTime(particleTS, &time) >> checkError;
    time += particles->timeInitial;

    // Create a vector of the current solution
    Vec exactSolutionVec;
    VecDuplicate(u, &exactSolutionVec) >> checkError;
    PetscScalar *exactSolutionArray;
    VecGetArrayWrite(exactSolutionVec, &exactSolutionArray) >> checkError;

    // Also store the exact location separately
    DMSwarmVectorDefineField(particles->swarmDm, ParticleInitialLocation) >> checkError;
    Vec exactLocationVec;
    DMGetGlobalVector(particles->swarmDm, &exactLocationVec);
    PetscScalar *exactLocationArray;
    VecGetArrayWrite(exactLocationVec, &exactLocationArray) >> checkError;

    // exact the exact solution from the initial location
    PetscInt np;
    DMSwarmGetLocalSize(particles->swarmDm, &np) >> checkError;
    const PetscInt dim = particles->ndims;

    // Calculate the size of solution field
    PetscInt solutionFieldSize = 0;
    for (const auto &field : particles->particleSolutionFieldDescriptors) {
        solutionFieldSize += field.components.size();
    }

    // get the initial location array
    const PetscScalar *initialParticleLocationArray;
    DMSwarmGetField(particles->swarmDm, ParticleInitialLocation, NULL, NULL, (void **)&initialParticleLocationArray) >> checkError;

    // extract the petsc function for fast update
    void *functionContext = particles->exactSolution->GetContext();
    ablate::mathFunctions::PetscFunction functionPointer = particles->exactSolution->GetPetscFunction();

    // for each local particle, get the exact location and other variables
    for (PetscInt p = 0; p < np; ++p) {
        // compute the array offset
        const PetscInt initialPositionOffset = p * dim;
        const PetscInt fieldOffset = p * solutionFieldSize;

        // Call the update function
        functionPointer(dim, time, initialParticleLocationArray + initialPositionOffset, solutionFieldSize, exactSolutionArray + fieldOffset, functionContext) >> checkError;

        // copy over the first dim to the exact solution array
        for (PetscInt d = 0; d < dim; ++d) {
            exactLocationArray[initialPositionOffset + d] = exactSolutionArray[fieldOffset + d];
        }
    }
    VecRestoreArrayWrite(exactSolutionVec, &exactSolutionArray) >> checkError;
    VecRestoreArrayWrite(exactLocationVec, &exactLocationArray) >> checkError;

    // Get all points still in this mesh
    DM flowDM = particles->subDomain->GetDM();
    PetscSF cellSF = NULL;
    DMLocatePoints(flowDM, exactLocationVec, DM_POINTLOCATION_NONE, &cellSF) >> checkError;
    const PetscSFNode *cells;
    PetscSFGetGraph(cellSF, NULL, NULL, NULL, &cells) >> checkError;

    // compute the difference between exact and u
    VecWAXPY(errorVec, -1, exactSolutionVec, u);

    // zero out the error if any particle moves outside of the domain
    for (PetscInt p = 0; p < np; ++p) {
        if (cells[p].index == DMLOCATEPOINT_POINT_NOT_FOUND) {
            for (PetscInt c = 0; c < solutionFieldSize; ++c) {
                VecSetValue(errorVec, p * solutionFieldSize + c, 0.0, INSERT_VALUES) >> checkError;
            }
        }
    }
    VecAssemblyBegin(errorVec) >> checkError;
    VecAssemblyEnd(errorVec) >> checkError;

    // restore all of the vecs/fields
    PetscSFDestroy(&cellSF) >> checkError;

    // cleanup
    DMSwarmRestoreField(particles->swarmDm, ParticleInitialLocation, NULL, NULL, (void **)&initialParticleLocationArray) >> checkError;
    VecDestroy(&exactSolutionVec) >> checkError;
    DMRestoreGlobalVector(particles->swarmDm, &exactLocationVec) >> checkError;

    PetscFunctionReturn(0);
}

void ablate::particles::Particles::SwarmMigrate() {
    // current number of local/global particles
    PetscInt numberLocal;
    PetscInt numberGlobal;

    // Get the current size
    DMSwarmGetLocalSize(swarmDm, &numberLocal) >> checkError;
    DMSwarmGetSize(swarmDm, &numberGlobal) >> checkError;

    // Migrate any particles that have moved
    DMSwarmMigrate(swarmDm, PETSC_TRUE) >> checkError;

    // get the new sizes
    PetscInt newNumberLocal;
    PetscInt newNumberGlobal;

    // Get the updated size
    DMSwarmGetLocalSize(swarmDm, &newNumberLocal) >> checkError;
    DMSwarmGetSize(swarmDm, &newNumberGlobal) >> checkError;

    // Check to see if any of the ranks changed size after migration
    PetscMPIInt dmChangedLocal = newNumberGlobal != numberGlobal || newNumberLocal != numberLocal;
    MPI_Comm comm;
    PetscObjectGetComm((PetscObject)particleTs, &comm) >> checkError;
    PetscMPIInt dmChangedAll = PETSC_FALSE;
    MPI_Allreduce(&dmChangedLocal, &dmChangedAll, 1, MPIU_INT, MPIU_MAX, comm) >> checkMpiError;
    dmChanged = dmChangedAll == PETSC_TRUE;
}

/**
 * Support function to project the math function onto a particle field
 * @param field
 * @param mathFunction
 */
void ablate::particles::Particles::ProjectFunction(const std::string &field, ablate::mathFunctions::MathFunction &mathFunction) {
    // Get the local number of particles
    PetscInt np;
    DMSwarmGetLocalSize(swarmDm, &np) >> checkError;

    // Get the raw access to position and update field
    PetscInt dim;
    PetscReal *positionData;
    DMSwarmGetField(swarmDm, DMSwarmPICField_coor, &dim, NULL, (void **)&positionData) >> checkError;

    PetscInt fieldComponents;
    PetscDataType fieldType;
    PetscReal *fieldData;
    DMSwarmGetField(swarmDm, field.c_str(), &fieldComponents, &fieldType, (void **)&fieldData) >> checkError;

    if (fieldType != PETSC_REAL) {
        throw std::invalid_argument("ProjectFunction only supports PETSC_REAL");
    }

    // extract the petsc function for fast update
    void *functionContext = mathFunction.GetContext();
    ablate::mathFunctions::PetscFunction functionPointer = mathFunction.GetPetscFunction();

    // Iterate over each local particle
    for (PetscInt p = 0; p < np; ++p) {
        // compute the position offset
        const PetscInt positionOffset = p * dim;

        // Compute the field offset
        const PetscInt fieldOffset = p * fieldComponents;

        // Call the update function
        functionPointer(dim, 0.0, positionData + positionOffset, fieldComponents, fieldData + fieldOffset, functionContext) >> checkError;
    }
    DMSwarmRestoreField(swarmDm, DMSwarmPICField_coor, NULL, NULL, (void **)&positionData);
    DMSwarmRestoreField(swarmDm, field.c_str(), NULL, NULL, (void **)&fieldData);
}

Vec ablate::particles::Particles::GetPackedSolutionVector() {
    const PetscInt nf = particleSolutionFieldDescriptors.size();

    // If there is more than one field, pack up the data
    if (nf > 1) {
        // Get the local number of particles
        PetscInt np;
        DMSwarmGetLocalSize(swarmDm, &np) >> checkError;

        // Get a vector of pointers to the solution fields
        std::vector<PetscReal *> fieldDatas(nf);
        std::vector<PetscInt> fieldSizes(nf);

        // get raw access to each array
        for (auto f = 0; f < nf; f++) {
            DMSwarmGetField(swarmDm, particleSolutionFieldDescriptors[f].name.c_str(), &fieldSizes[f], NULL, (void **)&fieldDatas[f]) >> checkError;
        }

        // and raw access to the combined array
        PetscInt solutionComponents;
        PetscReal *solutionFieldData;
        DMSwarmGetField(swarmDm, PackedSolution, &solutionComponents, NULL, (void **)&solutionFieldData) >> checkError;

        for (PetscInt p = 0; p < np; ++p) {
            PetscInt offset = 0;
            // March over each field
            for (PetscInt f = 0; f < nf; ++f) {
                for (PetscInt c = 0; c < fieldSizes[f]; c++) {
                    solutionFieldData[p * solutionComponents + offset++] = fieldDatas[f][p * fieldSizes[f] + c];
                }
            }
        }

        // return raw access
        for (auto f = 0; f < nf; f++) {
            DMSwarmRestoreField(swarmDm, particleSolutionFieldDescriptors[f].name.c_str(), NULL, NULL, (void **)&fieldDatas[f]) >> checkError;
        }
        DMSwarmRestoreField(swarmDm, PackedSolution, NULL, NULL, (void **)&solutionFieldData) >> checkError;
    }

    // get the updated values as a vec
    Vec packedVector;
    DMSwarmCreateGlobalVectorFromField(swarmDm, GetSolutionVectorName(), &packedVector) >> checkError;
    return packedVector;
}

void ablate::particles::Particles::RestorePackedSolutionVector(Vec solutionVector) {
    const PetscInt nf = particleSolutionFieldDescriptors.size();

    DMSwarmDestroyGlobalVectorFromField(swarmDm, GetSolutionVectorName(), &solutionVector) >> checkError;
    // If there is more than one field, unpack the data
    if (nf > 1) {
        // Get the local number of particle
        PetscInt np;
        DMSwarmGetLocalSize(swarmDm, &np) >> checkError;

        // Get a vector of pointers to the solution fields
        std::vector<PetscReal *> fieldDatas(nf);
        std::vector<PetscInt> fieldSizes(nf);

        // get raw access to each array
        for (auto f = 0; f < nf; f++) {
            DMSwarmGetField(swarmDm, particleSolutionFieldDescriptors[f].name.c_str(), &fieldSizes[f], NULL, (void **)&fieldDatas[f]) >> checkError;
        }

        // and raw access to the combined array
        PetscInt solutionComponents;
        PetscReal *solutionFieldData;
        DMSwarmGetField(swarmDm, PackedSolution, &solutionComponents, NULL, (void **)&solutionFieldData) >> checkError;

        for (PetscInt p = 0; p < np; ++p) {
            PetscInt offset = 0;
            // March over each field
            for (PetscInt f = 0; f < nf; ++f) {
                for (PetscInt c = 0; c < fieldSizes[f]; c++) {
                    fieldDatas[f][p * fieldSizes[f] + c] = solutionFieldData[p * solutionComponents + offset++];
                }
            }
        }

        // return raw access
        for (auto f = 0; f < nf; f++) {
            DMSwarmRestoreField(swarmDm, particleSolutionFieldDescriptors[f].name.c_str(), NULL, NULL, (void **)&fieldDatas[f]) >> checkError;
        }
        DMSwarmRestoreField(swarmDm, PackedSolution, NULL, NULL, (void **)&solutionFieldData) >> checkError;
    }
}

void ablate::particles::Particles::AdvectParticles(TS flowTS) {
    PetscReal time;

    // if the dm has changed size (new particles, particles moved between ranks, particles deleted) reset the ts
    if (dmChanged) {
        TSReset(particleTs) >> checkError;
        dmChanged = PETSC_FALSE;
    }

    // Get the position, velocity and Kinematics vector
    Vec solutionVector = GetPackedSolutionVector();
    // get the particle time step
    PetscReal dtInitial;
    TSGetTimeStep(particleTs, &dtInitial) >> checkError;

    // Set the max end time based upon the flow end time
    TSGetTime(flowTS, &time) >> checkError;
    TSSetMaxTime(particleTs, time) >> checkError;
    timeFinal = time;

    // take the needed timesteps to get to the flow time
    TSSolve(particleTs, solutionVector) >> checkError;
    timeInitial = timeFinal;

    // get the updated time step, and reset if it has gone down
    PetscReal dtUpdated;
    TSGetTimeStep(particleTs, &dtUpdated) >> checkError;
    if (dtUpdated < dtInitial) {
        TSSetTimeStep(particleTs, dtInitial) >> checkError;
    }

    RestorePackedSolutionVector(solutionVector);

    // Migrate any particles that have moved
    SwarmMigrate();
}

static PetscErrorCode DMSequenceViewTimeHDF5(DM dm, PetscViewer viewer) {
    Vec stamp;
    PetscMPIInt rank;
    PetscErrorCode ierr;

    PetscFunctionBegin;

    // get the seqnum and value from the dm
    PetscInt seqnum;
    PetscReal value;
    ierr = DMGetOutputSequenceNumber(dm, &seqnum, &value);
    CHKERRMPI(ierr);

    if (seqnum < 0) {
        PetscFunctionReturn(0);
    }
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)viewer), &rank);
    CHKERRMPI(ierr);
    ierr = VecCreateMPI(PetscObjectComm((PetscObject)viewer), rank ? 0 : 1, 1, &stamp);
    CHKERRQ(ierr);
    ierr = VecSetBlockSize(stamp, 1);
    CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)stamp, "time");
    CHKERRQ(ierr);
    if (!rank) {
        ierr = VecSetValue(stamp, 0, value, INSERT_VALUES);
        CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(stamp);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(stamp);
    CHKERRQ(ierr);
    ierr = PetscViewerHDF5PushGroup(viewer, "/");
    CHKERRQ(ierr);
    ierr = PetscViewerHDF5SetTimestep(viewer, seqnum);
    CHKERRQ(ierr);
    ierr = VecView(stamp, viewer);
    CHKERRQ(ierr);
    ierr = PetscViewerHDF5PopGroup(viewer);
    CHKERRQ(ierr);
    ierr = VecDestroy(&stamp);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

void ablate::particles::Particles::Save(PetscViewer viewer, PetscInt steps, PetscReal time) {
    DMSetOutputSequenceNumber(GetParticleDM(), steps, time) >> checkError;
    Vec particleVector;

    // if this is an hdf5Viewer
    PetscBool ishdf5;
    PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERHDF5, &ishdf5) >> checkError;
    if (ishdf5) {
        PetscBool isInTimestepping;
        PetscViewerHDF5IsTimestepping(viewer, &isInTimestepping) >> checkError;
        if (!isInTimestepping) {
            PetscViewerHDF5PushTimestepping(viewer) >> checkError;
        }
    }

    for (auto const &field : particleFieldDescriptors) {
        if (field.dataType == PETSC_REAL) {
            DMSwarmCreateGlobalVectorFromField(GetParticleDM(), field.name.c_str(), &particleVector) >> checkError;
            PetscObjectSetName((PetscObject)particleVector, field.name.c_str()) >> checkError;
            VecView(particleVector, viewer) >> checkError;
            DMSwarmDestroyGlobalVectorFromField(GetParticleDM(), field.name.c_str(), &particleVector) >> checkError;
        }
    }

    // Get the particle info
    int rank;
    MPI_Comm_rank(PetscObjectComm((PetscObject)GetParticleDM()), &rank) >> checkMpiError;

    // get the local number of particles
    PetscInt globalSize;
    DMSwarmGetSize(GetParticleDM(), &globalSize) >> checkMpiError;

    // record the number of particles per rank
    Vec particleCountVec;
    VecCreateMPI(PetscObjectComm((PetscObject)GetParticleDM()), PETSC_DECIDE, 1, &particleCountVec) >> checkError;
    PetscObjectSetName((PetscObject)particleCountVec, "particleCount") >> checkError;
    VecSetValue(particleCountVec, 0, globalSize, INSERT_VALUES) >> checkError;
    VecAssemblyBegin(particleCountVec) >> checkError;
    VecAssemblyEnd(particleCountVec) >> checkError;
    VecView(particleCountVec, viewer);
    VecDestroy(&particleCountVec) >> checkError;

    if (ishdf5) {
        DMSequenceViewTimeHDF5(GetParticleDM(), viewer) >> checkError;
    }
}

void ablate::particles::Particles::Restore(PetscViewer viewer, PetscInt sequenceNumber, PetscReal time) {
    DMSetOutputSequenceNumber(GetParticleDM(), sequenceNumber, time) >> checkError;

    // Update the ts with the current values
    TSSetTime(particleTs, time) >> checkError;
    timeInitial = time;

    // There is not a hdf5 specific swarm vec load, so that needs to be in this code
    PetscBool ishdf5;
    PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERHDF5, &ishdf5) >> checkError;
    if (ishdf5) {
        PetscViewerHDF5PushTimestepping(viewer) >> checkError;
        PetscViewerHDF5SetTimestep(viewer, sequenceNumber) >> checkError;
    }

    // load in the global particle size
    Vec particleCountVec;
    VecCreateSeq(PETSC_COMM_SELF, 1, &particleCountVec) >> checkError;
    PetscObjectSetName((PetscObject)particleCountVec, "particleCount") >> checkError;
    VecLoad(particleCountVec, viewer) >> checkError;

    PetscReal globalSize;
    PetscInt index[1] = {0};
    VecGetValues(particleCountVec, 1, index, &globalSize) >> checkError;
    VecDestroy(&particleCountVec) >> checkError;

    // Get the particle mpi, info
    int rank, size;
    MPI_Comm_rank(PetscObjectComm((PetscObject)GetParticleDM()), &rank) >> checkMpiError;
    MPI_Comm_size(PetscObjectComm((PetscObject)GetParticleDM()), &size) >> checkMpiError;

    // distribute the number of particles across all ranks
    PetscInt localSize = ((PetscInt)globalSize) / size;

    // Use the first rank to hold any left over
    if (rank == 0) {
        localSize = globalSize - (localSize * (size - 1));
    }

    // Set the local swarm size
    DMSwarmSetLocalSizes(GetParticleDM(), localSize, 0) >> checkError;

    // Move in the hdf5 to the right group
    if (ishdf5) {
        PetscViewerHDF5PushGroup(viewer, "/particle_fields") >> checkError;
        PetscViewerHDF5SetTimestep(viewer, sequenceNumber) >> checkError;
    }

    for (auto const &field : particleFieldDescriptors) {
        if (field.dataType == PETSC_REAL) {
            Vec particleVector;
            Vec particleVectorLoad;
            DMSwarmCreateGlobalVectorFromField(swarmDm, field.name.c_str(), &particleVector) >> checkError;

            // A copy of this vector is needed, because vec load breaks the memory linkage between the swarm and vec
            VecDuplicate(particleVector, &particleVectorLoad) >> checkError;

            // Load the vector
            PetscObjectSetName((PetscObject)particleVectorLoad, field.name.c_str()) >> checkError;
            VecLoad(particleVectorLoad, viewer) >> checkError;

            // Copy the data over
            VecCopy(particleVectorLoad, particleVector) >> checkError;

            DMSwarmDestroyGlobalVectorFromField(swarmDm, field.name.c_str(), &particleVector) >> checkError;
            VecDestroy(&particleVectorLoad) >> checkError;
        }
    }

    if (ishdf5) {
        PetscViewerHDF5PopGroup(viewer) >> checkError;
        PetscViewerHDF5PopTimestepping(viewer) >> checkError;
    }

    // Migrate the particle to the correct rank for the dmPlex
    DMSwarmMigrate(swarmDm, PETSC_TRUE) >> checkError;
    dmChanged = true;
}
std::vector<std::string> ablate::particles::Particles::CreateDimensionVector(const std::string &prefix, int dim) {
    std::vector<std::string> vector;

    for (int i = 0; i < dim; i++) {
        vector.emplace_back(prefix + std::to_string(dim));
    }

    return vector;
}