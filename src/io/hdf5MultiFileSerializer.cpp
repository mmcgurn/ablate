#include "hdf5MultiFileSerializer.hpp"
#include <petscviewerhdf5.h>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <utility>
#include "environment/runEnvironment.hpp"
#include "generators.hpp"
#include "utilities/petscOptions.hpp"

ablate::io::Hdf5MultiFileSerializer::Hdf5MultiFileSerializer(std::shared_ptr<ablate::io::interval::Interval> interval, std::shared_ptr<parameters::Parameters> options)
    : interval(std::move(interval)) {
    // Load the metadata from the file is available, otherwise set to 0
    auto restartFilePath = environment::RunEnvironment::Get().GetOutputDirectory() / "restart.rst";

    if (std::filesystem::exists(restartFilePath)) {
        resumed = true;
        auto yaml = YAML::LoadFile(restartFilePath);
        time = yaml["time"].as<PetscReal>();
        dt = yaml["dt"].as<PetscReal>();
        timeStep = yaml["timeStep"].as<PetscInt>();
        sequenceNumber = yaml["sequenceNumber"].as<PetscInt>();
    } else {
        resumed = false;
        time = NAN;
        dt = NAN;
        timeStep = -1;
        sequenceNumber = -1;
    }

    // setup petsc options if provided
    if (options) {
        PetscOptionsCreate(&petscOptions) >> checkError;
        options->Fill(petscOptions);
    }
}

ablate::io::Hdf5MultiFileSerializer::~Hdf5MultiFileSerializer() {
    // save each serializer
    for (std::string id : postProcessesIds) {
        std::vector<std::filesystem::path> inputFilePaths;

        auto directoryPath = GetOutputDirectoryPath(id);
        for (const auto& file : std::filesystem::directory_iterator(directoryPath)) {
            if (file.path().extension() == ".hdf5") {
                inputFilePaths.push_back(file.path());
            }
        }

        // sort the paths
        std::sort(inputFilePaths.begin(), inputFilePaths.end());

        // run the convert function
        std::filesystem::path outputFile = directoryPath / (id + ".xmf");
        xdmfGenerator::Generate(inputFilePaths, outputFile);
    }

    if (petscOptions) {
        ablate::utilities::PetscOptionsDestroyAndCheck("ablate::io::Hdf5MultiFileSerializer::Hdf5MultiFileSerializer", &petscOptions);
    }
}

void ablate::io::Hdf5MultiFileSerializer::Register(std::weak_ptr<Serializable> serializable) {
    serializables.push_back(serializable);

    // Mark this to clean up
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank) >> checkError;

    if (auto serializableObject = serializable.lock()) {
        // resume if needed
        if (resumed) {
            auto filePath = GetOutputFilePath(serializableObject->GetId());

            PetscViewer petscViewer = nullptr;
            StartEvent("PetscViewerHDF5Open");
            PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.string().c_str(), FILE_MODE_UPDATE, &petscViewer) >> checkError;
            EndEvent();

            // set the petsc options if provided
            PetscObjectSetOptions((PetscObject)petscViewer, petscOptions) >> checkError;
            PetscViewerSetFromOptions(petscViewer) >> checkError;
            PetscViewerViewFromOptions(petscViewer, nullptr, "-hdf5ViewerView") >> checkError;

            // Restore the simulation
            StartEvent("Restore");
            // NOTE: as far as the output file the sequence number is always zero because it is a new file
            serializableObject->Restore(petscViewer, 0, time);
            EndEvent();

            StartEvent("PetscViewerHDF5Destroy");
            PetscViewerDestroy(&petscViewer) >> checkError;
            EndEvent();
        } else {
            // Create an output directory
            auto outputDirectory = GetOutputDirectoryPath(serializableObject->GetId());
            if (rank == 0) {
                std::filesystem::create_directory(outputDirectory);
            }
            MPI_Barrier(PETSC_COMM_WORLD);
        }

        if (rank == 0) {
            postProcessesIds.push_back(serializableObject->GetId());
        }
    }
}

PetscErrorCode ablate::io::Hdf5MultiFileSerializer::Hdf5MultiFileSerializerSaveStateFunction(TS ts, PetscInt steps, PetscReal time, Vec, void* ctx) {
    PetscFunctionBeginUser;
    auto hdf5Serializer = (Hdf5MultiFileSerializer*)ctx;

    // Make sure that the same timeStep is not output more than once (this can be from a restart)
    if (steps <= hdf5Serializer->timeStep) {
        PetscFunctionReturn(0);
    }

    if (hdf5Serializer->interval->Check(PetscObjectComm((PetscObject)ts), steps, time)) {
        // Update all metadata
        hdf5Serializer->time = time;
        hdf5Serializer->timeStep = steps;
        hdf5Serializer->sequenceNumber++;
        TSGetTimeStep(ts, &(hdf5Serializer->dt)) >> checkError;

        // Save this to a file
        hdf5Serializer->SaveMetadata(ts);

        try {
            // save each serializer
            for (auto& serializablePtr : hdf5Serializer->serializables) {
                if (auto serializableObject = serializablePtr.lock()) {
                    // Create an output path
                    auto filePath = hdf5Serializer->GetOutputFilePath(serializableObject->GetId());

                    PetscViewer petscViewer = nullptr;
                    hdf5Serializer->StartEvent("PetscViewerHDF5Open");
                    PetscViewerHDF5Open(PETSC_COMM_WORLD, filePath.string().c_str(), FILE_MODE_WRITE, &petscViewer) >> checkError;
                    hdf5Serializer->EndEvent();

                    // set the petsc options if provided
                    PetscObjectSetOptions((PetscObject)petscViewer, hdf5Serializer->petscOptions) >> checkError;
                    PetscViewerSetFromOptions(petscViewer) >> checkError;
                    PetscViewerViewFromOptions(petscViewer, nullptr, "-hdf5ViewerView") >> checkError;

                    hdf5Serializer->StartEvent("Save");
                    // NOTE: as far as the output file the sequence number is always zero because it is a new file
                    serializableObject->Save(petscViewer, 0, time);
                    hdf5Serializer->EndEvent();

                    hdf5Serializer->StartEvent("PetscViewerHDF5Destroy");
                    PetscViewerDestroy(&petscViewer) >> checkError;
                    hdf5Serializer->EndEvent();
                }
            }
        } catch (std::exception& exception) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB, "%s", exception.what());
        }
    }
    PetscFunctionReturn(0);
}

void ablate::io::Hdf5MultiFileSerializer::SaveMetadata(TS ts) const {
    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "time";
    out << YAML::Value << time;
    out << YAML::Key << "dt";
    out << YAML::Value << dt;
    out << YAML::Key << "timeStep";
    out << YAML::Value << timeStep;
    out << YAML::Key << "sequenceNumber";
    out << YAML::Value << sequenceNumber;
    out << YAML::EndMap;

    int rank;
    MPI_Comm_rank(PetscObjectComm((PetscObject)ts), &rank) >> checkError;
    if (rank == 0) {
        auto restartFilePath = environment::RunEnvironment::Get().GetOutputDirectory() / "restart.rst";
        std::ofstream restartFile;
        restartFile.open(restartFilePath);
        restartFile << out.c_str();
        restartFile.close();
    }
}

void ablate::io::Hdf5MultiFileSerializer::RestoreTS(TS ts) {
    if (resumed) {
        TSSetStepNumber(ts, timeStep);
        TSSetTime(ts, time);
        TSSetTimeStep(ts, dt);
    }
}

std::filesystem::path ablate::io::Hdf5MultiFileSerializer::GetOutputFilePath(const std::string& objectId) const {
    std::stringstream sequenceNumberOutputStream;
    sequenceNumberOutputStream << std::setw(5) << std::setfill('0') << sequenceNumber;
    auto sequenceNumberOutputString = "." + sequenceNumberOutputStream.str();
    return GetOutputDirectoryPath(objectId) / (objectId + sequenceNumberOutputString + extension);
}

std::filesystem::path ablate::io::Hdf5MultiFileSerializer::GetOutputDirectoryPath(const std::string& objectId) { return environment::RunEnvironment::Get().GetOutputDirectory() / objectId; }

#include "registrar.hpp"
REGISTER(ablate::io::Serializer, ablate::io::Hdf5MultiFileSerializer, "serializer for IO that writes each time to a separate hdf5 file",
         ARG(ablate::io::interval::Interval, "interval", "The interval object used to determine write interval."),
         OPT(ablate::parameters::Parameters, "options", "options for the viewer passed directly to PETSc including (hdf5ViewerView, viewer_hdf5_collective, viewer_hdf5_sp_output"));
