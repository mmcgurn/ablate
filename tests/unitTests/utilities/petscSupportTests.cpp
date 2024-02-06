#include <petsc.h>
#include <memory>
#include "domain/boxMesh.hpp"
#include "domain/modifiers/distributeWithGhostCells.hpp"
#include "environment/runEnvironment.hpp"
#include "gtest/gtest.h"
#include "mpiTestFixture.hpp"
#include "petscTestErrorChecker.hpp"
#include "utilities/petscSupport.hpp"
#include "utilities/petscUtilities.hpp"

using namespace ablate;

/********************   Begin unit tests for DMPlexGetContainingCell    *************************/

struct RBFSupportParameters_ReturnID {
    testingResources::MpiTestParameter mpiTestParameter;
    std::vector<int> meshFaces;
    std::vector<double> meshStart;
    std::vector<double> meshEnd;
    std::vector<std::shared_ptr<domain::modifiers::Modifier>> meshModifiers;
    bool meshSimplex;
    std::vector<PetscScalar> xyz;
    std::vector<PetscInt> expectedCell;
};

class RBFSupportTestFixture_ReturnID : public testingResources::MpiTestFixture, public ::testing::WithParamInterface<RBFSupportParameters_ReturnID> {
   public:
    void SetUp() override { SetMpiParameters(GetParam().mpiTestParameter); }
};

TEST_P(RBFSupportTestFixture_ReturnID, ShouldReturnCellIDs) {
    StartWithMPI
        {
            // initialize petsc and mpi
            ablate::environment::RunEnvironment::Initialize(argc, argv);
            ablate::utilities::PetscUtilities::Initialize();

            auto testingParam = GetParam();

            // Create the mesh
            // Note that using -dm_view :mesh.tex:ascii_latex -dm_plex_view_scale 10 -dm_plex_view_numbers_depth 1,0,1 will create a mesh, changing numbers_depth as appropriate
            auto mesh = std::make_shared<domain::BoxMesh>("mesh",
                                                          std::vector<std::shared_ptr<domain::FieldDescriptor>>{},
                                                          testingParam.meshModifiers,
                                                          testingParam.meshFaces,
                                                          testingParam.meshStart,
                                                          testingParam.meshEnd,
                                                          std::vector<std::string>{},
                                                          testingParam.meshSimplex);

            PetscInt cell = -2;
            DMPlexGetContainingCell(mesh->GetDM(), &testingParam.xyz[0], &cell) >> utilities::PetscUtilities::checkError;

            PetscMPIInt rank;
            MPI_Comm_rank(PetscObjectComm((PetscObject)mesh->GetDM()), &rank);
            ASSERT_EQ(cell, testingParam.expectedCell[rank]);
        }
        ablate::environment::RunEnvironment::Finalize();
    EndWithMPI
}

INSTANTIATE_TEST_SUITE_P(MeshTests, RBFSupportTestFixture_ReturnID,
                         testing::Values((RBFSupportParameters_ReturnID){.mpiTestParameter = testingResources::MpiTestParameter("2DQuad"),
                                                                         .meshFaces = {10, 5},
                                                                         .meshStart = {0.0, 0.0},
                                                                         .meshEnd = {1.0, 1.0},
                                                                         .meshModifiers = {},
                                                                         .meshSimplex = false,
                                                                         .xyz = {0.55, 0.25},
                                                                         .expectedCell = {15}},
                                         (RBFSupportParameters_ReturnID){.mpiTestParameter = testingResources::MpiTestParameter("2DSimplex"),
                                                                         .meshFaces = {10, 5},
                                                                         .meshStart = {0.0, 0.0},
                                                                         .meshEnd = {1.0, 1.0},
                                                                         .meshModifiers = {},
                                                                         .meshSimplex = true,
                                                                         .xyz = {0.55, 0.25},
                                                                         .expectedCell = {49}},
                                         (RBFSupportParameters_ReturnID){.mpiTestParameter = testingResources::MpiTestParameter("3DQuad"),
                                                                         .meshFaces = {2, 2, 2},
                                                                         .meshStart = {0.0, 0.0, 0.0},
                                                                         .meshEnd = {1.0, 1.0, 1.0},
                                                                         .meshModifiers = {},
                                                                         .meshSimplex = false,
                                                                         .xyz = {0.6, 0.42, 0.8},
                                                                         .expectedCell = {5}},
                                         (RBFSupportParameters_ReturnID){.mpiTestParameter = testingResources::MpiTestParameter("3DSimplex"),
                                                                         .meshFaces = {1, 1, 1},
                                                                         .meshStart = {0.0, 0.0, 0.0},
                                                                         .meshEnd = {2.0, 1.0, 1.0},
                                                                         .meshModifiers = {},
                                                                         .meshSimplex = true,
                                                                         .xyz = {0.1, 0.9, 0.9},
                                                                         .expectedCell = {4}},
                                         (RBFSupportParameters_ReturnID){.mpiTestParameter = testingResources::MpiTestParameter("3DSimplexFail"),
                                                                         .meshFaces = {1, 1, 1},
                                                                         .meshStart = {0.0, 0.0, 0.0},
                                                                         .meshEnd = {2.0, 1.0, 1.0},
                                                                         .meshModifiers = {},
                                                                         .meshSimplex = true,
                                                                         .xyz = {2.1, 0.9, 0.9},
                                                                         .expectedCell = {-1}},
                                         (RBFSupportParameters_ReturnID){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadMPI", 2),
                                                                         .meshFaces = {10, 10},
                                                                         .meshStart = {0.0, 0.0},
                                                                         .meshEnd = {1.0, 1.0},
                                                                         .meshModifiers = {},
                                                                         .meshSimplex = false,
                                                                         .xyz = {0.55, 0.25},
                                                                         .expectedCell = {10, -1}},
                                         (RBFSupportParameters_ReturnID){
                                             .mpiTestParameter = testingResources::MpiTestParameter("2DQuadMPIMod",
                                                                                                    2),  // This is mainly here to check if there is ever a change in how DMLocatePoints functions
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(1)},
                                             .meshSimplex = false,
                                             .xyz = {0.55, 0.25},
                                             .expectedCell = {10, -1}}),
                         [](const testing::TestParamInfo<RBFSupportParameters_ReturnID>& info) { return info.param.mpiTestParameter.getTestName(); });

/********************   End unit tests for DMPlexGetContainingCell    *************************/

/********************   Begin unit tests for DMPlexGetNeighbors    *************************/

struct RBFSupportParameters_NeighborCells {
    testingResources::MpiTestParameter mpiTestParameter;
    std::vector<int> meshFaces;
    std::vector<double> meshStart;
    std::vector<double> meshEnd;
    std::vector<std::shared_ptr<domain::modifiers::Modifier>> meshModifiers;
    bool meshSimplex;
    std::vector<PetscInt> centerCell;
    PetscInt numLevels;
    PetscReal maxDistance;
    PetscInt minNumberCells;
    PetscBool useCells;
    PetscBool returnNeighborVertices;
    std::vector<PetscInt> expectedSizeOfList;
    std::vector<std::vector<PetscInt>> expectedList;
};

class RBFSupportTestFixture_NeighborCells : public testingResources::MpiTestFixture, public ::testing::WithParamInterface<RBFSupportParameters_NeighborCells> {
   public:
    void SetUp() override { SetMpiParameters(GetParam().mpiTestParameter); }
};

TEST_P(RBFSupportTestFixture_NeighborCells, ShouldReturnNeighborCells) {
    StartWithMPI
        {
            // initialize petsc and mpi
            ablate::environment::RunEnvironment::Initialize(argc, argv);
            ablate::utilities::PetscUtilities::Initialize();

            auto testingParam = GetParam();



            // Create the mesh
            // Note that using -dm_view :mesh.tex:ascii_latex -dm_plex_view_scale 10 -dm_plex_view_numbers_depth 1,0,1 will create a mesh, changing numbers_depth as appropriate
            //
            // For debugging testing code use the following:
            //PetscOptionsSetValue(NULL, "-dm_plex_view_scale", "20");
            //PetscOptionsSetValue(NULL, "-dm_plex_view_numbers_depth", "1,0,1");
            //PetscViewer viewer;
            //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "/path/to/file/mesh.tex", &viewer);
            //PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_LATEX);
            //DMView(mesh->GetDM(), viewer);
            //PetscViewerPopFormat(viewer);
            //PetscViewerDestroy(&viewer);
            auto mesh = std::make_shared<domain::BoxMesh>("mesh",
                                                          std::vector<std::shared_ptr<domain::FieldDescriptor>>{},
                                                          testingParam.meshModifiers,
                                                          testingParam.meshFaces,
                                                          testingParam.meshStart,
                                                          testingParam.meshEnd,
                                                          std::vector<std::string>{},
                                                          testingParam.meshSimplex);



            PetscInt nCells, *cells;
            PetscMPIInt rank;
            MPI_Comm_rank(PetscObjectComm((PetscObject)mesh->GetDM()), &rank);

            DMPlexGetNeighbors(mesh->GetDM(),
                               testingParam.centerCell[rank],
                               testingParam.numLevels,
                               testingParam.maxDistance,
                               testingParam.minNumberCells,
                               testingParam.useCells,
                               testingParam.returnNeighborVertices,
                               &nCells,
                               &cells) >>
                utilities::PetscUtilities::checkError;

            PetscSortInt(nCells, cells);
            PetscSortInt(testingParam.expectedSizeOfList[rank], testingParam.expectedList[rank].data());  // Should probably enter as a sorted list. Leaving for later.

            ASSERT_EQ(nCells, testingParam.expectedSizeOfList[rank]);



            // There may be a better way of doing this, but with DMPlexGetNeighbors sticking with C-only code there may not be.
            // Also note that as cells is a dynamically allocated array there is not way (that I know of) to get the number of elements.
            for (int i = 0; i < nCells; ++i) {
                ASSERT_EQ(cells[i], testingParam.expectedList[rank][i]);
            }

            // Restore the neighbors
            DMPlexRestoreNeighbors(mesh->GetDM(),
                                   testingParam.centerCell[rank],
                                   testingParam.numLevels,
                                   testingParam.maxDistance,
                                   testingParam.minNumberCells,
                                   testingParam.useCells,
                                   testingParam.returnNeighborVertices,
                                   &nCells,
                                   &cells) >>
                utilities::PetscUtilities::checkError;

            PetscFree(cells);
        }
        ablate::environment::RunEnvironment::Finalize();
    EndWithMPI
}

INSTANTIATE_TEST_SUITE_P(
    MeshTests, RBFSupportTestFixture_NeighborCells,
    testing::Values(
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadVert"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {25},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 25,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {25},
                                             .expectedList = {{25, 24, 35, 15, 26, 14, 34, 36, 16, 23, 45, 5, 27, 13, 33, 44, 46, 4, 37, 6, 17, 43, 3, 47, 7}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadVertCorner"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {0},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 25,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {25},
                                             .expectedList = {{0, 10, 1, 11, 2, 20, 12, 21, 22, 30, 3, 31, 13, 23, 32, 4, 40, 14, 41, 33, 42, 24, 43, 34, 44}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriVert"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {199},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 25,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {25},
                                             .expectedList = {{199, 76, 159, 79, 150, 80, 149, 98, 111, 73, 40, 78, 158, 112, 75, 109, 45, 81, 154, 95, 151, 82, 156, 152, 72}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriVertCorner"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {0},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 25,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {25},
                                             .expectedList = {{0, 6, 1, 4, 2, 3, 7, 19, 5, 9, 21, 12, 22, 8, 18, 25, 14, 11, 23, 13, 24, 47, 10, 30, 27}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriVertNoOverlap", 2),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {56, 19},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 10,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {10, 10},
                                             .expectedList = {{56, 60, 57, 40, 55, 45, 58, 71, 41, 54}, {19, 21, 17, 22, 16, 23, 35, 18, 20, 102}}},
        (RBFSupportParameters_NeighborCells){
            .mpiTestParameter = testingResources::MpiTestParameter("2DQuadVertOverlap", 4),
            .meshFaces = {10, 10},
            .meshStart = {0.0, 0.0},
            .meshEnd = {1.0, 1.0},
            .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(1)},
            .meshSimplex = false,
            .centerCell = {24, 4, 20, 0},
            .numLevels = -1,
            .maxDistance = -1.0,
            .minNumberCells = 9,
            .useCells = PETSC_FALSE,
            .returnNeighborVertices = PETSC_FALSE,
            .expectedSizeOfList = {9, 9, 9, 9},
            .expectedList = {{24, 23, 19, 29, 35, 18, 28, 34, 30}, {4, 9, 31, 3, 29, 8, 28, 30, 32}, {20, 21, 31, 15, 29, 16, 28, 30, 32}, {0, 31, 26, 5, 1, 6, 27, 32, 25}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadEdge", 1),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = false,
                                             .centerCell = {54},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 9,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {9},
                                             .expectedList = {{54, 64, 53, 55, 44, 63, 43, 45, 65}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriEdge", 1),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = true,
                                             .centerCell = {199},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 9,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {9},
                                             .expectedList = {{199, 76, 159, 79, 150, 149, 80, 98, 111}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadEdgeOverlap", 2),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(1)},
                                             .meshSimplex = false,
                                             .centerCell = {11, 34},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 9,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {9, 9},
                                             .expectedList = {{11, 12, 6, 16, 10, 17, 7, 15, 5}, {34, 33, 29, 39, 56, 38, 28, 57, 55}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("3DQuadFace", 1),
                                             .meshFaces = {4, 4, 4},
                                             .meshStart = {0.0, 0.0, 0.0},
                                             .meshEnd = {1.0, 1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = false,
                                             .centerCell = {25},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 20,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {20},
                                             .expectedList = {{25, 24, 9, 26, 21, 29, 41, 8, 13, 28, 10, 5, 20, 40, 22, 45, 30, 37, 42, 27}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("3DTriFace", 1),
                                             .meshFaces = {4, 4, 4},
                                             .meshStart = {0.0, 0.0, 0.0},
                                             .meshEnd = {1.0, 1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = true,
                                             .centerCell = {25},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 20,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {20},
                                             .expectedList = {{25, 41, 38, 33, 28, 51, 39, 46, 56, 17, 15, 4, 32, 198, 95, 57, 201, 123, 40, 74}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadDistanceEdge"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = false,
                                             .centerCell = {55},
                                             .numLevels = -1,
                                             .maxDistance = 0.28,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {21},
                                             .expectedList = {{34, 35, 36, 43, 44, 45, 46, 47, 53, 54, 55, 56, 57, 63, 64, 65, 66, 67, 74, 75, 76}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadDistanceVert"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = false,
                                             .centerCell = {55},
                                             .numLevels = -1,
                                             .maxDistance = 0.28,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {21},
                                             .expectedList = {{34, 35, 36, 43, 44, 45, 46, 47, 53, 54, 55, 56, 57, 63, 64, 65, 66, 67, 74, 75, 76}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadDistanceEdgeMPI", 2),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(3)},
                                             .meshSimplex = false,
                                             .centerCell = {25, 29},
                                             .numLevels = -1,
                                             .maxDistance = 0.28,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {21, 21},
                                             .expectedList = {{15, 16, 20, 21, 22, 25, 26, 27, 30, 31, 32, 35, 36, 61, 63, 64, 66, 67, 69, 70, 73},
                                                              {18, 19, 22, 23, 24, 27, 28, 29, 32, 33, 34, 38, 39, 59, 62, 63, 65, 66, 68, 69, 71}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadDistanceVertMPI", 2),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(3)},
                                             .meshSimplex = false,
                                             .centerCell = {25, 29},
                                             .numLevels = -1,
                                             .maxDistance = 0.28,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {21, 21},
                                             .expectedList = {{15, 16, 20, 21, 22, 25, 26, 27, 30, 31, 32, 35, 36, 61, 63, 64, 66, 67, 69, 70, 73},
                                                              {18, 19, 22, 23, 24, 27, 28, 29, 32, 33, 34, 38, 39, 59, 62, 63, 65, 66, 68, 69, 71}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriDistanceEdge"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = true,
                                             .centerCell = {199},
                                             .numLevels = -1,
                                             .maxDistance = 0.14,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {11},
                                             .expectedList = {{40, 73, 76, 79, 80, 98, 111, 149, 150, 159, 199}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriDistanceVert"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = true,
                                             .centerCell = {199},
                                             .numLevels = -1,
                                             .maxDistance = 0.14,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {11},
                                             .expectedList = {{40, 73, 76, 79, 80, 98, 111, 149, 150, 159, 199}}},

        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriDistanceEdgeMPI", 2),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(3)},
                                             .meshSimplex = true,
                                             .centerCell = {60, 102},
                                             .numLevels = -1,
                                             .maxDistance = 0.14,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {12, 11},
                                             .expectedList = {{38, 40, 41, 45, 56, 58, 60, 71, 135, 136, 154, 163}, {24, 38, 52, 53, 62, 102, 112, 123, 124, 125, 128}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriDistanceVertMPI", 2),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(3)},
                                             .meshSimplex = true,
                                             .centerCell = {60, 102},
                                             .numLevels = -1,
                                             .maxDistance = 0.14,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {12, 11},
                                             .expectedList = {{38, 40, 41, 45, 56, 58, 60, 71, 135, 136, 154, 163}, {24, 38, 52, 53, 62, 102, 112, 123, 124, 125, 128}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadLevelVert"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(3)},
                                             .meshSimplex = false,
                                             .centerCell = {55},
                                             .numLevels = 2,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {25},
                                             .expectedList = {{33, 34, 35, 36, 37, 43, 44, 45, 46, 47, 53, 54, 55, 56, 57, 63, 64, 65, 66, 67, 73, 74, 75, 76, 77}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadLevelEdge"),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(3)},
                                             .meshSimplex = false,
                                             .centerCell = {55},
                                             .numLevels = 2,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_FALSE,
                                             .expectedSizeOfList = {13},
                                             .expectedList = {{35, 44, 45, 46, 53, 54, 55, 56, 57, 64, 65, 66, 75}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadCenterLevelCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {12},
                                             .numLevels = 1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {16},
                                             .expectedList = {{32, 33, 34, 35, 38, 39, 40, 41, 44, 45, 46, 47, 50, 51, 52, 53}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadCenterLevelRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {12},
                                             .numLevels = 1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {12},
                                             .expectedList = {{33, 34, 38, 39, 40, 41, 44, 45, 46, 47, 51, 52}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadLevelCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {6},
                                             .numLevels = 2,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {25},
                                             .expectedList = {{25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadLevelRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {6},
                                             .numLevels = 2,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {20},
                                             .expectedList = {{25, 26, 27, 28, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 43, 44, 45, 46, 50, 51}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadCornerLevelCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {24},
                                             .numLevels = 1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {9},
                                             .expectedList = {{46, 47, 48, 52, 53, 54, 58, 59, 60}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadCornerLevelRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {24},
                                             .numLevels = 1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {8},
                                             .expectedList = {{47, 48, 52, 53, 54, 58, 59, 60}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {11},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 16,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {16},
                                             .expectedList = {{31, 32, 33, 34, 37, 38, 39, 40, 43, 44, 45, 46, 49, 50, 51, 52}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {3},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 10,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {10},
                                             .expectedList = {{27, 28, 29, 30, 33, 34, 35, 36, 40, 41}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadCornerDistanceCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {20},
                                             .numLevels = -1,
                                             .maxDistance = 0.4,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {8},
                                             .expectedList = {{43, 44, 49, 50, 51, 55, 56, 57}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadCenterDistanceCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = false,
                                             .centerCell = {12},
                                             .numLevels = -1,
                                             .maxDistance = 0.3,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {4},
                                             .expectedList = {{39, 40, 45, 46}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriLevelRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {24},
                                             .numLevels = 1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {13},
                                             .expectedList = {{58, 59, 60, 63, 64, 65, 66, 69, 70, 71, 72, 75, 76}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriLevelCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {24},
                                             .numLevels = 1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {13},
                                             .expectedList = {{58, 59, 60, 63, 64, 65, 66, 69, 70, 71, 72, 75, 76}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriEdgeLevelCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {13},
                                             .numLevels = 2,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {19},
                                             .expectedList = {{50, 51, 52, 56, 57, 58, 62, 63, 64, 65, 68, 69, 70, 71, 74, 75, 76, 80, 81}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriCornerDistanceRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {5.0, 5.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {0},
                                             .numLevels = -1,
                                             .maxDistance = 2.5,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {9},
                                             .expectedList = {{50, 51, 52, 56, 57, 58, 62, 63, 64}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriCellRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {5.0, 5.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {41},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 13,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {13},
                                             .expectedList = {{65, 66, 70, 71, 72, 73, 76, 77, 78, 79, 83, 84, 85}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriCornerRV"),
                                             .meshFaces = {5, 5},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {5.0, 5.0},
                                             .meshModifiers = {},
                                             .meshSimplex = true,
                                             .centerCell = {31},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 10,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {10},
                                             .expectedList = {{52, 53, 54, 55, 59, 60, 61, 65, 66, 67}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("3DQuadCornerRV"),
                                             .meshFaces = {2, 2, 2},
                                             .meshStart = {0.0, 0.0, 0.0},
                                             .meshEnd = {2.0, 2.0, 2.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = false,
                                             .centerCell = {0},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 20,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {20},
                                             .expectedList = {{8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 29, 30}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadEdgeOverlapRVMPI", 2),
                                             .meshFaces = {8, 8},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {1.0, 1.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{},
                                             .meshSimplex = false,
                                             .centerCell = {12, 15},
                                             .numLevels = 1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {10, 10},
                                             .expectedList = {{42, 43, 47, 48, 49, 52, 53, 54, 57, 58}, {45, 46, 49, 50, 51, 54, 55, 56, 60, 61}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadLevelEdgeMPI", 2),
                                             .meshFaces = {8, 8},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {8.0, 8.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(2)},
                                             .meshSimplex = false,
                                             .centerCell = {6, 5},
                                             .numLevels = 1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {12, 12},
                                             .expectedList = {{49, 50, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62}, {49, 50, 53, 54, 55, 56, 58, 59, 60, 61, 64, 65}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadDistanceEdgeRVMPI", 2),
                                             .meshFaces = {10, 10},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {10.0, 10.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(3)},
                                             .meshSimplex = false,
                                             .centerCell = {18, 31},
                                             .numLevels = -1,
                                             .maxDistance = 2.8,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_FALSE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {22, 22},
                                             .expectedList = {{87, 88, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 106, 107, 108, 109, 112, 113},
                                                              {105, 106, 110, 111, 112, 113, 116, 117, 118, 119, 120, 122, 123, 124, 125, 126, 128, 129, 130, 131, 135, 136}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadDistanceCellRVMPI", 2),
                                             .meshFaces = {8, 8},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {8.0, 8.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(1)},
                                             .meshSimplex = false,
                                             .centerCell = {26, 5},
                                             .numLevels = -1,
                                             .maxDistance = 2.3,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {16, 16},
                                             .expectedList = {{60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75}, {40, 41, 42, 43, 45, 46, 47, 48, 50, 51, 52, 53, 55, 56, 57, 58}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DQuadCellRVMPI", 2),
                                             .meshFaces = {8, 8},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {8.0, 8.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(1)},
                                             .meshSimplex = false,
                                             .centerCell = {26, 5},
                                             .numLevels = -1,
                                             .maxDistance = -1.0,
                                             .minNumberCells = 16,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {16, 16},
                                             .expectedList = {{60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75}, {40, 41, 42, 43, 45, 46, 47, 48, 50, 51, 52, 53, 55, 56, 57, 58}}},
        (RBFSupportParameters_NeighborCells){.mpiTestParameter = testingResources::MpiTestParameter("2DTriDistanceEdgeRVMPI", 2),
                                             .meshFaces = {6, 6},
                                             .meshStart = {0.0, 0.0},
                                             .meshEnd = {6.0, 6.0},
                                             .meshModifiers = std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<domain::modifiers::DistributeWithGhostCells>(1)},
                                             .meshSimplex = true,
                                             .centerCell = {17, 15},
                                             .numLevels = 1,
                                             .maxDistance = -1,
                                             .minNumberCells = -1,
                                             .useCells = PETSC_TRUE,
                                             .returnNeighborVertices = PETSC_TRUE,
                                             .expectedSizeOfList = {11, 11},
                                             .expectedList = {{51, 52, 53, 54, 58, 59, 60, 65, 66, 81, 82}, {51, 52, 57, 58, 59, 63, 64, 65, 66, 71, 72}}}),
    [](const testing::TestParamInfo<RBFSupportParameters_NeighborCells>& info) { return info.param.mpiTestParameter.getTestName(); });

struct RBFSupportParameters_ErrorChecking {
    testingResources::MpiTestParameter mpiTestParameter;
};

class RBFSupportTestFixture_ErrorChecking : public testingResources::MpiTestFixture, public ::testing::WithParamInterface<RBFSupportParameters_ErrorChecking> {
   public:
    void SetUp() override { SetMpiParameters(GetParam().mpiTestParameter); }
};

TEST_P(RBFSupportTestFixture_ErrorChecking, ShouldThrowErrorForTooManyInputs) {
    StartWithMPI
        {
            // initialize petsc and mpi
            ablate::environment::RunEnvironment::Initialize(argc, argv);
            ablate::utilities::PetscUtilities::Initialize();

            // Create the mesh
            auto mesh = std::make_shared<domain::BoxMesh>("mesh",
                                                          std::vector<std::shared_ptr<domain::FieldDescriptor>>{},
                                                          std::vector<std::shared_ptr<domain::modifiers::Modifier>>{},
                                                          std::vector<int>{2, 2},
                                                          std::vector<double>{0.0, 0.0},
                                                          std::vector<double>{1.0, 1.0});

            PetscInt nCells, *cells;

            EXPECT_ANY_THROW(DMPlexGetNeighbors(mesh->GetDM(), 0, 1, 1.0, 1, PETSC_TRUE, PETSC_FALSE, &nCells, &cells) >> utilities::PetscUtilities::checkError);
        }
        ablate::environment::RunEnvironment::Finalize();
    EndWithMPI
}

INSTANTIATE_TEST_SUITE_P(MeshTests, RBFSupportTestFixture_ErrorChecking,
                         testing::Values((RBFSupportParameters_ErrorChecking){.mpiTestParameter = testingResources::MpiTestParameter("SingleProc")},
                                         (RBFSupportParameters_ErrorChecking){.mpiTestParameter = testingResources::MpiTestParameter("DualProcs", 2)}),
                         [](const testing::TestParamInfo<RBFSupportParameters_ErrorChecking>& info) { return info.param.mpiTestParameter.getTestName(); });

/********************   End unit tests for DMPlexGetNeighbors    *************************/
