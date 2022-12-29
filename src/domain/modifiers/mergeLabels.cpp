#include "mergeLabels.hpp"
#include "utilities/petscUtilities.hpp"

ablate::domain::modifiers::MergeLabels::MergeLabels(std::shared_ptr<domain::Region> mergedRegion, std::vector<std::shared_ptr<domain::Region>> regions)
    : mergedRegion(mergedRegion), regions(regions) {}

void ablate::domain::modifiers::MergeLabels::Modify(DM& dm) {
    // Create a new label for the merged region
    DMCreateLabel(dm, mergedRegion->GetName().c_str()) >> utilities::PetscUtilities::checkError;

    // Get the label for each region
    std::vector<IS> regionISs(regions.size(), nullptr);

    // Get the data
    for (std::size_t r = 0; r < regions.size(); r++) {
        auto& regionIS = regionISs[r];
        DMGetStratumIS(dm, regions[r]->GetName().c_str(), regions[r]->GetValue(), &regionIS) >> utilities::PetscUtilities::checkError;
    }
    // Create Concatenate IS
    IS mergedIS;
    ISConcatenate(PETSC_COMM_SELF, regionISs.size(), &regionISs[0], &mergedIS) >> utilities::PetscUtilities::checkError;
    ISSortRemoveDups(mergedIS) >> utilities::PetscUtilities::checkError;

    // cleanup
    for (auto& is : regionISs) {
        ISDestroy(&is) >> utilities::PetscUtilities::checkError;
    }

    DMLabel mergedLabel;
    DMGetLabel(dm, mergedRegion->GetName().c_str(), &mergedLabel) >> utilities::PetscUtilities::checkError;
    DMLabelSetStratumIS(mergedLabel, mergedRegion->GetValue(), mergedIS) >> utilities::PetscUtilities::checkError;

    // cleanup
    ISDestroy(&mergedIS) >> utilities::PetscUtilities::checkError;
    DMPlexLabelComplete(dm, mergedLabel) >> utilities::PetscUtilities::checkError;
}

#include "registrar.hpp"
REGISTER(ablate::domain::modifiers::Modifier, ablate::domain::modifiers::MergeLabels, "Creates a new label for all faces on the outside of the boundary",
         ARG(ablate::domain::Region, "mergedRegion", "the merged region to create"), ARG(std::vector<ablate::domain::Region>, "regions", "the regions to include in the new merged region"));
