#ifndef ABLATELIBRARY_ADAPTPHYSICS_HPP
#define ABLATELIBRARY_ADAPTPHYSICS_HPP

#include <petsc.h>
#include <petsc/private/tsimpl.h>
#include "utilities/petscUtilities.hpp"
namespace ablate::solver {

/**
 * This is a static class that holds the required code needed by petsc for the adapt physics time step constraint
 */
class AdaptPhysics {
   private:
    static inline const char name[] = "physics";

    /**
     * The physics based implementation of TSAdaptChoose
     * @param adapt
     * @param ts
     * @param h
     * @param next_sc
     * @param next_h
     * @param accept
     * @param wlte
     * @param wltea
     * @param wlter
     * @return
     */
    static PetscErrorCode TSAdaptChoose(TSAdapt adapt, TS ts, PetscReal h, PetscInt *next_sc, PetscReal *next_h, PetscBool *accept, PetscReal *wlte, PetscReal *wltea, PetscReal *wlter);

    /**
     * Static call to modify the adapt petsc object into an adapt physics implementation
     * @param adapt
     * @return
     */
    static PetscErrorCode TSAdaptCreate(TSAdapt adapt) {
        PetscFunctionBegin;
        adapt->ops->choose = TSAdaptChoose;

        // Set default behavior needed by the physics adapt
        PetscCall(TSAdaptSetAlwaysAccept(adapt, PETSC_TRUE));

        PetscFunctionReturn(0);
    }

   public:
    /**
     * Function to register the ts adapt with petsc
     */
    static void Register() { TSAdaptRegister(name, TSAdaptCreate) >> ablate::utilities::PetscUtilities::checkError; }

    /**
     * Prevent this class from being used in an non static way
     */
    AdaptPhysics() = delete;
};
}  // namespace ablate::solver

#endif  // ABLATELIBRARY_ADAPTPHYSICS_HPP
