#ifndef ABLATELIBRARY_RBF_IMQ_HPP
#define ABLATELIBRARY_RBF_IMQ_HPP

#include "rbf.hpp"

#define __RBF_IMQ_DEFAULT_PARAM 0.1

namespace ablate::domain::rbf {

class IMQ: public RBF {
  private:
    const PetscReal scale = -1;

  public:

    std::string_view type() const override { return "IMQ"; }

    IMQ(PetscInt p = 4, PetscReal scale = 0.1, bool hasDerivatives = false, bool hasInterpolation = false);

    PetscReal RBFVal(PetscReal x[], PetscReal y[]) override;
    PetscReal RBFDer(PetscReal x[], PetscInt dx, PetscInt dy, PetscInt dz) override;

};

}  // namespace ablate::domain::RBF

#endif  // ABLATELIBRARY_RBF_IMQ_HPP
