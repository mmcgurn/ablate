#include "rbfV2.hpp"
//#include <petsc/private/dmpleximpl.h>
//#include "utilities/mpiError.hpp"
//#include "utilities/petscError.hpp"


template <typename Enumeration>
constexpr auto as_integer(Enumeration const value)
    -> typename std::underlying_type<Enumeration>::type
{
    static_assert(std::is_enum<Enumeration>::value, "parameter is not of type enum or enum class");
    return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}




static const std::map<std::string, ablate::radialBasisV2::RBFType> stringToRBFType = {{"", ablate::radialBasisV2::RBFType::mq},
                                                                                           {"phs", ablate::radialBasisV2::RBFType::phs},
                                                                                           {"PHS", ablate::radialBasisV2::RBFType::phs},
                                                                                           {"mq", ablate::radialBasisV2::RBFType::mq},
                                                                                           {"MQ", ablate::radialBasisV2::RBFType::mq},
                                                                                           {"imq", ablate::radialBasisV2::RBFType::imq},
                                                                                           {"IMQ", ablate::radialBasisV2::RBFType::imq},
                                                                                           {"ga", ablate::radialBasisV2::RBFType::ga},
                                                                                           {"GA", ablate::radialBasisV2::RBFType::ga}};

std::istream& ablate::radialBasisV2::operator>>(std::istream& is, ablate::radialBasisV2::RBFType& v) {
    std::string enumString;
    is >> enumString;
    v = stringToRBFType.at(enumString);
    return is;
}


using namespace ablate::radialBasisV2;


// Distance between two points
static PetscReal DistanceSquared(PetscInt dim, PetscReal x[], PetscReal y[]){
  PetscReal r = 0.0;
  for (PetscInt d = 0; d < dim; ++d) {
    r += (x[d] - y[d]) * (x[d] - y[d]);
  }

  return r;
}

// Distance between point and the origin
static PetscReal DistanceSquared(PetscInt dim, PetscReal x[]){
  PetscReal r = 0.0;
  for (PetscInt d = 0; d < dim; ++d) {
    r += x[d] * x[d];
  }
  return r;
}


/************ Begin Polyharmonic  **********************/

 //Polyharmonic spline: r^m
static PetscReal phsVal(const PetscInt dim, PetscReal x[], PetscReal y[], const PetscReal m) {
  PetscReal r = DistanceSquared(dim, x, y);

  return PetscPowReal(r, 0.5*m);
}

// Derivatives of Polyharmonic spline at a location.
static PetscReal phsDer(const PetscInt dim, PetscReal x[], PetscInt dx, PetscInt dy, PetscInt dz, const PetscReal m) {
  PetscReal r = DistanceSquared(dim, x);

  r = PetscSqrtReal(r);

  switch (dx + 10*dy + 100*dz) {
    case 0:
      r = PetscPowReal(r, m);
      break;
    case 1: // x
      r = -m*x[0]*PetscPowReal(r, (m-2));
      break;
    case 2: // xx
      r = m*PetscPowReal(r, (m-2)) + m*(m-2)*x[0]*x[0]*PetscPowReal(r, (m-4));
      break;
    case 10: // y
      r = -m*x[1]*PetscPowReal(r, (m-2));
      break;
    case 20: // yy
      r = m*PetscPowReal(r, (m-2)) + m*(m-2)*x[1]*x[1]*PetscPowReal(r, (m-4));
      break;
    case 100: // z
      r = -m*x[2]*PetscPowReal(r, (m-2));
      break;
    case 200: // zz
      r = m*PetscPowReal(r, (m-2)) + m*(m-2)*x[2]*x[2]*PetscPowReal(r, (m-4));
      break;
    case 11: // xy
      r = m*(m-2)*x[0]*x[1]*PetscPowReal(r, (m-4));
      break;
    case 101: // xz
      r = m*(m-2)*x[0]*x[2]*PetscPowReal(r, (m-4));
      break;
    case 110: // yz
      r = m*(m-2)*x[1]*x[2]*PetscPowReal(r, (m-4));
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown derivative!\n");
  }

  return r;
}
/************ End Polyharmonic  **********************/



/************ Begin Multiquadric **********************/

// Multiquadric: sqrt(1+(er)^2)
static PetscReal mqVal(const PetscInt dim, PetscReal x[], PetscReal y[], const PetscReal h) {

  PetscReal e = 1.0/h;
  PetscReal r = DistanceSquared(dim, x, y);

  return PetscSqrtReal(1.0 + e*e*r);
}

// Derivatives of Multiquadric spline at a location.
static PetscReal mqDer(const PetscInt dim, PetscReal x[], PetscInt dx, PetscInt dy, PetscInt dz, const PetscReal h) {

  PetscReal e = 1.0/h;
  PetscReal r = DistanceSquared(dim, x);

  r = PetscSqrtReal(1.0 + e*e*r);

  switch (dx + 10*dy + 100*dz) {
    case 0:
      // Do nothing
      break;
    case 1: // x
      r = -e*e*x[0]/r;
      break;
    case 2: // xx
      r = e*e*(1.0 + e*e*(x[1]*x[1] + x[2]*x[2]))/PetscPowReal(r, 3.0);
      break;
    case 10: // y
      r = -e*e*x[1]/r;
      break;
    case 20: // yy
      r = e*e*(1.0 + e*e*(x[0]*x[0] + x[2]*x[2]))/PetscPowReal(r, 3.0);
      break;
    case 100: // z
      r = -e*e*x[2]/r;
      break;
    case 200: // zz
      r = e*e*(1.0 + e*e*(x[0]*x[0] + x[1]*x[1]))/PetscPowReal(r, 3.0);
      break;
    case 11: // xy
      r = -PetscSqr(e*e)*x[0]*x[1]/PetscPowReal(r, 3.0);
      break;
    case 101: // xz
      r = -PetscSqr(e*e)*x[0]*x[2]/PetscPowReal(r, 3.0);
      break;
    case 110: // yz
      r = -PetscSqr(e*e)*x[1]*x[2]/PetscPowReal(r, 3.0);
      break;
    default:
      throw std::invalid_argument("Derivative of (" + std::to_string(dx) + ", " + std::to_string(dy) + ", " + std::to_string(dz) + ") is not setup.");
  }

  return r;
}

/************ End Multiquadric **********************/

/************ Begin Inverse Multiquadric **********************/
// Multiquadric: sqrt(1+(er)^2)
static PetscReal imqVal(const PetscInt dim, PetscReal x[], PetscReal y[], const PetscReal h) {

  PetscReal e = 1.0/h;
  PetscReal r = DistanceSquared(dim, x, y);

  return 1.0/PetscSqrtReal(1.0 + e*e*r);
}

// Derivatives of Multiquadric spline at a location.
static PetscReal imqDer(const PetscInt dim, PetscReal x[], PetscInt dx, PetscInt dy, PetscInt dz, const PetscReal h) {

  PetscReal e = 1.0/h;
  PetscReal r = DistanceSquared(dim, x);

  r = PetscSqrtReal(1.0 + e*e*r);

  switch (dx + 10*dy + 100*dz) {
    case 0:
      r = 1.0/r;
      break;
    case 1: // x
      r = -e*e*x[0]/PetscPowReal(r, 3.0);
      break;
    case 2: // xx
      r = -e*e*(1.0 + e*e*(-2.0*x[0]*x[0] + x[1]*x[1] + x[2]*x[2]))/PetscPowReal(r, 5.0);
      break;
    case 10: // y
      r = -e*e*x[1]/PetscPowReal(r, 3.0);
      break;
    case 20: // yy
      r = -e*e*(1.0 + e*e*(x[0]*x[0] - 2.0*x[1]*x[1] + x[2]*x[2]))/PetscPowReal(r, 5.0);
      break;
    case 100: // z
      r = -e*e*x[2]/PetscPowReal(r, 3.0);
      break;
    case 200: // zz
      r = -e*e*(1.0 + e*e*(x[0]*x[0] + x[1]*x[1] - 2.0*x[2]*x[2]))/PetscPowReal(r, 5.0);
      break;
    case 11: // xy
      r = 3.0*PetscSqr(e*e)*x[0]*x[1]/PetscPowReal(r, 5.0);
      break;
    case 101: // xz
      r = 3.0*PetscSqr(e*e)*x[0]*x[2]/PetscPowReal(r, 5.0);
      break;
    case 110: // yz
      r = 3.0*PetscSqr(e*e)*x[1]*x[2]/PetscPowReal(r, 5.0);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown derivative!\n");
  }

  return r;
}
/************ End Inverse Multiquadric Derived Class **********************/

/************ Begin Gaussian **********************/
// Gaussian: r^m
static PetscReal gaVal(const PetscInt dim, PetscReal x[], PetscReal y[], const PetscReal h) {

  PetscReal e2 = 1.0/(h*h);
  PetscReal r = DistanceSquared(dim, x, y);

  return PetscExpReal(-r*e2);
}

// Derivatives of Gaussian spline at a location.
static PetscReal gaDer(const PetscInt dim, PetscReal x[], PetscInt dx, PetscInt dy, PetscInt dz, const PetscReal h) {

  PetscReal e2 = 1.0/(h*h);
  PetscReal r = DistanceSquared(dim, x);

  r = PetscExpReal(-r*e2);

  switch (dx + 10*dy + 100*dz) {
    case 0:
      // Do nothing
      break;
    case 1: // x
      r *= -2.0*e2*x[0];
      break;
    case 2: // xx
      r *= 2.0*e2*(2.0*e2*x[0]*x[0]-1.0);
      break;
    case 10: // x[1]
      r *= -2.0*e2*x[1];
      break;
    case 20: // yy
      r *= 2.0*e2*(2.0*e2*x[1]*x[1]-1.0);
      break;
    case 100: // x[2]
      r *= -2.0*e2*x[2];
      break;
    case 200: // zz
      r *= 2.0*e2*(2.0*e2*x[2]*x[2]-1.0);
      break;
    case 11: // xy
      r *= 4.0*e2*e2*x[0]*x[1];
      break;
    case 101: // xz
      r *= 4.0*e2*e2*x[0]*x[2];
      break;
    case 110: // yz
      r *= 4.0*e2*e2*x[1]*x[2];
      break;
    case 111:
      r *= 8.0*e2*e2*e2*x[1]*x[2];
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown derivative!\n");
  }

  return r;
}
/************ End Gaussian **********************/

static PetscInt fac[11] =  {1,1,2,6,24,120,720,5040,40320,362880,3628800}; // Pre-computed factorials

// Compute the LU-factorization of the augmented RBF matrix
// c - Location in cellRange
// xCenters - Shifted centers of all the cells (output)
// A - The LU-factorization of the augmented RBF matrix (output)
void RBF::Matrix(const PetscInt c, PetscReal **xCenters, Mat *LUA) {

  PetscInt              i, j, d, px, py, pz, matSize;
  PetscInt              nCells, *list;
  const PetscInt        dim = subDomain->GetDimensions();
  const PetscInt        nPoly = RBF::nPoly, p = RBF::polyOrder, p1 = p + 1;
  const PetscReal       param = RBF::rbfParam;
  Mat                   A;
  PetscReal             *x, *vals;
  DM                    dm = subDomain->GetDM(); // The underlying mesh
  PetscReal             x0[dim];        // Center of the cell of interest
  PetscReal             *xp; // Powers of the cell centers
  ablate::solver::Range cellRange;

  GetCellRange(cellRange);

  const PetscInt cell = cellRange.points ? cellRange.points[c] : c;  //The cell in PETSc ordering

  // Get the list of neighbor cells
  DMPlexGetNeighborCells(dm, cell, -1, -1.0, RBF::minNumberCells, RBF::useVertices, &nCells, &list);
  RBF::nStencil[c] = nCells;
  RBF::stencilList[c] = list;

  PetscMalloc1(nCells*dim*p1, &xp) >> ablate::checkError;


  if(nPoly>=nCells){
    throw std::invalid_argument("Number of surrounding cells, " + std::to_string(nCells) + ", can not support a requested polynomial order of " + std::to_string(p) + " which requires " + std::to_string(nPoly) + " number of cells.");
  }

  // Get the cell center
  DMPlexComputeCellGeometryFVM(dm, cell, NULL, x0, NULL) >> ablate::checkError;

  // Shifted cell-centers of neighbor cells
  PetscMalloc1(nCells*dim, &x);
  for (i = 0; i < nCells; ++i) {
    DMPlexComputeCellGeometryFVM(dm, list[i], NULL, &x[i*dim], NULL) >> ablate::checkError;
    for (d = 0; d < dim; ++d) {
      x[i*dim+d] -= x0[d];
      // Precompute the powers for later use
      xp[(i*dim+d)*p1 + 0] = 1.0;
      for (px = 1; px < p+1; ++px) {
        xp[(i*dim+d)*p1 + px] = xp[(i*dim+d)*p1 + (px-1)]*x[i*dim+d];
      }
    }
  }

  matSize = nCells + nPoly;

  // Create the matrix
  MatCreateSeqDense(PETSC_COMM_SELF, matSize, matSize, NULL, &A) >> ablate::checkError;
  PetscObjectSetName((PetscObject)A,"ablate::radialBasisV2::RBF::A") >> ablate::checkError;
  MatZeroEntries(A) >> ablate::checkError;
  MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE) >> ablate::checkError;

  MatDenseGetArrayWrite(A, &vals) >> ablate::checkError;

  // RBF contributions to the matrix
  for (i = 0; i < nCells; ++i) {
    for (j = i; j < nCells; ++j) {
      vals[i*matSize + j] = vals[j*matSize + i] = RBFVal(dim, &x[i*dim], &x[j*dim], param);
    }
  }

  // Augmented polynomial contributions
  if (dim == 2) {
    for (i = 0; i < nCells; ++i) {
      j = nCells;
      for (py = 0; py < p1; ++py) {
        for (px = 0; px < p1-py; ++px) {
          vals[i*matSize + j] = vals[j*matSize + i] = xp[(i*dim+0)*p1 + px] * xp[(i*dim+1)*p1 + py];
          ++j;
        }
      }
    }
  } else {
    for (i = 0; i < nCells; ++i) {
      j = nCells;
      for (pz = 0; pz < p1; ++pz) {
        for (py = 0; py < p1-pz; ++py) {
          for (px = 0; px < p1-py-pz; ++px ){
            vals[i*matSize + j] = vals[j*matSize + i] = xp[(i*dim+0)*p1 + px] * xp[(i*dim+1)*p1 + py] * xp[(i*dim+2)*p1 + pz];
            ++j;
          }
        }
      }
    }
  }
  MatDenseRestoreArrayWrite(A, &vals) >> ablate::checkError;
  MatViewFromOptions(A,NULL,"-ablate::radialBasisV2::RBF::A_view") >> ablate::checkError;

  // Factor the matrix
  MatLUFactor(A, NULL, NULL, NULL) >> ablate::checkError;

  PetscFree(xp) >> ablate::checkError;

  // Assign output
  *xCenters = x;
  *LUA = A;

}


/************ Begin Derivative Code **********************/


// Set the derivatives to use
// nDer - Number of derivatives to set
// dx, dy, dz - Lists of length nDer indicating the derivatives
void RBF::SetDerivatives(PetscInt nDer, PetscInt dx[], PetscInt dy[], PetscInt dz[], PetscBool useVertices){

  if (nDer > 0) {

    ablate::solver::Range cellRange;
    GetCellRange(cellRange);

    PetscInt n = cellRange.end - cellRange.start;

    RBF::hasDerivativeInformation = PETSC_TRUE;
    RBF::useVertices = useVertices;
    RBF::nDer = nDer;

    PetscMalloc2(n, &(RBF::stencilWeights), 3*nDer, &(RBF::dxyz)) >> ablate::checkError;

    for (PetscInt c = cellRange.start; c < cellRange.end; ++c) {
      RBF::stencilWeights[c] = nullptr;
    }

    RestoreRange(cellRange);

    // Store the derivatives
    for (PetscInt i = 0; i < nDer; ++i) {
      RBF::dxyz[i*3 + 0] = dx[i];
      RBF::dxyz[i*3 + 1] = dy[i];
      RBF::dxyz[i*3 + 2] = dz[i];
    }

  }


}

// Set derivatives, defaulting to using vertices
void RBF::SetDerivatives(PetscInt nDer, PetscInt dx[], PetscInt dy[], PetscInt dz[]){
  RBF::SetDerivatives(nDer, dx, dy, dz, PETSC_TRUE);
}




// Compute the RBF weights at the cell center of p using a cell-list
// c - The center cell in cellRange ordering
void RBF::SetupDerivativeStencils(PetscInt c) {

  Mat       A = RBF::RBFMatrix[c], B = nullptr;
  PetscReal *x = nullptr;   // Shifted cell-centers of the neigbor list
  // This also computes nStencil[c]
  if (A == nullptr) {
    RBF::Matrix(c, &x, &A);
  }

  PetscInt    dim = subDomain->GetDimensions();
  PetscInt    nCells = RBF::nStencil[c];
  PetscInt    matSize = nCells + RBF::nPoly;
  PetscInt    nDer = RBF::nDer;
  PetscInt    *dxyz = RBF::dxyz;
  const PetscReal param = RBF::rbfParam;
  PetscInt    i, j;
  PetscInt    px, py, pz, p1 = RBF::polyOrder+1;
  PetscScalar *vals = nullptr;

  //Create the RHS
  MatCreateSeqDense(PETSC_COMM_SELF, matSize, nDer, NULL, &B) >> ablate::checkError;
  PetscObjectSetName((PetscObject)B,"ablate::radialBasis::RBF::rhs") >> ablate::checkError;
  MatZeroEntries(B) >> ablate::checkError;
  MatDenseGetArrayWrite(B, &vals) >> ablate::checkError;

  // Derivatives of the RBF
  for (i = 0; i < nCells; ++i) {
    for (j = 0; j < nDer; ++j) {
      vals[i + j*matSize] = RBFDer(dim, &x[i*dim], dxyz[j*3 + 0], dxyz[j*3 + 1], dxyz[j*3 + 2], param);
    }
  }

  // Derivatives of the augmented polynomials
  if (dim == 2) {
    for (j = 0; j < nDer; ++j) {
      i = nCells;
      for (py = 0; py < p1; ++py) {
        for (px = 0; px < p1-py; ++px ){
          if(dxyz[j*3 + 0] == px && dxyz[j*3 + 1] == py) {
            vals[i + j*matSize] = fac[px]*fac[py];
          }
          ++i;
        }
      }

    }
  } else {
    for (j = 0; j < nDer; ++j) {
      i = nCells;
      for (pz = 0; pz < p1; ++pz) {
        for (py = 0; py < p1-pz; ++py) {
          for (px = 0; px < p1-py-pz; ++px ){
            if(dxyz[j*3 + 0] == px && dxyz[j*3 + 1] == py && dxyz[j*3 + 2] == pz) {
              vals[i + j*matSize] = fac[px]*fac[py]*fac[pz];
            }
            ++i;
          }
        }
      }
    }
  }

  MatDenseRestoreArrayWrite(B, &vals) >> ablate::checkError;

  MatViewFromOptions(B,NULL,"-ablate::radialBasis::RBF::rhs_view") >> ablate::checkError;

  MatMatSolve(A, B, B) >> ablate::checkError;

  MatViewFromOptions(B,NULL,"-ablate::radialBasis::RBF::sol_view") >> ablate::checkError;

  // Now populate the output
  PetscMalloc1(nDer*nCells, &(RBF::stencilWeights[c])) >> ablate::checkError;
  PetscReal *wt = RBF::stencilWeights[c];
  MatDenseGetArrayWrite(B, &vals) >> ablate::checkError;
  for (i = 0; i < nCells; ++i) {
    for (j = 0; j < nDer; ++j) {
      wt[i*nDer + j] = vals[i + j*matSize];
    }
  }
  MatDenseGetArrayWrite(B, &vals) >> ablate::checkError;

  if (RBF::hasInterpolation) {
    RBF::RBFMatrix[c] = A;
    RBF::stencilXLocs[c] = x;
  }
  else {
    MatDestroy(&A) >> ablate::checkError;
    PetscFree(x) >> ablate::checkError;
  }
  MatDestroy(&B) >> ablate::checkError;

}



void RBF::SetupDerivativeStencils() {
  PetscInt c;
  ablate::solver::Range cellRange;

  GetCellRange(cellRange);

  for (c = cellRange.start; c < cellRange.end; ++c) {
    RBF::SetupDerivativeStencils(c);
  }

  RestoreRange(cellRange);

}



// Return the requested derivative
// field - The field to take the derivative of
// c - The location in ablate::solver::Range
// dx, dy, dz - The derivatives
PetscReal RBF::EvalDer(ablate::domain::Field field, PetscInt c, PetscInt dx, PetscInt dy, PetscInt dz){

  PetscInt  derID = -1, nDer = RBF::nDer, *dxyz = RBF::dxyz;
  // Search for the particular index. Probably need to do something different in the future to avoid re-doing the same calculation many times
  while ((dxyz[++derID*3 + 0] != dx || dxyz[derID*3 + 1] != dy || dxyz[derID*3 + 2] != dz) && derID<nDer){ }

  if (derID==nDer) {
    throw std::invalid_argument("Derivative of (" + std::to_string(dx) + ", " + std::to_string(dy) + ", " + std::to_string(dz) + ") is not setup.");
  }

  // If the stencil hasn't been setup yet do so
  if (RBF::nStencil[c] < 1) {
    RBF::SetupDerivativeStencils(c);
  }

  PetscReal         *wt = RBF::stencilWeights[c];
  PetscInt          nStencil = RBF::nStencil[c], *lst = RBF::stencilList[c];
  Vec               vec = subDomain->GetVec(field);
  DM                dm  = subDomain->GetFieldDM(field);
  PetscScalar       val = 0.0, f;
  const PetscScalar *array;
  VecGetArrayRead(vec, &array) >> checkError;

  for (PetscInt i = 0; i < nStencil; ++i) {
    DMPlexPointLocalFieldRead(dm, lst[i], field.id, array, &f) >> checkError;
    val += wt[i*nDer + derID]*f;
  }

  VecRestoreArrayRead(vec, &array);

  return val;
}



/************ End Derivative Code **********************/




/************ Begin Interpolation Code **********************/

void RBF::SetInterpolation(PetscBool hasInterpolation) {
  RBF::hasInterpolation = hasInterpolation;
}



// Return the interpolation of a field at a given location
// field - The field to interpolate
// c - The location in ablate::solver::Range
// xEval - The location to interpolate at
PetscReal RBF::Interpolate(ablate::domain::Field field, PetscInt c, PetscReal xEval[3]) {

  Mat         A = RBF::RBFMatrix[c];
  Vec         weights, rhs;
  PetscInt    nCells, *lst;
  PetscInt    i;
  PetscScalar *vals;
  const PetscScalar *fvals;
  PetscReal   *x = RBF::stencilXLocs[c], x0[3];
  Vec         f = subDomain->GetVec(field);
  DM          dm  = subDomain->GetFieldDM(field);

  if (A == nullptr) {
    RBF::Matrix(c, &x, &A);
    RBF::RBFMatrix[c] = A;
    RBF::stencilXLocs[c] = x;
  }

  nCells = RBF::nStencil[c];
  lst = RBF::stencilList[c];

  MatCreateVecs(A, &weights, &rhs) >> ablate::checkError;
  VecZeroEntries(weights) >> ablate::checkError;
  VecZeroEntries(rhs) >> ablate::checkError;

  // The function values
  VecGetArrayRead(f, &fvals) >> ablate::checkError;
  VecGetArray(rhs, &vals) >> ablate::checkError;
  for (i = 0; i < nCells; ++i) {
    DMPlexPointLocalFieldRead(dm, lst[i], field.id, fvals, &vals[i]) >> checkError;
  }
  VecRestoreArrayRead(f, &fvals) >> ablate::checkError;
  VecRestoreArray(rhs, &vals) >> ablate::checkError;

  MatSolve(A, rhs, weights) >> ablate::checkError;

  VecDestroy(&rhs) >> ablate::checkError;

  // Now do the actual interpolation

  // Get the cell center
  DMPlexComputeCellGeometryFVM(dm, c, NULL, x0, NULL) >> ablate::checkError;

  PetscInt p1 = RBF::polyOrder + 1, dim = subDomain->GetDimensions();
  PetscInt px, py, pz;
  PetscReal *xp;
  const PetscReal param = RBF::rbfParam;
  PetscMalloc1(dim*p1, &xp) >> ablate::checkError;

  for (PetscInt d = 0; d < dim; ++d) {
    x0[d] = xEval[d] - x0[d]; // Shifted center

    // precompute powers
    xp[d*p1 + 0] = 1.0;
    for (px = 1; px < p1; ++px) {
      xp[d*p1 + px] = xp[d*p1 + (px-1)]*x0[d];
    }
  }


  PetscReal   interpVal = 0.0;
  VecGetArray(weights, &vals) >> ablate::checkError;
  for (i = 0; i < nCells; ++i) {
    interpVal += vals[i]*RBFVal(dim, x0, &x[i*dim], param);
  }

  // Augmented polynomial contributions
  if (dim == 2) {
    for (py = 0; py < p1; ++py) {
      for (px = 0; px < p1-py; ++px) {
        interpVal += vals[i++] * xp[0*p1 + px] * xp[1*p1 + py];
      }
    }
  } else {
    for (pz = 0; pz < p1; ++pz) {
      for (py = 0; py < p1-pz; ++py) {
        for (px = 0; px < p1-py-pz; ++px ){
          interpVal += vals[i++] * xp[0*p1 + px] * xp[1*p1 + py] * xp[2*p1 + pz];
        }
      }
    }
  }
  VecRestoreArray(weights, &vals) >> ablate::checkError;
  VecDestroy(&weights) >> ablate::checkError;
  PetscFree(xp) >> ablate::checkError;

  return interpVal;

}



/************ End Interpolation Code **********************/


/************ Constructor, Setup, Registration, and Initialization Code **********************/

RBF::RBF(std::string solverId, std::shared_ptr<ablate::domain::Region> region, std::shared_ptr<ablate::parameters::Parameters> options,
  ablate::radialBasisV2::RBFType rbfType,
  PetscInt polyOrder,
  PetscReal rbfParam) :
    ablate::solver::Solver(solverId, region, options),
    rbfType(rbfType),
    polyOrder(polyOrder == 0 ? __RBF_DEFAULT_POLYORDER : polyOrder),
    rbfParam(rbfParam == 0 ? __RBF_DEFAULT_PARAM : rbfParam) { }

RBF::~RBF() {}


// This is done once
void RBF::Setup() {

  // Set the value and derivative functions
  switch (rbfType) {
    case RBFType::ga:
      RBFVal = &gaVal;
      RBFDer = &gaDer;
      break;
    case RBFType::imq:
      RBFVal = &imqVal;
      RBFDer = &imqDer;
      break;
    case RBFType::mq:
      RBFVal = &mqVal;
      RBFDer = &mqDer;
      break;
    case RBFType::phs:
      RBFVal = &phsVal;
      RBFDer = &phsDer;
      break;
    default:
      throw std::runtime_error("ablate::radialBasisV2::RBF has been passed an unknown type.");
  }

  const PetscInt dim = subDomain->GetDimensions();


  // The number of polynomial values is (p+2)(p+1)/2 in 2D and (p+3)(p+2)(p+1)/6 in 3D
  PetscInt p = RBF::polyOrder;
  if (dim == 2) {
    RBF::nPoly = (p+2)*(p+1)/2;
  } else {
    RBF::nPoly = (p+3)*(p+2)*(p+1)/6;
  }

  // Set the minimum number of cells to get compute the RBF matrix
  RBF::minNumberCells = (PetscInt)floor(2*(RBF::nPoly));

  // Now setup the derivatives required for curvature/normal calculations. This should probably move over to user-option
  PetscInt nDer = 0;
  PetscInt dx[10], dy[10], dz[10];


  nDer = ( dim == 2 ) ? 5 : 10;
  PetscInt i = 0;
  dx[i] = 1; dy[i] = 0; dz[i++] = 0;
  dx[i] = 0; dy[i] = 1; dz[i++] = 0;
  dx[i] = 2; dy[i] = 0; dz[i++] = 0;
  dx[i] = 0; dy[i] = 2; dz[i++] = 0;
  dx[i] = 1; dy[i] = 1; dz[i++] = 0;
  if( dim == 3) {
    dx[i] = 0; dy[i] = 0; dz[i++] = 1;
    dx[i] = 0; dy[i] = 0; dz[i++] = 2;
    dx[i] = 1; dy[i] = 0; dz[i++] = 1;
    dx[i] = 0; dy[i] = 1; dz[i++] = 1;
    dx[i] = 1; dy[i] = 1; dz[i++] = 1;
  }
  SetDerivatives(nDer, dx, dy, dz);

  // Let the RBF know that there will also be interpolation. This should probably move over to user-option
  SetInterpolation(PETSC_TRUE);


}

void RBF::Initialize() {

  // If this is called due to a grid change then release the old memory
  for (PetscInt c = 0; c < RBF::nCells; ++c ){
    PetscFree(RBF::stencilList[c]);
    if(RBF::RBFMatrix[c]) MatDestroy(&(RBF::RBFMatrix[c]));
    PetscFree(RBF::stencilXLocs[c]);
  }
  PetscFree4(RBF::nStencil, RBF::stencilList, RBF::RBFMatrix, RBF::stencilXLocs) >> ablate::checkError;


  ablate::solver::Range cellRange;
  GetCellRange(cellRange);

  // Both interpolation and derivatives need the list of points
  RBF::nCells = cellRange.end - cellRange.start;
  PetscMalloc4(RBF::nCells, &(RBF::nStencil), RBF::nCells, &(RBF::stencilList), RBF::nCells, &(RBF::RBFMatrix), RBF::nCells, &(RBF::stencilXLocs)) >> ablate::checkError;

  for (PetscInt c = cellRange.start; c < cellRange.end; ++c) {
    RBF::nStencil[c] = -1;
    RBF::stencilList[c] = nullptr;

    RBF::RBFMatrix[c] = nullptr;
    RBF::stencilXLocs[c] = nullptr;
  }

  RestoreRange(cellRange);
}

void RBF::Register(std::shared_ptr<ablate::domain::SubDomain> subDomain) { ablate::solver::Solver::Register(subDomain); }


//void ablate::radiation::VolumeRadiation::Setup() {
//    ablate::solver::CellSolver::Setup();
//    ablate::radiation::Radiation::Setup();
//}

//void ablate::radiation::VolumeRadiation::Register(std::shared_ptr<ablate::domain::SubDomain> subDomain) {
//    ablate::solver::Solver::Register(subDomain);
//    ablate::radiation::Radiation::Register(subDomain);
//}

//void ablate::radiation::VolumeRadiation::Initialize() {
//    solver::Range cellRange;
//    GetCellRange(cellRange);  //!< Gets the cell range that should be applied to the radiation solver

//    ablate::radiation::Radiation::Initialize(cellRange);  //!< Get the range of cells that the solver occupies in order for the radiation solver to give energy to the finite volume

//    RestoreRange(cellRange);
//}

#include "registrar.hpp"
REGISTER(ablate::solver::Solver, ablate::radialBasisV2::RBF, "Radial Basis Function",
         ARG(std::string, "id", "The name of the RBF."),
         OPT(ablate::domain::Region, "region", "The region to apply this solver.  Default is entire domain."),
         OPT(ablate::parameters::Parameters, "options", "The options passed to PETSC."),
         OPT(EnumWrapper<ablate::radialBasisV2::RBFType>, "rbfType", "Type of RBF. Default is MQ."),
         OPT(PetscInt, "polyOrder", "Order of the augmenting RBF polynomial. Must be >= 1. Default is 4."),
         OPT(PetscReal, "rbfParam", "Parameter required for the particular RBF. Default is 0.1.")
         );

