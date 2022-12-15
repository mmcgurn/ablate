#include "levelSetSolver.hpp"
#include "LS-VOF.hpp"
////#include <petsc/private/dmpleximpl.h>
////#include "utilities/mpiError.hpp"
////#include "utilities/petscError.hpp"

using namespace ablate::levelSet;



////ablate::levelSet::LevelSetSolver::LevelSetSolver(std::string solverId, std::shared_ptr<domain::Region> region, std::shared_ptr<parameters::Parameters> options)
////    : Solver(std::move(solverId), std::move(region), std::move(options)) {}

LevelSetSolver::LevelSetSolver(std::string solverId, std::shared_ptr<ablate::domain::Region> region, std::shared_ptr<ablate::parameters::Parameters> options) : Solver(solverId, region, options) { }


// This is done once
void LevelSetSolver::Setup() {


// Make sure that the level set field has been created in the YAML file.
  if (!(subDomain->ContainsField(LevelSetFields::LEVELSET_FIELD))) {
    throw std::runtime_error("ablate::levelSet::LevelSetSolver expects a level set field to be defined.");
  }
  if (!(subDomain->ContainsField(LevelSetFields::CURVATURE_FIELD))) {
    throw std::runtime_error("ablate::levelSet::LevelSetSolver expects a curvature field to be defined.");
  }
  if (!(subDomain->ContainsField(LevelSetFields::NORMAL_FIELD))) {
    throw std::runtime_error("ablate::levelSet::LevelSetSolver expects a normal field to be defined.");
  }

//  // Make sure that the RBF is setup.
//  subDomain->SetupRBF(subDomain);

  LevelSetSolver::lsField = &(subDomain->GetField(LevelSetFields::LEVELSET_FIELD));
  LevelSetSolver::curvField = &(subDomain->GetField(LevelSetFields::CURVATURE_FIELD));
  LevelSetSolver::normalField = &(subDomain->GetField(LevelSetFields::NORMAL_FIELD));

//  // Save the dimension
//  LevelSetSolver::dim = subDomain->GetDimensions();



printf("All Done!\n");
PetscFinalize();


}

void LevelSetSolver::Initialize() {

}

///*************   Begin Curvature and Normal Vector functions ******************/

void LevelSetSolver::Normal2D(PetscInt c, PetscScalar *n) {

  PetscReal   cx = 0.0, cy = 0.0, g = 0.0;
  std::shared_ptr<ablate::domain::RBF> rbf = subDomain->GetRBF();

  cx = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 1, 0, 0);
  cy = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 1, 0);
  g = PetscSqrtReal(cx*cx + cy*cy);

  n[0] = cx/g;
  n[1] = cy/g;


}

void LevelSetSolver::Normal3D(PetscInt c, PetscReal *n) {

  PetscReal   cx = 0.0, cy = 0.0, cz = 0.0, g = 0.0;
  std::shared_ptr<ablate::domain::RBF> rbf = subDomain->GetRBF();

  cx = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 1, 0, 0);
  cy = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 1, 0);
  cz = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 0, 1);
  g = sqrt(cx*cx + cy*cy + cz*cz);

  n[0] = cx/g;
  n[1] = cy/g;
  n[2] = cz/g;
}

PetscReal LevelSetSolver::Curvature2D(PetscInt c) {

  PetscReal k = 0.0;
  PetscReal cx, cy, cxx, cyy, cxy;
  std::shared_ptr<ablate::domain::RBF> rbf = subDomain->GetRBF();

  cx = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 1, 0, 0);
  cy = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 1, 0);
  cxx = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 2, 0, 0);
  cyy = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 2, 0);
  cxy = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 1, 1, 0);

  k = (cxx*cy*cy + cyy*cx*cx - 2.0*cxy*cx*cy)/pow(cx*cx+cy*cy,1.5);

  return k;
}

PetscReal LevelSetSolver::Curvature3D(PetscInt c) {

  PetscReal k = 0.0;
  PetscReal cx, cy, cz;
  PetscReal cxx, cyy, czz;
  PetscReal cxy, cxz, cyz;
  std::shared_ptr<ablate::domain::RBF> rbf = subDomain->GetRBF();

  cx = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 1, 0, 0);
  cy = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 1, 0);
  cz = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 0, 1);
  cxx = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 2, 0, 0);
  cyy = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 2, 0);
  czz = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 0, 2);
  cxy = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 1, 1, 0);
  cxz = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 1, 0, 1);
  cyz = LevelSetSolver::rbf->EvalDer(LevelSetSolver::lsField, c, 0, 1, 1);

  k = (cxx*(cy*cy + cz*cz) + cyy*(cx*cx + cz*cz) + czz*(cx*cx + cy*cy) - 2.0*(cxy*cx*cy + cxz*cx*cz + cyz*cy*cz))/pow(cx*cx+cy*cy+cz*cz,1.5);

  return k;
}

// There has to be a better way of doing this so that the curvature/normal function points directly to either 2D or 3D during setup.
PetscReal LevelSetSolver::Curvature(PetscInt c) {
  if (subDomain->GetDimensions()==2) {
    return LevelSetSolver::Curvature2D(c);
  }
  else {
    return LevelSetSolver::Curvature3D(c);
  }
}

void LevelSetSolver::Normal(PetscInt c, PetscReal *n) {
  if (subDomain->GetDimensions()==2) {
    return LevelSetSolver::Normal2D(c, n);
  }
  else {
    return LevelSetSolver::Normal3D(c, n);
  }
}


void LevelSetSolver::ComputeAllNormal() {
//  PetscScalar *val;
//  PetscInt    cStart, cEnd, c;

//  DMPlexGetHeightStratum(LevelSetField::dm, 0, &cStart, &cEnd) >> ablate::checkError;       // Range of cells

//  VecGetArray(LevelSetField::normal, &val) >> ablate::checkError;
//  for (c = cStart; c < cEnd; ++c) {
//    LevelSetField::Normal(c, &val[(c - cStart)*dim]);
//  }
//  VecRestoreArray(LevelSetField::normal, &val) >> ablate::checkError;
////  VecGhostUpdateBegin(LevelSetField::normal, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
////  VecGhostUpdateEnd(LevelSetField::normal, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;

}

void LevelSetSolver::ComputeAllCurvature() {
//  PetscScalar *val;
//  PetscInt    cStart, cEnd, c;

//  DMPlexGetHeightStratum(LevelSetField::dm, 0, &cStart, &cEnd) >> ablate::checkError;       // Range of cells

//  VecGetArray(LevelSetField::curv, &val) >> ablate::checkError;
//  for (c = cStart; c < cEnd; ++c) {
//    val[c - cStart] = LevelSetField::Curvature(c);
//  }
//  VecRestoreArray(LevelSetField::curv, &val) >> ablate::checkError;
//  VecGhostUpdateBegin(LevelSetField::curv, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
//  VecGhostUpdateEnd(LevelSetField::curv, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
}


/*************   End Curvature and Normal Vector functions ******************/







// Returns the VOF for a given cell. Refer to "Quadrature rules for triangular and tetrahedral elements with generalized functions"
//  by Holdych, Noble, and Secor, Int. J. Numer. Meth. Engng 2008; 73:1310-1327.
void LevelSetSolver::VOF(const PetscInt p, PetscReal *vof, PetscReal *area, PetscReal *vol) {

  DMPolytopeType    ct;
  DM                dm = subDomain->GetDM();
  const PetscInt    dim = subDomain->GetDimensions();
  PetscInt          Nc, nVerts, i, j;
  PetscReal         x0[3] = {0.0, 0.0, 0.0};
  PetscReal         *c = NULL, *coords = NULL;
  PetscScalar       c0, n[3] = {0.0, 0.0, 0.0}, g;
  const PetscScalar *array;
  PetscBool         isDG;
  Vec               solVec = subDomain->GetSolutionVector();  // For level set
  Vec               auxVec = subDomain->GetAuxVector();       // For normal vector




  // The cell center
  DMPlexComputeCellGeometryFVM(dm, p, NULL, x0, NULL) >> ablate::checkError;


  // Level-set value at cell-center
  VecGetArrayRead(solVec, &array) >> ablate::checkError;
  DMPlexPointLocalFieldRead(dm, p, LevelSetSolver::lsField->id, array, &c0) >> checkError;
  VecRestoreArrayRead(solVec, &array) >> ablate::checkError;

  // Normal vector
  VecGetArrayRead(auxVec, &array) >> ablate::checkError;
  DMPlexPointLocalFieldRead(dm, p, LevelSetSolver::normalField->id, array, n) >> checkError;
  VecRestoreArrayRead(auxVec, &array) >> ablate::checkError;

  g = 0.0;
  for (i = 0; i < dim; ++i) g += PetscSqr(n[i]);
  g = PetscSqrtReal(g);
  for (i = 0; i < dim; ++i) n[i] /= g;




  // Coordinates of the cell vertices
  DMPlexGetCellCoordinates(dm, p, &isDG, &Nc, &array, &coords) >> ablate::checkError;

  // Number of vertices
  nVerts = Nc/dim;

  PetscMalloc1(nVerts, &c) >> ablate::checkError;

  // The level set value of each vertex. This assumes that the interface is a line/plane
  //    with the given unit normal.
  for (i = 0; i < nVerts; ++i) {
    c[i] = c0;
    for (j = 0; j < dim; ++j) {
      c[i] += n[j]*(coords[i*dim + j] - x0[j]);
    }
  }

  // Get the cell type and call appropriate VOF function
  DMPlexGetCellType(dm, p, &ct) >> ablate::checkError;
  switch (ct) {
    case DM_POLYTOPE_SEGMENT:
      throw std::invalid_argument("No element geometry for cell " + std::to_string(p) + " with type " + DMPolytopeTypes[ct]);
    case DM_POLYTOPE_TRIANGLE:
      VOF_2D_Tri(coords, c, vof, area, vol);
      break;
    case DM_POLYTOPE_QUADRILATERAL:
      VOF_2D_Quad(coords, c, vof, area, vol);
      break;
    case DM_POLYTOPE_TETRAHEDRON:
      VOF_3D_Tetra(coords, c, vof, area, vol);
      break;
    case DM_POLYTOPE_HEXAHEDRON:
      VOF_3D_Hex(coords, c, vof, area, vol);
      break;
    default:
      throw std::invalid_argument("No element geometry for cell " + std::to_string(p) + " with type " + DMPolytopeTypes[ct]);
  }

  DMPlexRestoreCellCoordinates(dm, p, &isDG, &Nc, &array, &coords) >> ablate::checkError;
  PetscFree(c) >> ablate::checkError;

}











//// Reinitialize a level set field to make it a signed distance function and to match a target VOF for each cell
//void LevelSetSolver::Reinitialize(TS ts, ablate::solver::Solver &solver) {
////    // Get the solution vec and dm
////    auto dm = solver.GetSubDomain().GetDM();
////    auto solVec = solver.GetSubDomain().GetSolutionVector();
////    auto auxDm = solver.GetSubDomain().GetAuxDM();
////    auto auxVec = solver.GetSubDomain().GetAuxVector();

////    // Get the array vector
////    PetscScalar *solutionArray;
////    VecGetArray(solVec, &solutionArray) >> checkError;
////    PetscScalar *auxArray;
////    VecGetArray(auxVec, &auxArray) >> checkError;

////    // March over each cell in this domain
////    solver::Range cellRange;
////    solver.GetCellRange(cellRange);
////    auto dim = solver.GetSubDomain().GetDimensions();

////    for (PetscInt c = cellRange.start; c < cellRange.end; ++c) {
////        PetscInt cell = cellRange.points ? cellRange.points[c] : c;


////        // Get the euler and density field
////        const PetscScalar *euler = nullptr;
////        DMPlexPointGlobalFieldRef(dm, cell, eulerFieldInfo.id, solutionArray, &euler) >> checkError;
////        PetscScalar *densityYi;
////        DMPlexPointGlobalFieldRef(dm, cell, densityYiFieldInfo.id, solutionArray, &densityYi) >> checkError;
////        PetscScalar *yi;
////        DMPlexPointLocalFieldRead(auxDm, cell, yiFieldInfo.id, auxArray, &yi) >> checkError;
////        PetscFVCellGeom *cellGeom;
////        DMPlexPointLocalRead(cellGeomDm, cell, cellGeomArray, &cellGeom) >> checkError;

////        // compute the mass fractions on the boundary
////        massFractionsFunction(dim, time, cellGeom->centroid, yiFieldInfo.numberComponents, yi, massFractionsContext);

////        // Only update if in the global vector
////        if (euler) {
////            // Get density
////            const PetscScalar density = euler[finiteVolume::CompressibleFlowFields::RHO];

////            for (PetscInt sp = 0; sp < densityYiFieldInfo.numberComponents; sp++) {
////                densityYi[sp] = yi[sp] * density;
////            }
////        }
////    }

////    // cleanup
////    VecRestoreArrayRead(cellGeomVec, &cellGeomArray) >> checkError;
////    VecRestoreArray(auxVec, &auxArray) >> checkError;
////    VecRestoreArray(solVec, &solutionArray) >> checkError;
////    solver.RestoreRange(cellRange);









////  PetscInt          c, cStart, cEnd;
////  DM                dm = LevelSetField::dm;
////  const PetscScalar *vofVal;
////  PetscScalar       *phiVal;
////  PetscReal         vof, faceArea, cellVolume;
////  Vec               newPhi;

////  VecDuplicate(LevelSetField::phi, &newPhi);

////  VecDuplicate(newPhi, &newPhi);

////  VecGetArrayRead(VOF, &vofVal) >> ablate::checkError;
////  VecGetArray(newPhi, &phiVal) >> ablate::checkError;

//////Take a look at boundarySolver/physics/sublimation.cpp lines 233-246 + 276

//////Use DMPlexPointGlobalFieldRead to get field values
//////Stuff like const auto &eulerFieldInfo = solver.GetSubDomain().GetField(finiteVolume::CompressibleFlowFields::EULER_FIELD); will return the field info in the DM.
//////Make the level set a solution variable in the ablate solver


////  DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd) >> ablate::checkError;
////  for (c = cStart; c < cEnd; ++c) {
////    LevelSetField::VOF(c, &vof, &faceArea, &cellVolume);

////  }

////  VecRestoreArray(newPhi, &phiVal);
////  VecRestoreArrayRead(VOF, &vofVal) >> ablate::checkError;
////  VecDestroy(&newPhi);


//}


////std::string LevelSetSolver::GetRBFType() {

////  switch (LevelSetSolver::rbfType) {
////    case ablate::RBF::RBFType::PHS:
////      return("phs");
////    case ablate::RBF::RBFType::MQ:
////      return("mq");
////    case ablate::RBF::RBFType::IMQ:
////      return("imq");
////    case ablate::RBF::RBFType::GA:
////      return("ga");
////    default:
////      return("unknown");
////  }
////}



////bool LevelSetField::HasInterface(const PetscInt p) {
////  bool              hasInterface = false;
////  PetscInt          nCells = 0, *cells = NULL;
////  PetscInt          i, cStart;
////  Vec               phi = LevelSetField::phi;
////  const PetscScalar *array;
////  PetscScalar       c0;
////  DM                dm = LevelSetField::dm;

////  DMPlexGetHeightStratum(dm, 0, &cStart, NULL) >> ablate::checkError;

////  DMPlexGetNeighborCells(dm, p, 1, -1, -1, PETSC_TRUE, &nCells, &cells) >> ablate::checkError;

////  VecGetArrayRead(phi, &array) >> ablate::checkError;
////  c0 = array[p - cStart];

////  i = 0;
////  while (i < nCells && !hasInterface) {
////    hasInterface = ((c0 * array[cells[i] - cStart])<=0.0);
////    ++i;
////  }

////  VecRestoreArrayRead(phi, &array) >> ablate::checkError;
////  PetscFree(cells) >> ablate::checkError;

////  return hasInterface;

////}

/////* Sphere */
////PetscReal LevelSetField::Sphere(PetscReal pos[], PetscReal center[], PetscReal radius) {
////  PetscReal shiftedPos[] = {pos[0] - center[0], pos[1] - center[1], pos[2] - center[2]};
////  PetscReal phi = PetscSqrtReal(PetscSqr(shiftedPos[0]) + PetscSqr(shiftedPos[1]) + PetscSqr(shiftedPos[2])) - radius;
////  return phi;
////}

/////* Ellipse */
////PetscReal LevelSetField::Ellipse(PetscReal pos[], PetscReal center[], PetscReal radius) {
////  PetscReal shiftedPos[] = {pos[0] - center[0], pos[1] - center[1], pos[2] - center[2]};
////  PetscReal phi = PetscSqr(shiftedPos[0]/0.5) + PetscSqr(shiftedPos[1]/1.25) + PetscSqr(shiftedPos[2]) - radius;
////  return phi;
////}


/////* Star */
////PetscReal LevelSetField::Star(PetscReal pos[], PetscReal center[]) {
////  PetscReal shiftedPos[] = {pos[0] - center[0], pos[1] - center[1], pos[2] - center[2]};
////  PetscReal phi = 400.0*shiftedPos[0]*shiftedPos[0]*shiftedPos[1]*shiftedPos[1]-(1.0-0.5*shiftedPos[0]*shiftedPos[0]-0.5*shiftedPos[1]*shiftedPos[1]);
////  return phi;
////}

//PetscReal LevelSetSolver::Interpolate(PetscScalar xyz[3]) {
//  std::shared_ptr<ablate::domain::rbf::RBF>  rbf = LevelSetSolver::rbf;
//  DMInterpolationInfo   ctx;
//  DM                    dm = rbf->GetDM();
//  PetscInt              c = -1;
//  Vec                   phi = LevelSetField::phi;
//  PetscReal             val;

//  DMInterpolationCreate(PETSC_COMM_WORLD, &ctx) >> ablate::checkError;
//  DMInterpolationSetDim(ctx, LevelSetField::dim) >> ablate::checkError;
//  DMInterpolationAddPoints(ctx, 1, xyz) >> ablate::checkError;
//  DMInterpolationSetUp(ctx, dm, PETSC_FALSE, PETSC_FALSE) >> ablate::checkError;
//  c = ctx->cells[0];
//  DMInterpolationDestroy(&ctx) >> ablate::checkError;


//  PetscReal RBF::Interpolate(const ablate::domain::Field *field, PetscInt c, PetscReal xEval[3]) {

//  val = rbf->Interpolate(phi, c, xyz);

//  return val;
//}

//PetscReal LevelSetSolver::Interpolate(const PetscReal x, const double y, const double z) {

//  PetscReal xyz[3] = {x, y, z};
//  PetscReal val = LevelSetSolver::Interpolate(xyz);

//  return val;
//}



////void LevelSetField::Advect(Vec velocity, const PetscReal dt) {

////  Vec               phi = LevelSetField::phi, nextPhi = nullptr;
////  DM                dm = LevelSetField::dm;
////  PetscInt          dim = LevelSetField::dim;
////  PetscInt          cStart, cEnd, c, cShift;
////  PetscScalar       *newVal;
////  const PetscScalar *vel;
////  PetscReal         pos[3] = {0.0, 0.0, 0.0};


////  VecDuplicate(phi, &nextPhi);

////  DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd) >> ablate::checkError;       // Range of cells

////  VecGetArray(nextPhi, &newVal) >> ablate::checkError;
////  VecGetArrayRead(velocity, &vel) >> ablate::checkError;
////  for (c = cStart; c < cEnd; ++c) {
////    cShift = c - cStart;
////    // Cell center
////    DMPlexComputeCellGeometryFVM(dm, c, NULL, pos, NULL) >> ablate::checkError;

////    // Step backward
////    for (PetscInt d = 0; d < dim; ++d) {
////      pos[d] -= dt*vel[cShift*dim + d];
////    }

////    newVal[cShift] = LevelSetField::Interpolate(pos);
////  }
////  VecRestoreArrayRead(velocity, &vel) >> ablate::checkError;
////  VecRestoreArray(nextPhi, &newVal) >> ablate::checkError;

////  VecCopy(nextPhi, phi) >> ablate::checkError;
////  VecDestroy(&nextPhi) >> ablate::checkError;

////  VecGhostUpdateBegin(phi, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
////  VecGhostUpdateEnd(phi, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;


////  VecDestroy(&nextPhi);



////}



////// Reinitialize a level set field to make it a signed distance function and to match a target VOF for each cell
////void LevelSetField::Reinitialize(Vec VOF) {
////  PetscInt          c, cStart, cEnd;
////  DM                dm = LevelSetField::dm;
////  const PetscScalar *vofVal;
////  PetscScalar       *phiVal;
////  PetscReal         vof, faceArea, cellVolume;
////  Vec               newPhi;

////  VecDuplicate(LevelSetField::phi, &newPhi);

////  VecDuplicate(newPhi, &newPhi);

////  VecGetArrayRead(VOF, &vofVal) >> ablate::checkError;
////  VecGetArray(newPhi, &phiVal) >> ablate::checkError;

//////Take a look at boundarySolver/physics/sublimation.cpp lines 233-246 + 276

//////Use DMPlexPointGlobalFieldRead to get field values
//////Stuff like const auto &eulerFieldInfo = solver.GetSubDomain().GetField(finiteVolume::CompressibleFlowFields::EULER_FIELD); will return the field info in the DM.
//////Make the level set a solution variable in the ablate solver


////  DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd) >> ablate::checkError;
////  for (c = cStart; c < cEnd; ++c) {
////    LevelSetField::VOF(c, &vof, &faceArea, &cellVolume);

////  }

////  VecRestoreArray(newPhi, &phiVal);
////  VecRestoreArrayRead(VOF, &vofVal) >> ablate::checkError;
////  VecDestroy(&newPhi);


////}


#include "registrar.hpp"
REGISTER(ablate::solver::Solver, ablate::levelSet::LevelSetSolver, "level set solver",
         ARG(std::string, "id", "the name of the level set solver"),
         OPT(ablate::domain::Region, "region", "the region to apply this solver.  Default is entire domain"),
         OPT(ablate::parameters::Parameters, "options", "the options passed to PETSC for the solver"));

