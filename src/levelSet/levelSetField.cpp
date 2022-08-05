#include "levelSetField.hpp"


using namespace ablate::levelSet;

void LevelSetField::Normal2D(PetscInt c, PetscScalar *n) {

  PetscReal             cx, cy, g;
  Vec                   phi = LevelSetField::phi;
  std::shared_ptr<RBF>  rbf = LevelSetField::rbf;

  cx = rbf->EvalDer(phi, c, 1, 0, 0);
  cy = rbf->EvalDer(phi, c, 0, 1, 0);
  g = PetscSqrtReal(cx*cx + cy*cy);

  n[0] = cx/g;
  n[1] = cy/g;


}

void LevelSetField::Normal3D(PetscInt c, PetscReal *n) {

  PetscReal             cx, cy, cz, g;
  Vec                   phi = LevelSetField::phi;
  std::shared_ptr<RBF>  rbf = LevelSetField::rbf;

  cx = rbf->EvalDer(phi, c, 1, 0, 0);
  cy = rbf->EvalDer(phi, c, 0, 1, 0);
  cz = rbf->EvalDer(phi, c, 0, 0, 1);
  g = sqrt(cx*cx + cy*cy + cz*cz);

  n[0] = cx/g;
  n[1] = cy/g;
  n[2] = cz/g;
}

PetscReal LevelSetField::Curvature2D(PetscInt c) {

  PetscReal             k = 0.0;
  PetscReal             cx, cy, cxx, cyy, cxy;
  Vec                   phi = LevelSetField::phi;
  std::shared_ptr<RBF>  rbf = LevelSetField::rbf;

  cx = rbf->EvalDer(phi, c, 1, 0, 0);
  cy = rbf->EvalDer(phi, c, 0, 1, 0);
  cxx = rbf->EvalDer(phi, c, 2, 0, 0);
  cyy = rbf->EvalDer(phi, c, 0, 2, 0);
  cxy = rbf->EvalDer(phi, c, 1, 1, 0);

  k = (cxx*cy*cy + cyy*cx*cx - 2.0*cxy*cx*cy)/pow(cx*cx+cy*cy,1.5);

  return k;
}

PetscReal LevelSetField::Curvature3D(PetscInt c) {

  PetscReal             k = 0.0;
  PetscReal             cx, cy, cz;
  PetscReal             cxx, cyy, czz;
  PetscReal             cxy, cxz, cyz;
  Vec                   phi = LevelSetField::phi;
  std::shared_ptr<RBF>  rbf = LevelSetField::rbf;

  cx = rbf->EvalDer(phi, c, 1, 0, 0);
  cy = rbf->EvalDer(phi, c, 0, 1, 0);
  cz = rbf->EvalDer(phi, c, 0, 0, 1);
  cxx = rbf->EvalDer(phi, c, 2, 0, 0);
  cyy = rbf->EvalDer(phi, c, 0, 2, 0);
  czz = rbf->EvalDer(phi, c, 0, 0, 2);
  cxy = rbf->EvalDer(phi, c, 1, 1, 0);
  cxz = rbf->EvalDer(phi, c, 1, 0, 1);
  cyz = rbf->EvalDer(phi, c, 0, 1, 1);

  k = (cxx*(cy*cy + cz*cz) + cyy*(cx*cx + cz*cz) + czz*(cx*cx + cy*cy) - 2.0*(cxy*cx*cy + cxz*cx*cz + cyz*cy*cz))/pow(cx*cx+cy*cy+cz*cz,1.5);

  return k;
}

// There has to be a better way of doing this so that the curvature/normal function points directly to either 2D or 3D during setup.
PetscReal LevelSetField::Curvature(PetscInt c) {
  if (LevelSetField::dim==2) {
    return LevelSetField::Curvature2D(c);
  }
  else {
    return LevelSetField::Curvature3D(c);
  }
}

void LevelSetField::Normal(PetscInt c, PetscReal *n) {
  if (LevelSetField::dim==2) {
    return LevelSetField::Normal2D(c, n);
  }
  else {
    return LevelSetField::Normal3D(c, n);
  }
}


LevelSetField::LevelSetField(std::shared_ptr<domain::Region> region) : region(region) {}




LevelSetField::LevelSetField(std::shared_ptr<RBF> rbf, LevelSetField::levelSetShape shape) {

  LevelSetField::rbf = rbf;
  LevelSetField::dm = rbf->GetDM();

  PetscInt            d, dim;
  PetscInt            cStart, cEnd, c;
  PetscInt            nSet = 3;
  PetscScalar         *val;
  PetscReal           lo[] = {0.0, 0.0, 0.0}, hi[] = {0.0, 0.0, 0.0}, centroid[] = {0.0, 0.0, 0.0}, pos[] = {0.0, 0.0, 0.0};
  PetscReal           radius = 1.0;



  DMGetDimension(dm, &dim) >> ablate::checkError;
  LevelSetField::dim = dim;




  // Create the vectors
  DMCreateGlobalVector(dm, &(LevelSetField::phi)) >> ablate::checkError;
  DMCreateGlobalVector(dm, &(LevelSetField::curv)) >> ablate::checkError;

  PetscInt lsz, gsz;
  VecGetLocalSize(phi, &lsz) >> ablate::checkError;
  VecGetSize(phi, &gsz) >> ablate::checkError;
  VecCreateMPI(PETSC_COMM_WORLD, dim*lsz, dim*gsz, &(LevelSetField::normal)) >> ablate::checkError;
  VecSetBlockSize(LevelSetField::normal, dim) >> ablate::checkError;

  DMGetBoundingBox(dm, lo, hi) >> ablate::checkError;
  for (d = 0; d < dim; ++d) {
    centroid[d] = 0.5*(lo[d] + hi[d]);
  }

  PetscOptionsGetReal(NULL, NULL, "-radius", &(radius), NULL) >> ablate::checkError;
  PetscOptionsGetRealArray(NULL, NULL, "-centroid", centroid, &nSet, NULL) >> ablate::checkError;


  DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd) >> ablate::checkError;       // Range of cells
  VecGetArray(LevelSetField::phi, &val) >> ablate::checkError;
  // Set the initial shape.
  switch (shape) {
    case LevelSetField::levelSetShape::SPHERE:
     for (c = cStart; c < cEnd; ++c) {
       DMPlexComputeCellGeometryFVM(dm, c, NULL, pos, NULL) >> ablate::checkError;
       val[c - cStart] = LevelSetField::Sphere(pos, centroid, radius);
      }

      break;
    case LevelSetField::levelSetShape::ELLIPSE:
      for (c = cStart; c < cEnd; ++c) {
        DMPlexComputeCellGeometryFVM(dm, c, NULL, pos, NULL) >> ablate::checkError;
        val[c - cStart] = LevelSetField::Ellipse(pos, centroid, radius);
      }

      break;
    case LevelSetField::levelSetShape::STAR:
      for (c = cStart; c < cEnd; ++c) {
        DMPlexComputeCellGeometryFVM(dm, c, NULL, pos, NULL) >> ablate::checkError;
        val[c - cStart] = LevelSetField::Star(pos, centroid);
      }

      break;
    default:
      throw std::invalid_argument("Unknown level set shape shape");
  }

  VecRestoreArray(LevelSetField::phi, &val) >> ablate::checkError;

  VecGhostUpdateBegin(LevelSetField::phi, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
  VecGhostUpdateEnd(LevelSetField::phi, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;

  // Now setup the derivatives and set the curvature/normal calculations
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
  rbf->SetDerivatives(nDer, dx, dy, dz);

}

LevelSetField::~LevelSetField() {
  VecDestroy(&(LevelSetField::phi));
  VecDestroy(&(LevelSetField::normal));
  VecDestroy(&(LevelSetField::curv));


}

void LevelSetField::ComputeAllNormal() {
  PetscScalar *val;
  PetscInt    cStart, cEnd, c;

  DMPlexGetHeightStratum(LevelSetField::dm, 0, &cStart, &cEnd) >> ablate::checkError;       // Range of cells

  VecGetArray(LevelSetField::normal, &val) >> ablate::checkError;
  for (c = cStart; c < cEnd; ++c) {
    LevelSetField::Normal(c, &val[(c - cStart)*dim]);
  }
  VecRestoreArray(LevelSetField::normal, &val) >> ablate::checkError;
//  VecGhostUpdateBegin(LevelSetField::normal, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
//  VecGhostUpdateEnd(LevelSetField::normal, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;

}

void LevelSetField::ComputeAllCurvature() {
  PetscScalar *val;
  PetscInt    cStart, cEnd, c;

  DMPlexGetHeightStratum(LevelSetField::dm, 0, &cStart, &cEnd) >> ablate::checkError;       // Range of cells

  VecGetArray(LevelSetField::curv, &val) >> ablate::checkError;
  for (c = cStart; c < cEnd; ++c) {
    val[c - cStart] = LevelSetField::Curvature(c);
  }
  VecRestoreArray(LevelSetField::curv, &val) >> ablate::checkError;
  VecGhostUpdateBegin(LevelSetField::curv, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
  VecGhostUpdateEnd(LevelSetField::curv, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
}

std::vector<std::shared_ptr<ablate::domain::FieldDescription>> LevelSetField::GetFields() {
    std::vector<std::shared_ptr<ablate::domain::FieldDescription>> levelSetField{
        std::make_shared<domain::FieldDescription>("level set field", "phi", domain::FieldDescription::ONECOMPONENT, domain::FieldLocation::SOL, domain::FieldType::FVM, region)};

  return levelSetField;
}

PetscReal LevelSetField::Interpolate(PetscScalar xyz[3]) {
  std::shared_ptr<RBF>  rbf = LevelSetField::rbf;
  DMInterpolationInfo   ctx;
  DM                    dm = rbf->GetDM();
  PetscInt              c = -1;
  Vec                   phi = LevelSetField::phi;
  PetscReal             val;

  DMInterpolationCreate(PETSC_COMM_WORLD, &ctx) >> ablate::checkError;
  DMInterpolationSetDim(ctx, LevelSetField::dim) >> ablate::checkError;
  DMInterpolationAddPoints(ctx, 1, xyz) >> ablate::checkError;
  DMInterpolationSetUp(ctx, dm, PETSC_FALSE, PETSC_TRUE) >> ablate::checkError;
  c = ctx->cells[0];
  DMInterpolationDestroy(&ctx) >> ablate::checkError;

  val = rbf->Interpolate(phi, c, xyz);

  return val;
}

PetscReal LevelSetField::Interpolate(const PetscReal x, const double y, const double z) {

  PetscReal xyz[3] = {x, y, z};
  PetscReal val = LevelSetField::Interpolate(xyz);

  return val;
}



void LevelSetField::Advect(Vec velocity, const PetscReal dt) {

  Vec               phi = LevelSetField::phi, nextPhi = nullptr;
  DM                dm = LevelSetField::dm;
  PetscInt          dim = LevelSetField::dim;
  PetscInt          cStart, cEnd, c, cShift;
  PetscScalar       *newVal;
  const PetscScalar *vel;
  PetscReal         pos[3] = {0.0, 0.0, 0.0};


  VecDuplicate(phi, &nextPhi);

  DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd) >> ablate::checkError;       // Range of cells

  VecGetArray(nextPhi, &newVal) >> ablate::checkError;
  VecGetArrayRead(velocity, &vel) >> ablate::checkError;
  for (c = cStart; c < cEnd; ++c) {
    cShift = c - cStart;
    // Cell center
    DMPlexComputeCellGeometryFVM(dm, c, NULL, pos, NULL) >> ablate::checkError;

    // Step backward
    for (PetscInt d = 0; d < dim; ++d) {
      pos[d] -= dt*vel[cShift*dim + d];
    }

    newVal[cShift] = LevelSetField::Interpolate(pos);
  }
  VecRestoreArrayRead(velocity, &vel) >> ablate::checkError;
  VecRestoreArray(nextPhi, &newVal) >> ablate::checkError;

  VecCopy(nextPhi, phi) >> ablate::checkError;
  VecDestroy(&nextPhi) >> ablate::checkError;

  VecGhostUpdateBegin(phi, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;
  VecGhostUpdateEnd(phi, INSERT_VALUES, SCATTER_FORWARD) >> ablate::checkError;


  VecDestroy(&nextPhi);



}

bool LevelSetField::HasInterface(PetscInt p) {
  bool              hasInterface = false;
  PetscInt          nCells = 0, *cells = NULL;
  PetscInt          i, cStart;
  Vec               phi = LevelSetField::phi;
  const PetscScalar *array;
  PetscScalar       c0;
  DM                dm = LevelSetField::dm;

  DMPlexGetHeightStratum(dm, 0, &cStart, NULL) >> ablate::checkError;

  DMPlexGetNeighborCells(dm, p, 1, -1, -1, PETSC_TRUE, &nCells, &cells) >> ablate::checkError;

  VecGetArrayRead(phi, &array) >> ablate::checkError;
  c0 = array[p - cStart];

  i = 0;
  while (i < nCells && !hasInterface) {
    hasInterface = ((c0 * array[cells[i] - cStart])<=0.0);
    ++i;
  }

  VecRestoreArrayRead(phi, &array) >> ablate::checkError;
  PetscFree(cells) >> ablate::checkError;

  return hasInterface;

}

/* Sphere */
PetscReal LevelSetField::Sphere(PetscReal pos[], PetscReal center[], PetscReal radius) {
  PetscReal shiftedPos[] = {pos[0] - center[0], pos[1] - center[1], pos[2] - center[2]};
  PetscReal phi = PetscSqrtReal(PetscSqr(shiftedPos[0]) + PetscSqr(shiftedPos[1]) + PetscSqr(shiftedPos[2])) - radius;
  return phi;
}

/* Ellipse */
PetscReal LevelSetField::Ellipse(PetscReal pos[], PetscReal center[], PetscReal radius) {
  PetscReal shiftedPos[] = {pos[0] - center[0], pos[1] - center[1], pos[2] - center[2]};
  PetscReal phi = PetscSqr(shiftedPos[0]/0.5) + PetscSqr(shiftedPos[1]/1.25) + PetscSqr(shiftedPos[2]) - radius;
  return phi;
}


/* Star */
PetscReal LevelSetField::Star(PetscReal pos[], PetscReal center[]) {
  PetscReal shiftedPos[] = {pos[0] - center[0], pos[1] - center[1], pos[2] - center[2]};
  PetscReal phi = 400.0*shiftedPos[0]*shiftedPos[0]*shiftedPos[1]*shiftedPos[1]-(1.0-0.5*shiftedPos[0]*shiftedPos[0]-0.5*shiftedPos[1]*shiftedPos[1]);
  return phi;
}




#include "registrar.hpp"
REGISTER(ablate::domain::FieldDescriptor, LevelSetField, "Level Set fields need for interface tracking",
         OPT(ablate::domain::Region, "region", "the region for the level set (defaults to entire domain)"));
