#include "dpm.h"
#include <functional>

using namespace Eigen;
using namespace std;

/******************************

        C O N S T R U C T O R S  &

                D E S T R U C T O R

*******************************/

// Main constructor
dpm::dpm(int n, int ndim, int seed) {
  // local variables
  int d, i;

  // print to console
  cout << "** Instantiating configobj2D object, NCELLS = " << n << ",  ndim = " << ndim << ", seed = " << seed << " ..." << endl;

  // main variables
  NCELLS = n;
  NDIM = ndim;
  NNN = 4;

  // set scalars to default values
  dt = 0.0;

  ka = 0.0;
  kl = 0.0;
  kb = 0.0;
  kc = 0.0;

  maxwellRelaxationTime = INFINITY;
  taus = INFINITY;
  B = 0.0;

  l1 = 0.0;
  l2 = 0.0;

  // default boundary variables
  L.resize(NDIM);
  pbc.resize(NDIM);
  for (d = 0; d < NDIM; d++) {
    L[d] = 1.0;
    pbc[d] = 1;
  }

  // preferred area for each cell
  a0.resize(NCELLS);
  cellU.resize(NCELLS);

  // macroscopic stress vector
  stress.resize(NDIM * (NDIM + 1) / 2);
  for (i = 0; i < NDIM * (NDIM + 1) / 2; i++)
    stress.at(i) = 0.0;

  // contact network vector
  cij.resize(NCELLS * (NCELLS - 1) / 2);
  for (i = 0; i < NCELLS * (NCELLS - 1) / 2; i++)
    cij.at(i) = 0;

  // initialize nearest neighbor info
  NBX = -1;

  // seed random number generator
  srand48(seed);
}

/*// destructor
dpm::~dpm() {
  // clear all private vectors
  // should update this soon
  L.clear();
  pbc.clear();
  a0.clear();
  l0.clear();
  t0.clear();
  nv.clear();
  szList.clear();
  im1.clear();
  ip1.clear();
  r.clear();
  x.clear();
  v.clear();
  F.clear();
  stress.clear();
  sb.clear();
  lb.clear();
  for (int i = 0; i < nn.size(); i++) {
    // cout << i << "\t" << nn.size() << '\n';
    nn.at(i).clear();
  }
  nn.clear();
  head.clear();
  last.clear();
  list.clear();

  if (posout.is_open())
    posout.close();
}*/

/******************************

        C E L L   S H A P E

        G E T T E R S

*******************************/

// get global vertex index gi given input cell index ci and vertex index vi
int dpm::gindex(int ci, int vi) {
  return szList[ci] + vi;
}

// get cell index ci and vertex index
void dpm::cindices(int& ci, int& vi, int gi) {
  for (int i = NCELLS - 1; i >= 0; i--) {
    if (gi >= szList[i]) {
      ci = i;
      vi = gi - szList[ci];
      break;
    }
  }
}

// get cell area
double dpm::area(int ci) {
  // local variables
  int vi, gi, gip1, nvtmp;
  double dx, dy, xi, yi, xip1, yip1, areaVal = 0.0;

  // initial position: vi = 0
  nvtmp = nv.at(ci);
  gi = gindex(ci, 0);
  xi = x[NDIM * gi];
  yi = x[NDIM * gi + 1];

  // loop over vertices of cell ci, get area by shoe-string method
  for (vi = 0; vi < nvtmp; vi++) {
    // next vertex
    gip1 = ip1[gi];
    gi++;

    // get positions (check minimum images)
    dx = x[NDIM * gip1] - xi;
    if (pbc[0])
      dx -= L[0] * round(dx / L[0]);
    xip1 = xi + dx;

    dy = x[NDIM * gip1 + 1] - yi;
    if (pbc[1])
      dy -= L[1] * round(dy / L[1]);
    yip1 = yi + dy;

    // increment area
    areaVal += xi * yip1 - xip1 * yi;

    // set next coordinates
    xi = xip1;
    yi = yip1;
  }
  areaVal *= 0.5;

  return abs(areaVal);
}

// get cell perimeter
double dpm::perimeter(int ci) {
  // local variables
  int vi, gi, gip1, nvtmp;
  double dx, dy, xi, yi, xip1, yip1, l, perimVal = 0.0;

  // initial position: vi = 0
  nvtmp = nv.at(ci);
  gi = gindex(ci, 0);
  xi = x[NDIM * gi];
  yi = x[NDIM * gi + 1];

  // loop over vertices of cell ci, get perimeter
  for (vi = 0; vi < nvtmp; vi++) {
    // next vertex
    gip1 = ip1[gi];
    gi++;

    // get positions (check minimum images)
    dx = x[NDIM * gip1] - xi;
    if (pbc[0])
      dx -= L[0] * round(dx / L[0]);
    xip1 = xi + dx;

    dy = x[NDIM * gip1 + 1] - yi;
    if (pbc[1])
      dy -= L[1] * round(dy / L[1]);
    yip1 = yi + dy;

    // compute segment length
    l = sqrt(dx * dx + dy * dy);

    // add to perimeter
    perimVal += l;

    // update coordinates
    xi = xip1;
    yi = yip1;
  }

  // return perimeter
  return perimVal;
}

// get cell center of mass position
void dpm::com2D(int ci, double& cx, double& cy) {
  // local variables
  int vi, gi, gip1, nvtmp;
  double dx, dy, xi, yi, xip1, yip1;

  // initial position: vi = 0
  nvtmp = nv.at(ci);
  gi = gindex(ci, 0);
  xi = x[NDIM * gi];
  yi = x[NDIM * gi + 1];

  // initialize center of mass coordinates
  cx = xi;
  cy = yi;

  // loop over vertices of cell ci, get perimeter
  for (vi = 0; vi < nvtmp - 1; vi++) {
    // next vertex
    gip1 = ip1.at(gi);
    gi++;

    // get positions (check minimum images)
    dx = x[NDIM * gip1] - xi;
    if (pbc[0])
      dx -= L[0] * round(dx / L[0]);
    xip1 = xi + dx;

    dy = x[NDIM * gip1 + 1] - yi;
    if (pbc[1])
      dy -= L[1] * round(dy / L[1]);
    yip1 = yi + dy;

    // add to center of mass
    cx += xip1;
    cy += yip1;

    // update coordinates
    xi = xip1;
    yi = yip1;
  }

  // take average to get com
  cx /= nvtmp;
  cy /= nvtmp;
}

// get configuration packing fraction
double dpm::vertexPackingFraction2D() {
  int ci;
  double val, boxV, areaSum = 0.0;

  // numerator
  for (ci = 0; ci < NCELLS; ci++) {
    areaSum += area(ci) + 0.25 * PI * pow(2.0 * r.at(szList[ci]), 2.0) * (0.5 * nv.at(ci) - 1);
  }

  // denominator
  boxV = L[0] * L[1];

  // return packing fraction
  val = areaSum / boxV;
  return val;
}

// get configuration "preferred" packing fraction
double dpm::vertexPreferredPackingFraction2D() {
  int ci;
  double val, boxV, areaSum = 0.0;

  // numerator
  for (ci = 0; ci < NCELLS; ci++)
    areaSum += a0[ci] + 0.25 * PI * pow(2.0 * r.at(szList[ci]), 2.0) * (0.5 * nv.at(ci) - 1);

  // denominator
  boxV = L[0] * L[1];

  // return packing fraction
  val = areaSum / boxV;
  return val;
}

// get configuration "preferred" packing fraction with respect to polygon boundaries poly_bd_x[i] and poly_bd_y[i], where we loop over i to get all boundaries
double dpm::vertexPreferredPackingFraction2D_polygon() {
  if (poly_bd_x.size() == 0) {
    cout << "poly_bd is empty, exiting in vertexPreferredPackingFraction2D_polygon!\n";
    return 0.0;
  }
  int ci;
  double val, boxV, boxV_temp, areaSum = 0.0;

  // numerator
  for (ci = 0; ci < NCELLS; ci++) {
    if (fabs(a0[ci] - a0[0]) < 2 * a0[0])
      areaSum += a0[ci] + 0.25 * PI * pow(2.0 * r.at(szList[ci]), 2.0) * (0.5 * nv.at(ci) - 1);
    // else
    //  cout << "area of cell " << ci << "is much larger than cell 0, assuming it's a boundary and not counting it in the areaSum of vertexPreferredPackingFraction2D_polygon\n";
  }

  // denominator = boundary polygonal area via shoelace method
  boxV = 0.0;
  for (int i = 0; i < poly_bd_x.size(); i++) {
    boxV_temp = 0.0;
    std::vector<double> poly_x = poly_bd_x[i];
    std::vector<double> poly_y = poly_bd_y[i];
    int j = poly_x.size() - 1;
    for (int k = 0; k < poly_x.size(); k++) {
      boxV_temp += (poly_x[j] + poly_x[k]) * (poly_y[j] - poly_y[k]);
      j = k;
    }
    boxV_temp = abs(boxV_temp / 2.0);
    boxV += boxV_temp;
  }

  // return packing fraction
  val = areaSum / boxV;
  cout << "packing fraction = " << val << ", boxV = " << boxV << ", areaSum = " << areaSum << '\n';
  cout << "estimated circle constructed from polygon has area = " << 0.25 * PI * L[0] * L[0] << '\n';
  return val;
}

// get vertex kinetic energy
double dpm::vertexKineticEnergy() {
  double K = 0;

  for (int i = 0; i < vertDOF; i++)
    K += v[i] * v[i];
  K *= 0.5;

  return K;
}

// get number of vertex-vertex contacts
int dpm::vvContacts() {
  int nvv = 0;

  for (int ci = 0; ci < NCELLS; ci++) {
    for (int cj = ci + 1; cj < NCELLS; cj++)
      nvv += cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
  }

  return nvv;
}

// get number of cell-cell contacts
int dpm::ccContacts() {
  int ncc = 0;

  for (int ci = 0; ci < NCELLS; ci++) {
    for (int cj = ci + 1; cj < NCELLS; cj++) {
      if (cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2] > 0)
        ncc++;
    }
  }

  return ncc;
}

/******************************

        I N I T I A L -

                        I Z A T I O N

*******************************/

void dpm::initializeFieldStress() {
  // local stress vector
  fieldStress.resize(NVTOT);
  for (int i = 0; i < NVTOT; i++) {
    fieldStress[i].resize(NDIM * (NDIM + 1) / 2);
    for (int j = 0; j < NDIM * (NDIM + 1) / 2; j++)
      fieldStress[i][j] = 0.0;
  }
  fieldShapeStress.resize(NVTOT);
  for (int i = 0; i < NVTOT; i++) {
    fieldShapeStress[i].resize(NDIM * (NDIM + 1) / 2);
    for (int j = 0; j < NDIM * (NDIM + 1) / 2; j++)
      fieldShapeStress[i][j] = 0.0;
  }
  fieldStressCells.resize(NCELLS);
  for (int i = 0; i < NCELLS; i++) {
    fieldStressCells[i].resize(NDIM * (NDIM + 1) / 2);
    for (int j = 0; j < NDIM * (NDIM + 1) / 2; j++)
      fieldStressCells[i][j] = 0.0;
  }
  fieldShapeStressCells.resize(NCELLS);
  for (int i = 0; i < NCELLS; i++) {
    fieldShapeStressCells[i].resize(NDIM * (NDIM + 1) / 2);
    for (int j = 0; j < NDIM * (NDIM + 1) / 2; j++)
      fieldShapeStressCells[i][j] = 0.0;
  }
}

// initialize vertex indexing
void dpm::initializeVertexIndexing2D() {
  int gi, vi, vip1, vim1, ci;

  // check that vertDOF has been assigned
  if (NVTOT <= 0) {
    cerr << "	** ERROR: in initializeVertexIndexing2D, NVTOT not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
    exit(1);
  }
  if (vertDOF <= 0) {
    cerr << "	** ERROR: in initializeVertexIndexing2D, vertDOF not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
    exit(1);
  } else if (nv.size() == 0) {
    cerr << "	** ERROR: in initializeVertexIndexing2D, nv vector not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
    exit(1);
  }

  // save list of adjacent vertices
  im1.resize(NVTOT);
  ip1.resize(NVTOT);
  for (ci = 0; ci < NCELLS; ci++) {
    // vertex indexing
    for (vi = 0; vi < nv.at(ci); vi++) {
      // wrap local indices
      vim1 = (vi - 1 + nv.at(ci)) % nv.at(ci);
      vip1 = (vi + 1) % nv.at(ci);

      // get global wrapped indices
      gi = gindex(ci, vi);
      im1.at(gi) = gindex(ci, vim1);
      ip1.at(gi) = gindex(ci, vip1);
    }
  }

  // initialize vertex configuration vectors
  x.resize(vertDOF);
  v.resize(vertDOF);
  F.resize(vertDOF);
}

// initialize vertex shape parameters and (a0, l0, t0, r) based on nv (nref is the reference nv, smallest nv among the polydispersity)
// sets a0 to 1 if nvtmp=nref, which is true for monodisperse2D()
void dpm::initializeVertexShapeParameters(double calA0, int nref) {
  // local variables
  int gi, ci, vi, nvtmp;
  double rtmp, calA0tmp, calAntmp;

  // check that vertDOF has been assigned
  if (NVTOT <= 0) {
    cerr << "	** ERROR: in initializeVertexShapeParameters, NVTOT not assigned. Ending here." << endl;
    exit(1);
  }
  if (vertDOF <= 0) {
    cerr << "	** ERROR: in initializeVertexShapeParameters, vertDOF not assigned. Ending here." << endl;
    exit(1);
  } else if (nv.size() == 0) {
    cerr << "	** ERROR: in initializeVertexShapeParameters, nv vector not assigned. Ending here." << endl;
    exit(1);
  }

  // resize shape paramters
  l0.resize(NVTOT);
  l00.resize(NVTOT);
  vl0.resize(NVTOT);
  Fl0.resize(NVTOT);
  t0.resize(NVTOT);
  r.resize(NVTOT);

  // loop over cells, determine shape parameters
  for (ci = 0; ci < NCELLS; ci++) {
    // number of vertices on cell ci
    nvtmp = nv.at(ci);

    // a0 based on nv
    rtmp = (double)nvtmp / nref;
    a0.at(ci) = rtmp * rtmp;

    // shape parameter
    calAntmp = nvtmp * tan(PI / nvtmp) / PI;
    calA0tmp = calA0 * calAntmp;

    // l0 and vertex radii
    gi = szList.at(ci);
    for (vi = 0; vi < nv.at(ci); vi++) {
      l0.at(gi + vi) = 2.0 * sqrt(PI * calA0tmp * a0.at(ci)) / nvtmp;
      vl0.at(gi + vi) = 0.0;
      Fl0.at(gi + vi) = 0.0;
      t0.at(gi + vi) = 0.0;
      r.at(gi + vi) = 0.5 * l0.at(gi + vi);
    }
  }
}

// initialize vertex shape parameters based on nv (nref is the reference nv, smallest nv among the polydispersity)
void dpm::initializeVertexShapeParameters(std::vector<double> calA0, int nref) {
  // local variables
  int gi, ci, vi, nvtmp;
  double rtmp, calA0tmp, calAntmp;

  // check that vertDOF has been assigned
  if (NVTOT <= 0) {
    cerr << "	** ERROR: in initializeVertexShapeParameters, NVTOT not assigned. Ending here." << endl;
    exit(1);
  }
  if (vertDOF <= 0) {
    cerr << "	** ERROR: in initializeVertexShapeParameters, vertDOF not assigned. Ending here." << endl;
    exit(1);
  } else if (nv.size() == 0) {
    cerr << "	** ERROR: in initializeVertexShapeParameters, nv vector not assigned. Ending here." << endl;
    exit(1);
  }

  // resize shape paramters
  l0.resize(NVTOT);
  l00.resize(NVTOT);
  vl0.resize(NVTOT);
  Fl0.resize(NVTOT);
  t0.resize(NVTOT);
  r.resize(NVTOT);

  // loop over cells, determine shape parameters
  for (ci = 0; ci < NCELLS; ci++) {
    // number of vertices on cell ci
    nvtmp = nv.at(ci);

    // a0 based on nv
    rtmp = (double)nvtmp / nref;
    a0.at(ci) = rtmp * rtmp;

    // shape parameter
    calAntmp = nvtmp * tan(PI / nvtmp) / PI;
    calA0tmp = calA0[ci] * calAntmp;

    // l0 and vertex radii
    gi = szList.at(ci);
    for (vi = 0; vi < nv.at(ci); vi++) {
      l0.at(gi + vi) = 2.0 * sqrt(PI * calA0tmp * a0.at(ci)) / nvtmp;
      vl0.at(gi + vi) = 0.0;
      Fl0.at(gi + vi) = 0.0;
      t0.at(gi + vi) = 0.0;
      r.at(gi + vi) = 0.5 * l0.at(gi + vi);
    }
  }
}

// calculate smallest x and y values and then shift simulation so they're both positive
void dpm::moveSimulationToPositiveCoordinates(double xshift, double yshift) {
  std::vector<double> posX(NVTOT / 2), posY(NVTOT / 2);
  for (int i = 0; i < NVTOT; i++) {
    if (i % 2 == 0)
      posX[i / 2] = x[NDIM * i];
    else
      posY[(i - 1) / 2] = x[NDIM * (i - 1) + 1];
  }

  double xLow = *std::min_element(posX.begin(), posX.end());
  // double xHigh = *std::max_element(posX.begin(), posX.end());
  double yLow = *std::min_element(posY.begin(), posY.end());
  // double yHigh = *std::max_element(posY.begin(), posY.end());

  if (xLow < 0)
    for (int gi = 0; gi < NVTOT; gi++)
      x[NDIM * gi] += fabs(xLow) + xshift;
  if (yLow < 0)
    for (int gi = 0; gi < NVTOT; gi++)
      x[NDIM * gi + 1] += fabs(yLow) + yshift;
}

// initialize monodisperse cell system, single calA0
void dpm::monodisperse2D(double calA0, int n) {
  int ci;

  // print to console
  cout << "** initializing monodisperse DPM particles in 2D ..." << endl;

  // total number of vertices
  NVTOT = n * NCELLS;
  vertDOF = NDIM * NVTOT;

  // szList and nv (keep track of global vertex indices)
  nv.resize(NCELLS);
  szList.resize(NCELLS);

  nv.at(0) = n;
  for (ci = 1; ci < NCELLS; ci++) {
    nv.at(ci) = n;
    szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
  }

  // initialize vertex shape parameters
  initializeVertexShapeParameters(calA0, n);

  // initialize vertex indexing
  initializeVertexIndexing2D();

  numVertexContacts.resize(NVTOT, std::vector<int>(NVTOT, 0.0));
}

// initialize bidisperse cell system, single calA0
void dpm::bidisperse2D(double calA0, int nsmall, double smallfrac, double sizefrac) {
  // local variables
  int ci, nlarge, smallN, largeN, NVSMALL;

  // print to console
  cout << "** initializing bidisperse DPM particles in 2D ..." << endl;

  // number of vertices on large particles
  nlarge = round(sizefrac * nsmall);

  // total number of vertices
  smallN = round(smallfrac * NCELLS);
  largeN = NCELLS - smallN;
  NVSMALL = nsmall * smallN;
  NVTOT = NVSMALL + nlarge * largeN;
  vertDOF = NDIM * NVTOT;

  // szList and nv (keep track of global vertex indices)
  nv.resize(NCELLS);
  szList.resize(NCELLS);

  nv.at(0) = nsmall;
  for (ci = 1; ci < NCELLS; ci++) {
    if (ci < smallN) {
      nv.at(ci) = nsmall;
      szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
    } else {
      nv.at(ci) = nlarge;
      szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
    }
  }

  // initialize vertex shape parameters
  initializeVertexShapeParameters(calA0, nsmall);

  // initialize vertex indexing
  initializeVertexIndexing2D();
}

// initialize gaussian polydisperse cell system, single calA0
void dpm::gaussian2D(double dispersion, double calA0, int n1) {
  // local variables
  double r1, r2, grv;
  int ci, nvtmp;

  // print to console
  cout << "** initializing gaussian DPM particles in 2D with size dispersion " << dispersion << " ..." << endl;

  // szList and nv (keep track of global vertex indices)
  nv.resize(NCELLS);
  szList.resize(NCELLS);

  nv.at(0) = n1;
  NVTOT = n1;
  for (ci = 1; ci < NCELLS; ci++) {
    // use Box-Muller to generate polydisperse sample
    r1 = drand48();
    r2 = drand48();
    grv = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
    nvtmp = floor(dispersion * n1 * grv + n1);
    if (nvtmp < nvmin)
      nvtmp = nvmin;

    // store size of cell ci
    nv.at(ci) = nvtmp;
    szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);

    // add to total NV count
    NVTOT += nvtmp;
  }
  vertDOF = NDIM * NVTOT;

  // initialize vertex shape parameters
  initializeVertexShapeParameters(calA0, n1);

  // initialize vertex indexing
  initializeVertexIndexing2D();
}

// set sinusoidal preferred angle
void dpm::sinusoidalPreferredAngle(double thA, double thK) {
  int ci, vi, gi;
  double thR;

  // print to console
  cout << "** setting initial th0 values to sinusoids, thA = " << thA << ", thK = " << thK << " ..." << endl;

  // loop over cells
  gi = 0;
  for (ci = 0; ci < NCELLS; ci++) {
    thR = (2.0 * PI) / nv.at(ci);
    for (vi = 0; vi < nv.at(ci); vi++) {
      t0.at(gi) = thA * thR * sin(thR * thK * vi);
      gi++;
    }
  }
}

// initialize CoM positions of cells (i.e. use soft disks) using SP FIRE. setupCircularBoundaries enables polygonal walls
void dpm::initializePositions2D(double phi0, double Ftol, bool isFixedBoundary, double aspectRatio, bool setUpCircularBoundary) {
  // isFixedBoundary is an optional bool argument that tells cells to stay away from the boundary during initialization
  // aspectRatio is the ratio L[0] / L[1]
  int i, d, ci, cj, vi, gi, cellDOF = NDIM * NCELLS;
  int numEdges = 20;  // number of edges in the polygonal walls to approximate a circle
  double areaSum, xtra = 1.0;
  std::vector<double> aspects = {1.0 * aspectRatio, 1.0 / aspectRatio};

  // local disk vectors
  vector<double> drad(NCELLS, 0.0);
  vector<double> dpos(cellDOF, 0.0);
  vector<double> dv(cellDOF, 0.0);
  vector<double> dF(cellDOF, 0.0);

  // print to console
  cout << "** initializing particle positions using 2D SP model and FIRE relaxation ..." << endl;

  // initialize stress field
  initializeFieldStress();

  // initialize box size based on packing fraction
  areaSum = 0.0;
  for (ci = 0; ci < NCELLS; ci++)
    areaSum += a0.at(ci) + 0.25 * PI * pow(l0.at(ci), 2.0) * (0.5 * nv.at(ci) - 1);

  // set box size : phi_0 = areaSum / A => A = areaSum/phi_0 which gives us the following formulas for L
  for (d = 0; d < NDIM; d++) {
    L.at(d) = pow(areaSum / phi0, 1.0 / NDIM) * aspects[d];
    if (setUpCircularBoundary)
      L.at(d) = pow(4 / PI * areaSum / phi0, 1.0 / NDIM);
  }

  // initialize cell centers randomly
  if (!setUpCircularBoundary) {
    for (ci = 0; ci < cellDOF; ci += 2) {
      dpos.at(ci) = L[ci % 2] * drand48();
    }
    for (ci = cellDOF - 1; ci > 0; ci -= 2) {
      dpos.at(ci) = L[ci % 2] * drand48();
    }
  } else {
    cout << "setUpCircularBoundary is enabled, so initializing cell centers randomly but rejecting if further than R = L/2 from the center (which is (L/2,L/2))\n";
    double scale_radius = 1.0;                   // make the polygon radius slightly larger so that it encompasses the circle that points are initialized in
    poly_bd_x.push_back(std::vector<double>());  // make new data for generateCircularBoundary to write a polygon
    poly_bd_y.push_back(std::vector<double>());
    double cx = L[0] / 2, cy = L[1] / 2, tissueRadius = L[0] / 2;
    ofstream boundaryStream("polyBoundary.txt");  // clear boundary text file, for visualizing or for debugging
    generateCircularBoundary(numEdges, scale_radius * tissueRadius, cx, cy, poly_bd_x[0], poly_bd_y[0]);

    for (i = 0; i < NCELLS; i++) {
      double dpos_x, dpos_y;
      bool insidePolygon = false;
      while (!insidePolygon) {
        dpos_x = (tissueRadius - 4 * r[0]) * (2 * drand48() - 1) + cx;
        dpos_y = (tissueRadius - 4 * r[0]) * (2 * drand48() - 1) + cy;
        // Check if the generated position is inside the polygon
        insidePolygon = isInsidePolygon(dpos_x, dpos_y, poly_bd_x[0], poly_bd_y[0]);
      }
      dpos.at(i * NDIM) = dpos_x;
      dpos.at(i * NDIM + 1) = dpos_y;
    }
  }

  // set radii of SP disks
  for (ci = 0; ci < NCELLS; ci++) {
    if (setUpCircularBoundary)
      xtra = 1.0;  // disks should have radius similar to the final particle radius, or could modify vrad[i] condition in wall calculation later
    drad.at(ci) = xtra * sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin(2.0 * PI / nv.at(ci))));
  }

  // FIRE VARIABLES
  double P = 0;
  double fnorm = 0;
  double vnorm = 0;
  double alpha = alpha0;

  double dt0 = 0.01;
  double dtmax = 10 * dt0;
  double dtmin = 1e-8 * dt0;

  int npPos = 0;
  int npNeg = 0;

  int fireit = 0;
  double fcheck = 10 * Ftol;

  // interaction variables
  double rij, sij, dtmp, ftmp, vftmp;
  double dr[NDIM];

  // initial step size
  dt = dt0;

  // loop until force relaxes
  while ((fcheck > Ftol) && fireit < itmax) {
    // FIRE step 1. Compute P
    P = 0.0;
    for (i = 0; i < cellDOF; i++)
      P += dv[i] * dF[i];

    // FIRE step 2. adjust simulation based on net motion of degrees of freedom
    if (P > 0) {
      // increase positive counter
      npPos++;

      // reset negative counter
      npNeg = 0;

      // alter simulation if enough positive steps have been taken
      if (npPos > NMIN) {
        // change time step
        if (dt * finc < dtmax)
          dt *= finc;

        // decrease alpha
        alpha *= falpha;
      }
    } else {
      // reset positive counter
      npPos = 0;

      // increase negative counter
      npNeg++;

      // check if simulation is stuck
      if (npNeg > NNEGMAX) {
        cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
        exit(1);
      }

      // take half step backwards, reset velocities
      for (i = 0; i < cellDOF; i++) {
        // take half step backwards
        dpos[i] -= 0.5 * dt * dv[i];

        // reset velocities
        dv[i] = 0.0;
      }

      // decrease time step if past initial delay
      if (fireit > NDELAY) {
        // decrease time step
        if (dt * fdec > dtmin)
          dt *= fdec;

        // reset alpha
        alpha = alpha0;
      }
    }

    // FIRE step 3. First VV update
    for (i = 0; i < cellDOF; i++)
      dv[i] += 0.5 * dt * dF[i];

    // FIRE step 4. adjust velocity magnitude
    fnorm = 0.0;
    vnorm = 0.0;
    for (i = 0; i < cellDOF; i++) {
      fnorm += dF[i] * dF[i];
      vnorm += dv[i] * dv[i];
    }
    fnorm = sqrt(fnorm);
    vnorm = sqrt(vnorm);
    if (fnorm > 0) {
      for (i = 0; i < cellDOF; i++)
        dv[i] = (1 - alpha) * dv[i] + alpha * (vnorm / fnorm) * dF[i];
    }

    // FIRE step 4. Second VV update
    for (i = 0; i < cellDOF; i++) {
      dpos[i] += dt * dv[i];
      dF[i] = 0.0;
    }

    // FIRE step 5. Update forces
    for (ci = 0; ci < NCELLS; ci++) {
      for (cj = ci + 1; cj < NCELLS; cj++) {
        // contact distance
        sij = drad[ci] + drad[cj];

        // true distance
        rij = 0.0;
        for (d = 0; d < NDIM; d++) {
          // get distance element
          dtmp = dpos[NDIM * cj + d] - dpos[NDIM * ci + d];
          if (pbc[d])
            dtmp -= L[d] * round(dtmp / L[d]);

          // add to true distance
          rij += dtmp * dtmp;

          // save in distance array
          dr[d] = dtmp;
        }
        rij = sqrt(rij);

        // check distances
        if (rij < sij) {
          // force magnitude
          ftmp = kc * (1.0 - (rij / sij)) / sij;

          // add to vectorial force
          for (d = 0; d < NDIM; d++) {
            vftmp = ftmp * (dr[d] / rij);
            dF[NDIM * ci + d] -= vftmp;
            dF[NDIM * cj + d] += vftmp;
          }
        }
      }
    }
    // FIRE step 4.1 Compute wall forces
    if (setUpCircularBoundary) {
      for (int k = 0; k < poly_bd_x.size(); k++) {
        std::vector<double> poly_x = poly_bd_x[k];
        std::vector<double> poly_y = poly_bd_y[k];
        int n = poly_x.size();
        double distanceParticleWall, Rx, Ry, dw, K = 1;
        double bound_x1, bound_x2, bound_y1, bound_y2;
        // loop over boundary bars
        // loop over particles
        //  compute particle-boundary bar overlaps
        //  if overlap, Fx += K * dw * Rx/R, where K is a constant, dw = diameter/2 - R, Rx = x - px, R = sqrt(Rx^2 + Ry^2)
        for (int bound_i = 0; bound_i < n; bound_i++) {
          // use distanceLineAndPoint to get R, Rx, and Ry
          bound_x1 = poly_x[bound_i];
          bound_x2 = poly_x[(bound_i + 1) % n];
          bound_y1 = poly_y[bound_i];
          bound_y2 = poly_y[(bound_i + 1) % n];
          for (i = 0; i < cellDOF / NDIM; i++) {
            distanceParticleWall = distanceLinePointComponents(bound_x1, bound_y1, bound_x2, bound_y2, dpos[i * NDIM], dpos[i * NDIM + 1], Rx, Ry);
            dw = drad[i] - distanceParticleWall;
            if (distanceParticleWall <= drad[i]) {
              dF[i * NDIM] += K * dw * Rx / distanceParticleWall;
              dF[i * NDIM + 1] += K * dw * Ry / distanceParticleWall;
            }
          }
        }
      }
    } else if (isFixedBoundary) {
      for (i = 0; i < cellDOF; i++) {
        bool collideTopOrRight = dpos[i] > L[i % NDIM] - drad[i];
        bool collideBottomOrLeft = dpos[i] < drad[i];

        if (collideTopOrRight) {  // deflect particle down or left
          dF[i] += -1 * (drad[i] - L[i % NDIM] + dpos[i]);
        }
        if (collideBottomOrLeft) {
          dF[i] += 1 * (drad[i] - dpos[i]);
        }
      }
    }

    // FIRE step 5. Final VV update
    for (i = 0; i < cellDOF; i++)
      dv[i] += 0.5 * dt * dF[i];

    // update forces to check
    fcheck = 0.0;
    for (i = 0; i < cellDOF; i++)
      fcheck += dF[i] * dF[i];
    fcheck = sqrt(fcheck / NCELLS);

    // print to console
    if (fireit % NSKIP == 0) {
      cout << endl
           << endl;
      cout << "===========================================" << endl;
      cout << "		I N I T I A L  S P 			" << endl;
      cout << " 	F I R E 						" << endl;
      cout << "		M I N I M I Z A T I O N 	" << endl;
      cout << "===========================================" << endl;
      cout << endl;
      cout << "	** fireit = " << fireit << endl;
      cout << "	** fcheck = " << fcheck << endl;
      cout << "	** fnorm = " << fnorm << endl;
      cout << "	** vnorm = " << vnorm << endl;
      cout << "	** dt = " << dt << endl;
      cout << "	** P = " << P << endl;
      cout << "	** Pdir = " << P / (fnorm * vnorm) << endl;
      cout << "	** alpha = " << alpha << endl;
    }

    // update iterate
    fireit++;
  }
  // check if FIRE converged
  if (fireit == itmax) {
    cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
    exit(1);
  } else {
    cout << endl
         << endl;
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
    cout << "===========================================" << endl;
    cout << " 	F I R E 						" << endl;
    cout << "		M I N I M I Z A T I O N 	" << endl;
    cout << "	C O N V E R G E D! 				" << endl
         << endl;

    cout << "	(for initial disk minimization) " << endl;
    cout << "===========================================" << endl;
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
    cout << endl;
    cout << "	** fireit = " << fireit << endl;
    cout << "	** fcheck = " << fcheck << endl;
    cout << "	** vnorm = " << vnorm << endl;
    cout << "	** dt = " << dt << endl;
    cout << "	** P = " << P << endl;
    cout << "	** alpha = " << alpha << endl;
  }

  // initialize vertex positions based on cell centers
  for (ci = 0; ci < NCELLS; ci++) {
    for (vi = 0; vi < nv.at(ci); vi++) {
      // get global vertex index
      gi = gindex(ci, vi);

      // length from center to vertex minus 1/2 the vertex radius to prevent overlaps
      dtmp = sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin((2.0 * PI) / nv.at(ci)))) - r[gi] / 2.0;

      // set positions
      x.at(NDIM * gi) = dtmp * cos((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci) + 1e-2 * l0[gi] * drand48();
      x.at(NDIM * gi + 1) = dtmp * sin((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci + 1) + 1e-2 * l0[gi] * drand48();
    }
  }
}

void dpm::initializeFromConfigurationFile(std::string vertexPositionFile, double phi0) {
  // in case of variable calA0, nv, and positions, this function subsumes monodisperse2D

  std::ifstream positionFile(vertexPositionFile);
  int cellNum, vertNum, a, b;
  double vertLocX, vertLocY;
  std::vector<int> numVertsPerDPM;
  std::vector<double> vertexPositions;

  while (positionFile >> cellNum >> vertNum) {  // header line is cell ID and total # vertices
    numVertsPerDPM.push_back(vertNum);
    for (int i = 0; i < vertNum; i++) {
      positionFile >> a >> b >> vertLocX >> vertLocY;
      vertexPositions.push_back(vertLocX);
      vertexPositions.push_back(vertLocY);
    }
  }

  // total number of vertices
  NCELLS = cellNum;
  cout << "NCELLS = " << cellNum << '\n';

  NVTOT = 0;
  for (auto i : numVertsPerDPM) {
    NVTOT += i;
  }
  vertDOF = NDIM * NVTOT;

  // szList and nv keep track of global vertex indices
  nv.resize(NCELLS);
  szList.resize(NCELLS);

  nv.at(0) = numVertsPerDPM[0];
  for (int ci = 1; ci < NCELLS; ci++) {
    nv.at(ci) = numVertsPerDPM[ci];
    szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
  }

  // initialize connectivity between vertices in DPs
  initializeVertexIndexing2D();

  // minimum coordinates to subtract so particles start near origin
  double min_x = 1e10, min_y = 1e10;
  for (int i = 0; i < vertexPositions.size(); i += 2) {
    if (vertexPositions[NDIM * i] < min_x)
      min_x = vertexPositions[NDIM * i];
    if (vertexPositions[NDIM * i + 1] < min_y)
      min_y = vertexPositions[NDIM * i + 1];
  }

  int gi;
  // initialize vertex positions
  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv.at(ci); vi++) {
      gi = gindex(ci, vi);
      x.at(gi * NDIM) = vertexPositions[gi * NDIM] - min_x;
      x.at(gi * NDIM + 1) = vertexPositions[gi * NDIM + 1] - min_y;
    }
  }

  // once i've assigned cell IDs and vertex coordinates, use dpm functions to set the preferred lengths, radii

  // resize shape paramters
  l0.resize(NVTOT);
  l00.resize(NVTOT);
  vl0.resize(NVTOT);
  Fl0.resize(NVTOT);
  t0.resize(NVTOT);
  r.resize(NVTOT);

  double l0_all, r_all, d_all;  // variables for rescaling lengths
  // loop over cells, determine shape parameters
  for (int ci = 0; ci < NCELLS; ci++) {
    double l0_ci, r_ci;
    for (int vi = 0; vi < nv.at(ci); vi++) {
      gi = gindex(ci, vi);
      // l0.at(gi) = sqrt(pow(x.at(NDIM * ip1[gi]) - x.at(NDIM * gi), 2) + pow(x.at(NDIM * ip1[gi] + 1) - x.at(NDIM * gi + 1), 2));
      // save parameters from the first vertex of each cell to set l0, r
      if (vi == 0) {
        l0_ci = sqrt(pow(x.at(NDIM * ip1[gi]) - x.at(NDIM * gi), 2) + pow(x.at(NDIM * ip1[gi] + 1) - x.at(NDIM * gi + 1), 2));
        r_ci = 0.5 * l0_ci;
        if (ci == 0) {
          // save parameters from the first vertex of the first cell to set length renormalization scale
          l0_all = l0_ci;
          r_all = l0_ci / 2.0;
          d_all = l0_ci;
        }
      }
      l0.at(gi) = l0_ci / d_all;
      vl0.at(gi) = 0.0;
      Fl0.at(gi) = 0.0;
      t0.at(gi) = 0.0;
      r.at(gi) = r_ci / d_all;
      cout << "r[gi] at time of setting = " << r[gi] << '\n';
      x[NDIM * gi] /= d_all;
      x[NDIM * gi + 1] /= d_all;
    }
  }

  for (int d = 0; d < NDIM; d++) {
    L.at(d) = 1e10;  // set L to be large to not interfere with initial shape calculations, reset L later
  }

  double areaSum = 0.0;
  for (int ci = 0; ci < NCELLS; ci++) {
    areaSum += area(ci);
    a0.at(ci) = area(ci);
    cout << "new shape subject to length rescaling is = " << pow(perimeter(ci), 2) / (4 * 3.1415 * area(ci)) << '\n';
  }
  // todo: calculate area of all DP cells, then use input phi to give myself a square box. then run testInterpolatedConfiguration.m to see the initial configuration and judge its quality

  // set box size
  for (int d = 0; d < NDIM; d++) {
    // L.at(d) = pow(areaSum / phi0, 1.0 / NDIM);
    double max = *max_element(x.begin(), x.end());
    L.at(d) = 1.1 * max;
  }
}

// initialize vertices according to an input configuration file, and initialize/set relevant lists and degrees of freedom variables.
void dpm::initializeAllPositions(std::string vertexPositionFile, int nref) {
  // for every block in vertexPositionFile, place vertices according to a particle.
  // vertexPositionFile should have NCELLS blocks separated by headers
  // each block should have a header
  //      'nan nan'
  //      'N(particle label) NV(number of vertices, raw)'
  //
  // every successive line should have x and y positions for vertices belonging to the cell
  // every line (x, y) corresponds to a vertex position, spaced 1-pixel away from the previous line
  // start by calculating the perimeter, area, shape of the block

  double posx, posy;
  double pixelsPerParticle;
  int N = -1, numVertsCounter = 0, lineCounter = 0;
  std::vector<double> vertexPositions;
  std::vector<int> numVerts(NCELLS);
  std::vector<double> shapes;

  // read in file
  ifstream positionFile(vertexPositionFile);
  while (positionFile >> posx >> posy) {
    if (std::isnan(posx)) {  // found a header, so read in the next line and calculate how many vertices we want
      N++;
      positionFile >> posx >> posy;  // particle label and pixels in perimeter
      if (N == 0) {
        pixelsPerParticle = posy / double(nref);  // since N = 0 corresponds to the smallest perimeter, set this as baseline
      }
      if (N > 0) {  // need to store vertexPositions and other cell related quantities, since we are moving to a new cell
        numVerts.push_back(numVertsCounter);
      }
      // reset counters
      numVertsCounter = 0;
      lineCounter = 0;
      continue;
    }
    if (lineCounter == round(pixelsPerParticle)) {  // we have read in the coordinates of the vertex we want to initialize
      vertexPositions.push_back(posx);
      vertexPositions.push_back(posy);
      numVertsCounter++;
    }
    lineCounter++;
  }
  double perimeter = 0;
  double area = 0;

  for (int i = 0; i < N; i++) {
    // calculate perimeter
    // calculate area
    // shapes.push_back(perimeter^2/4pi*area);
  }

  // might need to set box length at some point here

  // set NVTOT, vertDOF, nv, szList, vertex lists
  // NCELLS = N;
  // NVTOT = accumulate(numVerts);
  // vertDOF = NVTOT * NDIM;
  // nv.resize(NCELLS);
  // szList.resize(NCELLS);
  // nv.at(0) = nref?;
  // for (ci = 1; ci < NCELLS; ci++){
  //  nv.at(ci) = nref;
  //  szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);
  // }
  //

  initializeVertexShapeParameters(shapes, nref);  // overloaded this function, need to check later if it compiles

  initializeVertexIndexing2D();

  int gi;
  // initialize vertex positions
  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv.at(ci); vi++) {
      gi = gindex(ci, vi);
      x.at(gi * NDIM) = vertexPositions[gi * NDIM];  // not sure about this, double check RHS
      x.at(gi * NDIM + 1) = vertexPositions[gi * NDIM + 1];
    }
  }
}

// initialize neighbor linked list
void dpm::initializeNeighborLinkedList2D(double boxLengthScale) {
  // local variables
  double llscale;
  int i, d, nntmp, scx;

  // print to console
  cout << "** initializing neighbor linked list, boxLengthScale = " << boxLengthScale << '\n';

  // get largest diameter times attraction shell (let buffer be larger than attraction range would ever be) as llscale
  double buffer = 1.5;
  llscale = buffer * 2 * (*max_element(r.begin(), r.end()));
  cout << "llscale = " << llscale << '\n';

  // initialize box length vectors
  NBX = 1;

  sb.resize(NDIM);
  lb.resize(NDIM);
  for (d = 0; d < NDIM; d++) {
    // determine number of cells along given dimension by rmax
    sb[d] = round(L[d] / (boxLengthScale * llscale));

    // just in case, if < 3, change to 3 so box neighbor checking will work
    if (sb[d] < 3)
      sb[d] = 3;

    // determine box length by number of cells
    lb[d] = L[d] / sb[d];

    // count total number of cells
    NBX *= sb[d];
  }

  // initialize list of box nearest neighbors
  scx = sb[0];
  nn.resize(NBX);

  // loop over cells, save forward neighbors for each box
  for (i = 0; i < NBX; i++) {
    // reshape entry
    nn[i].resize(NNN);

    // right neighbor (i+1)
    nn[i][0] = (i + 1) % NBX;

    // top neighbors (i,j+1), (i+1,j+1)
    if (pbc[1]) {
      // (i,j+1) w/ pbc
      nn[i][1] = (i + scx) % NBX;

      // (i+1,j+1) w/ pbc
      nn[i][2] = (nn[i][1] + 1) % NBX;
    } else {
      // if on top row, both = -1
      if (i >= NBX - scx) {
        nn[i][1] = -1;
        nn[i][2] = -1;
      }
      // if not on top row, still add
      else {
        nn[i][1] = i + scx;
        nn[i][2] = nn[i][1] + 1;
      }
    }

    // bottom neighbor w/ pbc (j-1)
    nntmp = (i + NBX - scx) % NBX;

    // bottom-right neighbor (i+1, j-1)
    if (pbc[1])
      nn[i][3] = nntmp + 1;
    else {
      // if on bottom row, skip
      if (i < scx)
        nn[i][3] = -1;
      // otherwise, set
      else
        nn[i][3] = nntmp + 1;
    }

    // right-hand bc (periodic)
    if ((i + 1) % scx == 0) {
      if (pbc[0]) {
        nn[i][0] = i - scx + 1;
        if (pbc[1]) {
          nn[i][2] = nn[i][1] - scx + 1;
          nn[i][3] = nntmp - scx + 1;
        }
      } else {
        nn[i][0] = -1;
        nn[i][2] = -1;
        nn[i][3] = -1;
      }
    }
  }

  // linked-list variables
  head.resize(NBX);
  last.resize(NBX);
  list.resize(NVTOT + 1);

  // print box info to console
  cout << ";  initially NBX = " << NBX << " ..." << endl;
}

// resize neighbor linked list based on current position of particles
void dpm::resizeNeighborLinkedList2D() {
  // local variables
  double llscale;
  int i, d, nntmp, scx;

  // print to console
  double boxLengthScale = 2.5;
  cout << "** resizing neighbor linked list, lengthscale = " << boxLengthScale << '\n';

  if (!pbc[0] && !pbc[1]) {
    double maxX = 0, maxY = 0;
    for (int i = 0; i < x.size(); i += 2) {
      if (x[i] > maxX)
        maxX = x[i];
      if (x[i + 1] > maxY)
        maxY = x[i + 1];
    }
    // double currentBoxLength = *std::max_element(x.begin(), x.end());
    cout << "assigning max particle dimensions to be the new box length to size the neighbor list properly\n";
    L[0] = maxX * 1.05;
    L[1] = maxY * 1.05;
  }

  // get largest diameter times attraction shell (let buffer be larger than attraction range would ever be) as llscale
  double buffer = 1.5;
  llscale = buffer * 2 * max((*max_element(r.begin(), r.end())), *max_element(l0.begin(), l0.end()));
  cout << "llscale = " << llscale << '\n';

  // initialize box length vectors
  NBX = 1;
  sb.resize(NDIM);
  lb.resize(NDIM);
  for (d = 0; d < NDIM; d++) {
    // determine number of cells along given dimension by rmax
    sb[d] = round(L[d] / (boxLengthScale * llscale));

    // just in case, if < 3, change to 3 so box neighbor checking will work
    if (sb[d] < 3)
      sb[d] = 3;

    // determine box length by number of cells
    lb[d] = L[d] / sb[d];

    // count total number of cells
    NBX *= sb[d];
  }

  // initialize list of box nearest neighbors
  scx = sb[0];
  nn.resize(NBX);

  // loop over cells, save forward neighbors for each box
  for (i = 0; i < NBX; i++) {
    // reshape entry
    nn[i].resize(NNN);

    // right neighbor (i+1)
    nn[i][0] = (i + 1) % NBX;

    // top neighbors (i,j+1), (i+1,j+1)
    if (pbc[1]) {
      // (i,j+1) w/ pbc
      nn[i][1] = (i + scx) % NBX;

      // (i+1,j+1) w/ pbc
      nn[i][2] = (nn[i][1] + 1) % NBX;
    } else {
      // if on top row, both = -1
      if (i >= NBX - scx) {
        nn[i][1] = -1;
        nn[i][2] = -1;
      }
      // if not on top row, still add
      else {
        nn[i][1] = i + scx;
        nn[i][2] = nn[i][1] + 1;
      }
    }

    // bottom neighbor w/ pbc (j-1)
    nntmp = (i + NBX - scx) % NBX;

    // bottom-right neighbor (i+1, j-1)
    if (pbc[1])
      nn[i][3] = nntmp + 1;
    else {
      // if on bottom row, skip
      if (i < scx)
        nn[i][3] = -1;
      // otherwise, set
      else
        nn[i][3] = nntmp + 1;
    }

    // right-hand bc (periodic)
    if ((i + 1) % scx == 0) {
      if (pbc[0]) {
        nn[i][0] = i - scx + 1;
        if (pbc[1]) {
          nn[i][2] = nn[i][1] - scx + 1;
          nn[i][3] = nntmp - scx + 1;
        }
      } else {
        nn[i][0] = -1;
        nn[i][2] = -1;
        nn[i][3] = -1;
      }
    }
  }

  // linked-list variables
  head.resize(NBX);
  last.resize(NBX);
  list.resize(NVTOT + 1);

  // print box info to console
  cout << ";  initially NBX = " << NBX << " ..." << endl;
}

/******************************

        E D I T I N G   &

                        U P D A T I N G

*******************************/

// sort vertices into neighbor linked list
void dpm::sortNeighborLinkedList2D() {
  // local variables
  int d, gi, boxid, sbtmp;
  double xtmp;

  /*cout << "before neighborLinkedList\n";
  cout << "list = \n";
  for (auto i : list)
    cout << i << '\t';
  cout << "\nhead = \n";
  for (auto i : head)
    cout << i << '\t';
  cout << "\nlast = \n";
  for (auto i : last)
    cout << i << '\t';*/

  // reset linked list info
  fill(list.begin(), list.end(), 0);
  fill(head.begin(), head.end(), 0);
  fill(last.begin(), last.end(), 0);
  // sort vertices into linked list
  for (gi = 0; gi < NVTOT; gi++) {
    // 1. get cell id of current particle position
    boxid = 0;
    sbtmp = 1;
    for (d = 0; d < NDIM; d++) {
      // current location
      xtmp = x[NDIM * gi + d];
      // check out-of-bounds
      if (xtmp < 0) {
        if (pbc[d])
          xtmp -= L[d] * floor(xtmp / L[d]);
        else
          xtmp = 0.00001;
      } else if (xtmp > L[d]) {
        if (pbc[d])
          xtmp -= L[d] * floor(xtmp / L[d]);
        else  // if pbc is off, whenever I call resizeNeighborLinkedList I also shift L to be the max dimensions of the system with a little wiggle room
          xtmp = 0.99999 * L[d];
      }

      // add d index to 1d list
      boxid += floor(xtmp / lb[d]) * sbtmp;
      if (boxid < -2147483600) {
        cout << "boxid = " << boxid << ", xtmp = " << xtmp << ", lb[d] = " << lb[d] << ", d = " << d << ", gi = " << gi << ", floor(xtmp / lb[d]) = " << floor(xtmp / lb[d]) << ", list.size = " << list.size() << '\n';
        cout << "x[NDIM * gi + d (before modification)] = " << x[NDIM * gi + d] << ", x.size() = " << x.size() << '\n';
        cout << "pbc[d] = " << pbc[d] << ", L[d] = " << L[d] << '\n';

        for (int k = 0; k < NVTOT; k++) {
          cout << "vert " << k << ": " << x[NDIM * k] << ", " << x[NDIM * k + 1] << '\n';
        }
      }

      // increment dimensional factor
      sbtmp *= sb[d];
    }

    // 2. add to head list or link within list
    // NOTE: particle ids are labelled starting from 1, setting to 0 means end of linked list
    if (head[boxid] == 0) {
      head[boxid] = gi + 1;
      last[boxid] = gi + 1;
    } else {
      if (boxid >= last.size()) {
        cout << "gi " << gi << ", boxid out of range : " << boxid << '\t' << last.size() << '\n';
        cout << "x out of range = " << x[NDIM * gi] << '\t' << x[NDIM * gi + 1] << '\n';
        cout << "L = " << L[0] << '\t' << L[1] << '\n';
      }
      if (last[boxid] >= list.size()) {
        cout << "gi " << gi << ", last[boxid] out of range : " << last[boxid] << '\t' << list.size() << '\n';
      }
      list[last[boxid]] = gi + 1;
      last[boxid] = gi + 1;
    }
  }
  /*cout << "after neighborLinkedList\n";
  cout << "list = \n";
  for (auto i : list)
    cout << i << '\t';
  cout << "\nhead = \n";
  for (auto i : head)
    cout << i << '\t';
  cout << "\nlast = \n";
  for (auto i : last)
    cout << i << '\t';*/
}

// change size of particles
void dpm::scaleParticleSizes2D(double scaleFactor) {
  // local variables
  int gi, ci, vi, xind, yind;
  double xi, yi, cx, cy, dx, dy;

  // loop over cells, scale
  for (ci = 0; ci < NCELLS; ci++) {
    // scale preferred area
    a0[ci] *= scaleFactor * scaleFactor;

    // first global index for ci
    gi = szList.at(ci);

    // compute cell center of mass
    xi = x[NDIM * gi];
    yi = x[NDIM * gi + 1];
    cx = xi;
    cy = yi;
    for (vi = 1; vi < nv.at(ci); vi++) {
      dx = x.at(NDIM * (gi + vi)) - xi;
      if (pbc[0])
        dx -= L[0] * round(dx / L[0]);

      dy = x.at(NDIM * (gi + vi) + 1) - yi;
      if (pbc[1])
        dy -= L[1] * round(dy / L[1]);

      xi += dx;
      yi += dy;

      cx += xi;
      cy += yi;
    }
    cx /= nv.at(ci);
    cy /= nv.at(ci);

    for (vi = 0; vi < nv.at(ci); vi++) {
      // x and y inds
      xind = NDIM * (gi + vi);
      yind = xind + 1;

      // closest relative position
      dx = x[xind] - cx;
      if (pbc[0])
        dx -= L[0] * round(dx / L[0]);

      dy = x[yind] - cy;
      if (pbc[1])
        dy -= L[1] * round(dy / L[1]);

      // update vertex positions
      x[xind] += (scaleFactor - 1.0) * dx;
      x[yind] += (scaleFactor - 1.0) * dy;

      // scale vertex radii
      r[gi + vi] *= scaleFactor;
      l0[gi + vi] *= scaleFactor;
    }
  }
}

// remove rattlers from contact network, return number of rattlers
int dpm::removeRattlers() {
  // local variables
  int ci, cj, ctmp, rvv, rcc, nr, nm = 1;

  // loop over rows, eliminate contacts to rattlers until nm = 0
  while (nm > 0) {
    // reset number of rattlers
    nr = 0;

    // number of "marginal" rattlers to be removed
    nm = 0;
    for (ci = 0; ci < NCELLS; ci++) {
      // get number of contacts on cell ci
      rvv = 0;
      rcc = 0;
      for (cj = 0; cj < NCELLS; cj++) {
        if (ci != cj) {
          if (ci > cj)
            ctmp = cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];
          else
            ctmp = cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
        } else
          ctmp = 0;

        rvv += ctmp;
        if (ctmp > 0)
          rcc++;
      }

      // check to see if particle should be removed from network
      if (rcc <= NDIM && rvv <= 3) {
        // increment # of rattlers
        nr++;

        // if in contact, remove contacts
        if (rvv > 0) {
          nm++;

          for (cj = 0; cj < NCELLS; cj++) {
            // delete contact between ci and cj
            if (ci != cj) {
              if (ci > cj)
                cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2] = 0;
              else
                cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2] = 0;
            }
          }
        }
      }
    }
  }

  // return total number of rattlers
  return nr;
}

// draw random velocities based on input temperature
void dpm::drawVelocities2D(double T) {
  // local variables
  int gi;
  double r1, r2, grv1, grv2, tscale = sqrt(T), vcomx = 0.0, vcomy = 0.0;

  // loop over velocities, draw from maxwell-boltzmann distribution
  for (gi = 0; gi < NVTOT; gi++) {
    // draw random numbers using Box-Muller
    r1 = drand48();
    r2 = drand48();
    grv1 = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
    grv2 = sqrt(-2.0 * log(r1)) * sin(2.0 * PI * r2);

    // assign to velocities
    v[NDIM * gi] = tscale * grv1;
    v[NDIM * gi + 1] = tscale * grv2;

    // add to center of mass
    vcomx += v[NDIM * gi];
    vcomy += v[NDIM * gi + 1];
  }
  vcomx = vcomx / NVTOT;
  vcomy = vcomy / NVTOT;

  // subtract off center of mass drift
  for (gi = 0; gi < NVTOT; gi++) {
    v[NDIM * gi] -= vcomx;
    v[NDIM * gi + 1] -= vcomy;
  }
}

double dpm::distanceLineAndPoint(double x1, double y1, double x2, double y2, double x0, double y0) {
  // get the distance from a line segment going through (x1,y1), (x2,y2) and a
  // point located at (x0,y0)
  double dx21 = x2 - x1, dy21 = y2 - y1, dx10 = x1 - x0, dy10 = y1 - y0;
  if (pbc[0]) {
    dx21 -= L[0] * round(dx21 / L[0]);
    dx10 -= L[0] * round(dx10 / L[0]);
  }
  if (pbc[1]) {
    dy21 -= L[1] * round(dy21 / L[1]);
    dy10 -= L[1] * round(dy10 / L[1]);
  }

  double l2 = pow(dx21, 2) + pow(dy21, 2);  // |(pt2 - pt1)|^2
  if (l2 == 0.0)                            // pt2 == pt1 case
    return sqrt(pow(dx10, 2) + pow(dy10, 2));

  double dot = (-dx10) * (dx21) +
               (-dy10) * (dy21);  // (pt0 - pt1) dot (pt2 - pt1)
  const double t = max(0.0, min(1.0, dot / l2));
  const double projectionx = x1 + t * (dx21);
  const double projectiony = y1 + t * (dy21);
  const double distance =
      sqrt(pow(x0 - projectionx, 2) + pow(y0 - projectiony, 2));
  return distance;
}

double dpm::distanceLinePointComponents(double x1, double y1, double x2, double y2, double x0, double y0, double& xcomp, double& ycomp) {
  // get the distance from a line segment going through (x1,y1), (x2,y2) and a
  // point located at (x0,y0), and extract x and y components of the distance
  double dx21 = x2 - x1, dy21 = y2 - y1, dx10 = x1 - x0, dy10 = y1 - y0;
  if (pbc[0]) {
    dx21 -= L[0] * round(dx21 / L[0]);
    dx10 -= L[0] * round(dx10 / L[0]);
  }
  if (pbc[1]) {
    dy21 -= L[1] * round(dy21 / L[1]);
    dy10 -= L[1] * round(dy10 / L[1]);
  }

  double l2 = pow(dx21, 2) + pow(dy21, 2);  // |(pt2 - pt1)|^2
  if (l2 == 0.0)                            // pt2 == pt1 case
    return sqrt(pow(dx10, 2) + pow(dy10, 2));

  double dot = (-dx10) * (dx21) +
               (-dy10) * (dy21);  // (pt0 - pt1) dot (pt2 - pt1)
  const double t = max(0.0, min(1.0, dot / l2));
  xcomp = x0 - x1 - t * (dx21);
  ycomp = y0 - y1 - t * (dy21);
  const double distance = sqrt(xcomp * xcomp + ycomp * ycomp);
  return distance;
}

double dpm::linePointDistancesAndProjection(double x1, double y1, double x2, double y2, double x0, double y0, double& xcomp, double& ycomp, double& projection) {
  // get the distance from a line SEGMENT (ending at its endpoints) going through (x1,y1), (x2,y2) and a
  // point located at (x0,y0), and extract x and y components of the distance
  // and determine whether the projected contact is vertex-vertex between pt 0 and pt 1: (projection < 0)
  //  or vertex-line-segment betw 0 and 1-2 (projection 0<x<1)
  //  or vertex-vertex betw 0 and 2 (projection >1)
  double dx21 = x2 - x1, dy21 = y2 - y1;

  double dx10 = x1 - x0;
  double dy10 = y1 - y0;
  if (pbc[0]) {
    dx21 -= L[0] * round(dx21 / L[0]);
    dx10 -= L[0] * round(dx10 / L[0]);
  }
  if (pbc[1]) {
    dy21 -= L[1] * round(dy21 / L[1]);
    dy10 -= L[1] * round(dy10 / L[1]);
  }

  double l2_sq = dx21 * dx21 + dy21 * dy21;  // |(pt2 - pt1)|^2
  if (l2_sq == 0.0) {                        // pt2 == pt1 case
    xcomp = dx10;
    ycomp = dy10;
    return sqrt(dx10 * dx10 + dy10 * dy10);
  }

  double dot = (-dx10) * (dx21) +
               (-dy10) * (dy21);  // (pt0 - pt1) dot (pt2 - pt1)
  projection = dot / l2_sq;

  const double t = max(0.0, min(1.0, projection));  // t is restricted to [0,1], which parametrizes the line segment (v + t (w - v))

  xcomp = x0 - x1 - t * (dx21);  // x0 - projection in x
  if (pbc[0]) {
    xcomp -= L[0] * round(xcomp / L[0]);
  }
  ycomp = y0 - y1 - t * (dy21);  // y0 - projection in y
  if (pbc[1]) {
    ycomp -= L[1] * round(ycomp / L[1]);
  }

  const double distance = sqrt(xcomp * xcomp + ycomp * ycomp);
  return distance;
}

void dpm::generateCircularBoundary(int numEdges, double radius, double cx, double cy, std::vector<double>& poly_x, std::vector<double>& poly_y) {
  // place a circular boundary, then fix L[0] and L[1] just outside the circle for plotting purposes.
  // in initializePositions2D, isCircle flag is chosen to ensure that particles fall within L/2 (radius) of L/2 (center of box)
  // cx cy are the center of the boundary polygon
  cout << "in generateCircularBoundary\n";
  cout << L[0] << '\t' << L[1] << '\n';
  generateCircle(numEdges, cx, cy, radius, poly_x, poly_y);

  cout << "after generating the polygon boundary, set L[0] and L[1] to be outside this boundary\n";
  L[0] = 1.0 * *max_element(std::begin(poly_x), std::end(poly_x));
  L[1] = 1.0 * *max_element(std::begin(poly_y), std::end(poly_y));
  cout << L[0] << '\t' << L[1] << '\n';

  // plot final polygonal boundary, make sure to clear the file when running a new simulation (should be run with ofstream polyBoundary.txt without app before generateCircularBoundary is called)
  ofstream boundaryStream("polyBoundary.txt", std::ios_base::app);
  for (int i = 0; i < poly_x.size(); i++) {
    boundaryStream << poly_x[i] << '\t' << poly_y[i] << '\t';
  }
  boundaryStream << '\n';
}

void dpm::generateCircle(int numEdges, double cx, double cy, double r, std::vector<double>& px, std::vector<double>& py) {
  // generate an n-gon at center cx,cy with radius r. Modifies data in px, py to give the n-gon
  double theta;
  std::vector<double> poly_x_temp, poly_y_temp;
  cout << "px, py (list) = \n";
  for (int i = 0; i < numEdges; i++) {
    theta = i * 2 * PI / numEdges;
    poly_x_temp.push_back(r * cos(theta) + cx);
    poly_y_temp.push_back(r * sin(theta) + cy);
    cout << r * cos(theta) + cx << '\t' << r * sin(theta) + cy << '\n';
  }
  px = poly_x_temp;
  py = poly_y_temp;
}

void dpm::generateRectangularBoundary(double radius, double cx, double cy, std::vector<double>& poly_x, std::vector<double>& poly_y) {
  // place a rectangular boundary, then fix L[0] and L[1] just outside the circle for plotting purposes.
  // in initializePositions2D, isCircle flag is chosen to ensure that particles fall within L/2 (radius) of L/2 (center of box)
  // cx cy are the center of the boundary polygon
  cout << "in generateRectangularBoundary\n";
  cout << L[0] << '\t' << L[1] << '\n';
  // generateCircle(numEdges, cx, cy, radius, poly_x, poly_y);
  double aspectRatio = 2;
  std::vector<double> x_pos = {cx + radius, cx - radius, cx - radius, cx + radius},
                      y_pos = {cy + radius * aspectRatio, cy + radius * aspectRatio, cy - radius * aspectRatio, cy - radius * aspectRatio};
  for (int i = 0; i < 4; i++) {
    poly_x.push_back(x_pos[i]);
    poly_y.push_back(y_pos[i]);
  }

  cout << "after generating the polygon boundary, set L[0] and L[1] to be outside this boundary\n";
  L[0] = 1.05 * *max_element(std::begin(poly_x), std::end(poly_x));
  L[1] = 1.05 * *max_element(std::begin(poly_y), std::end(poly_y));
  cout << L[0] << '\t' << L[1] << '\n';

  // plot final polygonal boundary, make sure to clear the file when running a new simulation (should be run with ofstream polyBoundary.txt without app before generateCircularBoundary is called)
  ofstream boundaryStream("polyBoundary.txt", std::ios_base::app);
  for (int i = 0; i < poly_x.size(); i++) {
    boundaryStream << poly_x[i] << '\t' << poly_y[i] << '\t';
  }
  boundaryStream << '\n';
}

void dpm::generateHorseshoeBoundary(double cx, double cy, std::vector<double>& poly_x, std::vector<double>& poly_y) {
  // place a rectangular boundary, then fix L[0] and L[1] just outside the circle for plotting purposes.
  // in initializePositions2D, isCircle flag is chosen to ensure that particles fall within L/2 (radius) of L/2 (center of box)
  // cx cy are the center of the boundary polygon
  // generateCircle(numEdges, cx, cy, radius, poly_x, poly_y);

  /*std::vector<double> x_pos = {cx + radius, cx - radius, cx - radius, cx + radius},
                      y_pos = {cy + radius * aspectRatio, cy + radius * aspectRatio, cy - radius * aspectRatio, cy - radius * aspectRatio};
  for (int i = 0; i < 4; i++) {
    poly_x.push_back(x_pos[i]);
    poly_y.push_back(y_pos[i]);
  }*/
  double radius = 10;
  double length = 2 * radius / 3;
  std::vector<double> x_pos = {0, 0, length, length, 2 * length, 2 * length, 3 * length, 3 * length},
                      y_pos = {radius, 3 * radius, 3 * radius, radius, radius, 3 * radius, 3 * radius, radius};
  // calculate coordinates of a semicircle with center at (1.5 * length, radius))
  int numEdges = 30;
  for (int n = 0; n < numEdges; n++) {
    double theta = -n * PI / (numEdges - 1);
    x_pos.push_back(1.5 * length + radius * cos(theta));
    y_pos.push_back(radius + radius * sin(theta));
  }

  // wrote the above in clockwise order, switching it to be counter clockwise
  std::reverse(x_pos.begin(), x_pos.end());
  std::reverse(y_pos.begin(), y_pos.end());

  for (int i = 0; i < x_pos.size(); i++) {
    cout << x_pos[i] << '\t' << y_pos[i] << '\n';
    poly_x.push_back(x_pos[i]);
    poly_y.push_back(y_pos[i]);
  }

  cout << "after generating the polygon boundary, set L[0] and L[1] to be outside this boundary\n";
  L[0] = 1.05 * *max_element(std::begin(poly_x), std::end(poly_x));
  L[1] = 1.05 * *max_element(std::begin(poly_y), std::end(poly_y));
  cout << L[0] << '\t' << L[1] << '\n';

  // plot final polygonal boundary, make sure to clear the file when running a new simulation (should be run with ofstream polyBoundary.txt without app before generateCircularBoundary is called)
  ofstream boundaryStream("polyBoundary.txt", std::ios_base::app);
  for (int i = 0; i < poly_x.size(); i++) {
    boundaryStream << poly_x[i] << '\t' << poly_y[i] << '\t';
  }
  boundaryStream << '\n';
}

void dpm::replaceCircularBoundary(int polyShapeID, double aspectRatio) {
  int numBoundaryElements = poly_bd_x.size();
  // calculate bounding box for the poly_bd
  if (polyShapeID == 1) {
    double top, bottom, left, right;
    top = 1.0 * *max_element(std::begin(poly_bd_y[0]), std::end(poly_bd_y[0])) + r[0];
    bottom = 1.0 * *min_element(std::begin(poly_bd_y[0]), std::end(poly_bd_y[0])) - r[0];
    left = 1.0 * *min_element(std::begin(poly_bd_x[0]), std::end(poly_bd_x[0])) - r[0];
    right = 1.0 * *max_element(std::begin(poly_bd_x[0]), std::end(poly_bd_x[0])) + r[0];
    double xshift = 0.0;
    double yshift = (top - bottom) / 2;
    double height = top - bottom;
    // extend the rectangle only in the y-direction. this leaves the cluster of cells in the bottom of the rectangle
    top += (aspectRatio - 1.0) * height;
    std::vector<double> x_pos = {right, left, left, right},
                        y_pos = {top, top, bottom, bottom};
    // clear existing poly_bd
    poly_bd_x[0].clear();
    poly_bd_y[0].clear();
    // assuming only one tissue, so use poly_bd[0]
    for (int i = 0; i < 4; i++) {
      poly_bd_x[0].push_back(x_pos[i]);
      poly_bd_y[0].push_back(y_pos[i]);
    }

    // move cells to the center of the rectangle
    moveSimulationToPositiveCoordinates(0.0, (top - bottom) / 2 - yshift);

    cout << "after generating the polygon boundary, set L[0] and L[1] to be outside this boundary\n";
    L[0] = right * 1.05;
    L[1] = top * 1.05;
    cout << L[0] << '\t' << L[1] << '\n';
  } else
    assert(false);
}

std::vector<double> dpm::resample_polygon(std::vector<double> px, std::vector<double> py, double perimeter, int numPoints) {
  cout << "resample_polygon with perimeter " << perimeter << "\t, numPoints = " << numPoints << "\t, px.size() = " << px.size() << '\n';
  double arc_length = perimeter / numPoints, t = 0.0, d_t = 0.0;
  double newx, newy, lerp;
  vector<double> result = {px[0], py[0]};

  for (int vi = 0; vi < px.size(); vi++) {
    d_t = sqrt(pow(px[vi] - px[(vi + 1) % px.size()], 2) + pow(py[vi] - py[(vi + 1) % py.size()], 2)) / arc_length;
    while (t + d_t >= result.size() / 2 && result.size() / 2 < numPoints) {
      lerp = (result.size() / 2 - t) / d_t;
      newx = (1 - lerp) * px[vi] + lerp * px[(vi + 1) % px.size()];
      newy = (1 - lerp) * py[vi] + lerp * py[(vi + 1) % py.size()];
      result.push_back(newx);
      result.push_back(newy);
      cout << "newx newy = " << newx << '\t' << newy << '\n';
    }
    t += d_t;
  }
  cout << "arc_length = " << arc_length << '\n';
  cout << "in resample_polygon, result.size() = " << result.size() << '\n';
  return result;
}

// Check if a point (x, y) is inside a polygon defined by vectors poly_x and poly_y
bool dpm::isInsidePolygon(double x, double y, const std::vector<double>& poly_x, const std::vector<double>& poly_y) {
  int i, j, nvert = poly_x.size();
  bool inside = false;

  for (i = 0, j = nvert - 1; i < nvert; j = i++) {
    if (((poly_y[i] > y) != (poly_y[j] > y)) &&
        (x < (poly_x[j] - poly_x[i]) * (y - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i])) {
      inside = !inside;
    }
  }

  return inside;
}

/******************************

        D P M  F O R C E

                        U P D A T E S

*******************************/

void dpm::resetForcesAndEnergy() {
  fill(F.begin(), F.end(), 0.0);
  fill(stress.begin(), stress.end(), 0.0);
  fill(fieldStress.begin(), fieldStress.end(), vector<double>(3));
  fill(fieldShapeStress.begin(), fieldShapeStress.end(), vector<double>(3));
  fill(fieldStressCells.begin(), fieldStressCells.end(), vector<double>(3));
  fill(fieldShapeStressCells.begin(), fieldShapeStressCells.end(), vector<double>(3));
  U = 0.0;
  fill(cellU.begin(), cellU.end(), 0.0);
}

void dpm::shapeForces2D() {
  // local variables
  int ci, gi, vi, nvtmp;
  double fa, fli, flim1, fb, cx, cy, xi, yi;
  double rho0, l0im1, l0i, a0tmp, atmp;
  double dx, dy, da, dli, dlim1, dtim1, dti, dtip1;
  double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
  double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
  double nim1x, nim1y, nix, niy, sinim1, sini, sinip1, cosim1, cosi, cosip1;
  double ddtim1, ddti;
  double forceX, forceY;
  double unwrappedX, unwrappedY;
  std::vector<double> li_vector(NVTOT, 0.0);

  // loop over vertices, add to force
  rho0 = sqrt(a0.at(0));
  ci = 0;
  for (gi = 0; gi < NVTOT; gi++) {
    // -- Area force (and get cell index ci)
    if (ci < NCELLS) {
      if (gi == szList[ci]) {
        // shape information
        nvtmp = nv[ci];
        a0tmp = a0[ci];

        // preferred segment length of last segment
        l0im1 = l0[im1[gi]];

        // compute area deviation
        atmp = area(ci);
        da = (atmp / a0tmp) - 1.0;

        // update potential energy
        U += 0.5 * ka * (da * da);
        cellU[ci] += 0.5 * ka * (da * da);

        // shape force parameters
        fa = ka * da * (rho0 / a0tmp);
        fb = kb * rho0;

        // compute cell center of mass
        xi = x[NDIM * gi];
        yi = x[NDIM * gi + 1];
        cx = xi;
        cy = yi;
        for (vi = 1; vi < nvtmp; vi++) {
          // get distances between vim1 and vi
          dx = x[NDIM * (gi + vi)] - xi;
          dy = x[NDIM * (gi + vi) + 1] - yi;
          if (pbc[0])
            dx -= L[0] * round(dx / L[0]);
          if (pbc[1])
            dy -= L[1] * round(dy / L[1]);

          // add to centers
          xi += dx;
          yi += dy;

          cx += xi;
          cy += yi;
        }
        cx /= nvtmp;
        cy /= nvtmp;

        // get coordinates relative to center of mass
        rix = x[NDIM * gi] - cx;
        riy = x[NDIM * gi + 1] - cy;

        // get prior adjacent vertices
        rim2x = x[NDIM * im1[im1[gi]]] - cx;
        rim2y = x[NDIM * im1[im1[gi]] + 1] - cy;
        if (pbc[0])
          rim2x -= L[0] * round(rim2x / L[0]);
        if (pbc[1])
          rim2y -= L[1] * round(rim2y / L[1]);

        rim1x = x[NDIM * im1[gi]] - cx;
        rim1y = x[NDIM * im1[gi] + 1] - cy;
        if (pbc[0])
          rim1x -= L[0] * round(rim1x / L[0]);
        if (pbc[1])
          rim1y -= L[1] * round(rim1y / L[1]);

        // get prior segment vectors
        lim2x = rim1x - rim2x;
        lim2y = rim1y - rim2y;

        lim1x = rix - rim1x;
        lim1y = riy - rim1y;

        // increment cell index
        ci++;
      }
    }

    if (ci >= NCELLS)
      ci--;

    // unwrapped vertex coordinate
    unwrappedX = cx + rix;
    unwrappedY = cy + riy;

    // preferred segment length
    l0i = l0[gi];

    // get next adjacent vertices
    rip1x = x[NDIM * ip1[gi]] - cx;
    rip1y = x[NDIM * ip1[gi] + 1] - cy;
    if (pbc[0])
      rip1x -= L[0] * round(rip1x / L[0]);
    if (pbc[1])
      rip1y -= L[1] * round(rip1y / L[1]);

    // -- Area force (comes from a cross product)
    forceX = 0.5 * fa * (rim1y - rip1y);
    forceY = 0.5 * fa * (rip1x - rim1x);

    F[NDIM * gi] += forceX;
    F[NDIM * gi + 1] += forceY;

    fieldShapeStress[gi][0] += unwrappedX * forceX;
    fieldShapeStress[gi][1] += unwrappedY * forceY;
    fieldShapeStress[gi][2] += unwrappedX * forceY;

    // -- Perimeter force
    lix = rip1x - rix;
    liy = rip1y - riy;

    // segment lengths
    lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);
    li = sqrt(lix * lix + liy * liy);
    li_vector[gi] = li;  // save li for Maxwell relaxation integration

    // segment deviations (note: m is prior vertex, p is next vertex i.e. gi - 1, gi + 1 mod the right number of vertices)
    dlim1 = (lim1 / l0im1) - 1.0;
    dli = (li / l0i) - 1.0;

    // segment forces
    flim1 = kl * (rho0 / l0im1);
    fli = kl * (rho0 / l0i);

    // add to forces
    forceX = (fli * dli * lix / li) - (flim1 * dlim1 * lim1x / lim1);
    forceY = (fli * dli * liy / li) - (flim1 * dlim1 * lim1y / lim1);
    F[NDIM * gi] += forceX;
    F[NDIM * gi + 1] += forceY;

    // note - Andrew here, confirmed that the shape stress matrix is diagonal as written
    fieldShapeStress[gi][0] += unwrappedX * forceX;
    fieldShapeStress[gi][1] += unwrappedY * forceY;
    fieldShapeStress[gi][2] += unwrappedX * forceY;

    // update potential energy
    U += 0.5 * kl * (dli * dli);
    cellU[ci] += 0.5 * kl * (dli * dli);

    // -- Bending force
    if (kb > 0) {
      // get ip2 for third angle
      rip2x = x[NDIM * ip1[ip1[gi]]] - cx;
      rip2y = x[NDIM * ip1[ip1[gi]] + 1] - cy;
      if (pbc[0])
        rip2x -= L[0] * round(rip2x / L[0]);
      if (pbc[1])
        rip2y -= L[1] * round(rip2y / L[1]);

      // get last segment length
      lip1x = rip2x - rip1x;
      lip1y = rip2y - rip1y;

      // get angles
      sinim1 = lim1x * lim2y - lim1y * lim2x;
      cosim1 = lim1x * lim2x + lim1y * lim2y;

      sini = lix * lim1y - liy * lim1x;
      cosi = lix * lim1x + liy * lim1y;

      sinip1 = lip1x * liy - lip1y * lix;
      cosip1 = lip1x * lix + lip1y * liy;

      // get normal vectors
      nim1x = lim1y;
      nim1y = -lim1x;

      nix = liy;
      niy = -lix;

      // get change in angles
      dtim1 = atan2(sinim1, cosim1) - t0[im1[gi]];
      dti = atan2(sini, cosi) - t0[gi];
      dtip1 = atan2(sinip1, cosip1) - t0[ip1[gi]];

      // get delta delta theta's
      ddtim1 = (dti - dtim1) / (lim1 * lim1);
      ddti = (dti - dtip1) / (li * li);

      // add to force
      F[NDIM * gi] += fb * (ddtim1 * nim1x + ddti * nix);
      F[NDIM * gi + 1] += fb * (ddtim1 * nim1y + ddti * niy);

      // update potential energy
      U += 0.5 * kb * (dti * dti);
      cellU[ci] += 0.5 * kb * (dti * dti);
    }

    // update old coordinates
    rim2x = rim1x;
    rim1x = rix;
    rix = rip1x;

    rim2y = rim1y;
    rim1y = riy;
    riy = rip1y;

    // update old segment vectors
    lim2x = lim1x;
    lim2y = lim1y;

    lim1x = lix;
    lim1y = liy;

    // update old preferred segment length
    l0im1 = l0i;
  }

  maxwellRelaxationRestLengths(li_vector);

  // normalize per-cell stress by preferred cell area
  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv[ci]; vi++) {
      int gi = gindex(ci, vi);
      fieldShapeStressCells[ci][0] += fieldShapeStress[gi][0];
      fieldShapeStressCells[ci][1] += fieldShapeStress[gi][1];
      fieldShapeStressCells[ci][2] += fieldShapeStress[gi][2];
    }
    // nondimensionalize the stress
    fieldShapeStressCells[ci][0] *= rho0 / a0[ci];
    fieldShapeStressCells[ci][1] *= rho0 / a0[ci];
    fieldShapeStressCells[ci][2] *= rho0 / a0[ci];
  }
}

void dpm::maxwellRelaxationRestLengths(std::vector<double>& l) {
  // we are integrating a 1D equation of motion for rest lengths
  //  assuming a Maxwell model for stress relaxation.
  //  k(l-l0) is the spring, -B*vl0 is the dashpot.
  double al0, al0_old;
  double li;
  for (int i = 0; i < l0.size(); i++) {
    li = l[i];

    // overdamped integration of rest length relaxation
    vl0[i] = kl / maxwellRelaxationTime * (li - l0[i]);
    vl0[i] += kl / taus * (l00[i] - li);
    l0[i] += vl0[i] * dt;
  }
}

void dpm::vertexRepulsiveForces2D() {
  // local variables
  int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
  double sij, rij, dx, dy, rho0;
  double ftmp, fx, fy;

  // sort particles
  sortNeighborLinkedList2D();

  // get fundamental length
  rho0 = sqrt(a0.at(0));

  // reset contact network
  fill(cij.begin(), cij.end(), 0);

  // loop over boxes in neighbor linked list
  for (bi = 0; bi < NBX; bi++) {
    // get start of list of vertices
    pi = head[bi];

    // loop over linked list
    while (pi > 0) {
      // real particle index
      gi = pi - 1;

      // next particle in list
      pj = list[pi];

      // loop down neighbors of pi in same cell
      while (pj > 0) {
        // real index of pj
        gj = pj - 1;

        if (gj == ip1[gi] || gj == im1[gi]) {
          pj = list[pj];
          continue;
        }

        // contact distance
        sij = r[gi] + r[gj];

        // particle distance
        dx = x[NDIM * gj] - x[NDIM * gi];
        if (pbc[0])
          dx -= L[0] * round(dx / L[0]);
        if (dx < sij) {
          dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
          if (pbc[1])
            dy -= L[1] * round(dy / L[1]);
          if (dy < sij) {
            rij = sqrt(dx * dx + dy * dy);
            if (rij < sij) {
              // force scale
              ftmp = kc * (1 - (rij / sij)) * (rho0 / sij);
              fx = ftmp * (dx / rij);
              fy = ftmp * (dy / rij);

              // add to forces
              F[NDIM * gi] -= fx;
              F[NDIM * gi + 1] -= fy;

              F[NDIM * gj] += fx;
              F[NDIM * gj + 1] += fy;

              // increase potential energy
              U += 0.5 * kc * pow((1 - (rij / sij)), 2.0);

              // add to virial stress
              stress[0] += dx * fx;
              stress[1] += dy * fy;
              stress[2] += 0.5 * (dx * fy + dy * fx);

              fieldStress[gi][0] += dx * fx;
              fieldStress[gi][1] += dy * fy;
              fieldStress[gi][2] += 0.5 * (dx * fy + dy * fx);

              // add to contacts
              cindices(ci, vi, gi);
              cindices(cj, vj, gj);
              cellU[ci] += 0.5 * kc * pow((1 - (rij / sij)), 2.0) / 2.0;
              cellU[cj] += 0.5 * kc * pow((1 - (rij / sij)), 2.0) / 2.0;

              if (ci > cj)
                cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
              else if (ci < cj)
                cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
            }
          }
        }

        // update pj
        pj = list[pj];
      }

      // test overlaps with forward neighboring cells
      for (bj = 0; bj < NNN; bj++) {
        // only check if boundaries permit
        if (nn[bi][bj] == -1)
          continue;

        // get first particle in neighboring cell
        pj = head[nn[bi][bj]];

        // loop down neighbors of pi in same cell
        while (pj > 0) {
          // real index of pj
          gj = pj - 1;

          if (gj == ip1[gi] || gj == im1[gi]) {
            pj = list[pj];
            continue;
          }
          // contact distance
          sij = r[gi] + r[gj];

          // particle distance
          dx = x[NDIM * gj] - x[NDIM * gi];
          if (pbc[0])
            dx -= L[0] * round(dx / L[0]);
          if (dx < sij) {
            dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
            if (pbc[1])
              dy -= L[1] * round(dy / L[1]);
            if (dy < sij) {
              rij = sqrt(dx * dx + dy * dy);
              if (rij < sij) {
                // force scale
                ftmp = kc * (1 - (rij / sij)) * (rho0 / sij);
                fx = ftmp * (dx / rij);
                fy = ftmp * (dy / rij);

                // add to forces
                F[NDIM * gi] -= fx;
                F[NDIM * gi + 1] -= fy;

                F[NDIM * gj] += fx;
                F[NDIM * gj + 1] += fy;

                // increase potential energy
                U += 0.5 * kc * pow((1 - (rij / sij)), 2.0);

                // add to virial stress
                stress[0] += dx * fx;
                stress[1] += dy * fy;
                stress[2] += 0.5 * (dx * fy + dy * fx);

                fieldStress[gi][0] += dx * fx;
                fieldStress[gi][1] += dy * fy;
                fieldStress[gi][2] += 0.5 * (dx * fy + dy * fx);

                // add to contacts
                cindices(ci, vi, gi);
                cindices(cj, vj, gj);

                cellU[ci] += 0.5 * kc * pow((1 - (rij / sij)), 2.0) / 2.0;
                cellU[cj] += 0.5 * kc * pow((1 - (rij / sij)), 2.0) / 2.0;

                if (ci > cj)
                  cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                else if (ci < cj)
                  cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
              }
            }
          }

          // update pj
          pj = list[pj];
        }
      }

      // update pi index to be next
      pi = list[pi];
    }
  }

  // normalize stress by box area, make dimensionless
  // units explanation (Andrew): pressure is force/area. up til now, stress has been force times distance.
  // so after dividing by area, we have force / distance. We need to multiply one last time by rho0, done.

  stress[0] *= (rho0 / (L[0] * L[1]));
  stress[1] *= (rho0 / (L[0] * L[1]));
  stress[2] *= (rho0 / (L[0] * L[1]));

  // normalize per-cell stress by preferred cell area
  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv[ci]; vi++) {
      int gi = gindex(ci, vi);
      fieldStressCells[ci][0] += fieldStress[gi][0];
      fieldStressCells[ci][1] += fieldStress[gi][1];
      fieldStressCells[ci][2] += fieldStress[gi][2];
    }
    fieldStressCells[ci][0] *= rho0 / a0[ci];
    fieldStressCells[ci][1] *= rho0 / a0[ci];
    fieldStressCells[ci][2] *= rho0 / a0[ci];
  }
}

void dpm::vertexAttractiveForces2D() {
  // local variables
  int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
  double sij, rij, dx, dy, rho0;
  double ftmp, fx, fy;

  // attraction shell parameters
  double shellij, cutij, xij, kint = (kc * l1) / (l2 - l1);

  // sort particles
  sortNeighborLinkedList2D();

  // get fundamental length
  rho0 = sqrt(a0[0]);

  // reset contact network
  fill(cij.begin(), cij.end(), 0);

  // loop over boxes in neighbor linked list
  for (bi = 0; bi < NBX; bi++) {
    // get start of list of vertices
    pi = head[bi];

    // loop over linked list
    while (pi > 0) {
      // real particle index
      gi = pi - 1;

      // cell index of gi
      cindices(ci, vi, gi);

      // next particle in list
      pj = list[pi];

      // loop down neighbors of pi in same cell
      while (pj > 0) {
        // real index of pj
        gj = pj - 1;

        cindices(cj, vj, gj);

        if (gj == ip1[gi] || gj == im1[gi]) {
          pj = list[pj];
          continue;
        }

        // contact distance
        sij = r[gi] + r[gj];

        // attraction distances
        shellij = (1.0 + l2) * sij;
        cutij = (1.0 + l1) * sij;

        // particle distance
        dx = x[NDIM * gj] - x[NDIM * gi];
        if (pbc[0])
          dx -= L[0] * round(dx / L[0]);
        if (dx < shellij) {
          dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
          if (pbc[1])
            dy -= L[1] * round(dy / L[1]);
          if (dy < shellij) {
            rij = sqrt(dx * dx + dy * dy);
            if (rij < shellij) {
              // scaled distance
              xij = rij / sij;

              // pick force based on vertex-vertex distance
              if (rij > cutij) {
                // force scale
                ftmp = kint * (xij - 1.0 - l2) / sij;

                // increase potential energy
                U += -0.5 * kint * pow(1.0 + l2 - xij, 2.0);
                cellU[ci] += -0.5 * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                cellU[cj] += -0.5 * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
              } else {
                // force scale
                ftmp = kc * (1 - xij) / sij;

                // increase potential energy
                U += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                cellU[ci] += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                cellU[cj] += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
              }

              // force elements
              fx = ftmp * (dx / rij);
              fy = ftmp * (dy / rij);

              // add to forces
              F[NDIM * gi] -= fx;
              F[NDIM * gi + 1] -= fy;

              F[NDIM * gj] += fx;
              F[NDIM * gj + 1] += fy;

              // add to virial stress
              stress[0] += dx * fx;
              stress[1] += dy * fy;
              stress[2] += 0.5 * (dx * fy + dy * fx);

              fieldStress[gi][0] += dx * fx;
              fieldStress[gi][1] += dy * fy;
              fieldStress[gi][2] += 0.5 * (dx * fy + dy * fx);

              // cindices(cj, vj, gj);
              //  add to contacts
              if (ci > cj)
                cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
              else if (ci < cj)
                cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
            }
          }
        }

        // update pj
        pj = list[pj];
      }

      // test overlaps with forward neighboring cells
      for (bj = 0; bj < NNN; bj++) {
        // only check if boundaries permit
        if (nn[bi][bj] == -1)
          continue;

        // get first particle in neighboring cell
        pj = head[nn[bi][bj]];

        // loop down neighbors of pi in same cell
        while (pj > 0) {
          // real index of pj
          gj = pj - 1;

          cindices(cj, vj, gj);

          if (gj == ip1[gi] || gj == im1[gi]) {
            pj = list[pj];
            continue;
          }

          // contact distance
          sij = r[gi] + r[gj];

          // attraction distances
          shellij = (1.0 + l2) * sij;
          cutij = (1.0 + l1) * sij;

          // particle distance
          dx = x[NDIM * gj] - x[NDIM * gi];
          if (pbc[0])
            dx -= L[0] * round(dx / L[0]);
          if (dx < shellij) {
            dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
            if (pbc[1])
              dy -= L[1] * round(dy / L[1]);
            if (dy < shellij) {
              rij = sqrt(dx * dx + dy * dy);
              if (rij < shellij) {
                // scaled distance
                xij = rij / sij;

                // pick force based on vertex-vertex distance
                if (rij > cutij) {
                  // force scale
                  ftmp = kint * (xij - 1.0 - l2) / sij;

                  // increase potential energy
                  U += -0.5 * kint * pow(1.0 + l2 - xij, 2.0);
                  cellU[ci] += -0.5 * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                  cellU[cj] += -0.5 * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                } else {
                  // force scale
                  ftmp = kc * (1 - xij) / sij;

                  // increase potential energy
                  U += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                  cellU[ci] += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                  cellU[cj] += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                }

                // force elements
                fx = ftmp * (dx / rij);
                fy = ftmp * (dy / rij);

                // add to forces
                F[NDIM * gi] -= fx;
                F[NDIM * gi + 1] -= fy;

                F[NDIM * gj] += fx;
                F[NDIM * gj + 1] += fy;

                // add to virial stress
                stress[0] += dx * fx;
                stress[1] += dy * fy;
                stress[2] += 0.5 * (dx * fy + dy * fx);

                fieldStress[gi][0] += dx * fx;
                fieldStress[gi][1] += dy * fy;
                fieldStress[gi][2] += 0.5 * (dx * fy + dy * fx);

                if (ci > cj)
                  cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                else if (ci < cj)
                  cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
              }
            }
          }

          // update pj
          pj = list[pj];
        }
      }

      // update pi index to be next
      pi = list[pi];
    }
  }

  // normalize stress by box area, make dimensionless
  stress[0] *= (rho0 / (L[0] * L[1]));
  stress[1] *= (rho0 / (L[0] * L[1]));
  stress[2] *= (rho0 / (L[0] * L[1]));

  // normalize per-cell stress by cell area
  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv[ci]; vi++) {
      int gi = gindex(ci, vi);
      fieldStressCells[ci][0] += fieldStress[gi][0];
      fieldStressCells[ci][1] += fieldStress[gi][1];
      fieldStressCells[ci][2] += fieldStress[gi][2];
    }
    fieldStressCells[ci][0] *= rho0 / a0[ci];
    fieldStressCells[ci][1] *= rho0 / a0[ci];
    fieldStressCells[ci][2] *= rho0 / a0[ci];
  }
}

// if there are multiple polygonal bds, need to throw evaluatePolygonalWallForces in a loop over poly_bd_x.size()
void dpm::evaluatePolygonalWallForces(const std::vector<double>& poly_x, const std::vector<double>& poly_y, bool attractionOn) {
  // evaluates particle-wall forces for a polygonal boundary specified by poly_x,poly_y. Does not compute stress yet.
  int n = poly_x.size();
  double wallThicknessBuffer = 1.0;  // used to make sure wall can't phase through a vertex when scaling vertex or wall positions
  double distanceParticleWall, scaledDistParticleWall, Rx, Ry, dw;
  double kint = (kc * l1) / (l2 - l1);
  double bound_x1, bound_x2, bound_y1, bound_y2;
  double shellij, cutij, ftmp;
  // loop over boundary bars
  // loop over particles
  //  compute particle-boundary bar overlaps
  //  if overlap, Fx += K * dw * Rx/R, where K is a constant, dw = diameter/2 - R, Rx = x - px, R = sqrt(Rx^2 + Ry^2)
  for (int bound_i = 0; bound_i < n; bound_i++) {
    // use distanceLineAndPoint to get R, Rx, and Ry
    bound_x1 = poly_x[bound_i];
    bound_x2 = poly_x[(bound_i + 1) % n];
    bound_y1 = poly_y[bound_i];
    bound_y2 = poly_y[(bound_i + 1) % n];
    for (int i = 0; i < NVTOT; i++) {
      distanceParticleWall = distanceLinePointComponents(bound_x1, bound_y1, bound_x2, bound_y2, x[i * NDIM], x[i * NDIM + 1], Rx, Ry);
      dw = fabs(r[i] - distanceParticleWall);
      if (attractionOn) {  // attractive, so choose a weaker K and check for attractive shell interaction
        scaledDistParticleWall = distanceParticleWall / r[i];
        // K = 1;
        shellij = (1.0 + l2) * r[i];
        cutij = (1.0 + l1) * r[i];
        if (distanceParticleWall <= shellij) {  // within attracting shell 2
          if (distanceParticleWall > cutij) {
            ftmp = kint * (scaledDistParticleWall - 1.0 - l2) / r[i];
          } else {  // within attracting shell 1, potentially within repulsion distance
            ftmp = kc * (1 - scaledDistParticleWall) / r[i];
          }
          F[i * NDIM] += ftmp * Rx / distanceParticleWall;
          F[i * NDIM + 1] += ftmp * Ry / distanceParticleWall;
        }
      } else if (distanceParticleWall <= wallThicknessBuffer * r[i]) {  // purely repulsive, only look for particle-wall overlap, with thickness buffer on wall
        F[i * NDIM] += kc * dw * Rx / distanceParticleWall;
        F[i * NDIM + 1] += kc * dw * Ry / distanceParticleWall;
        U += kc / 2 * pow(dw, 2);
      }
    }
  }
}

void dpm::repulsiveForceUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexRepulsiveForces2D();
}

void dpm::attractiveForceUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D();
}

/******************************

        D P M

                I N T E G R A T O R S

*******************************/

void dpm::setdt(double dt0) {
  // local variables
  int i = 0;
  double ta = 0, tl = 0, tb = 0, tmin = 0, rho0 = 0;

  // typical length
  rho0 = sqrt(a0.at(0));

  // set typical time scales
  ta = rho0 / sqrt(ka);
  tl = (rho0 * l0.at(0)) / sqrt(ka * kl);
  tb = (rho0 * l0.at(0)) / sqrt(ka * kb);

  // set main time scale as min
  tmin = 1e8;
  if (ta < tmin)
    tmin = ta;
  if (tl < tmin)
    tmin = tl;
  if (tb < tmin)
    tmin = tb;

  // set dt
  dt = dt0 * tmin;
}

void dpm::vertexFIRE2D(dpmMemFn forceCall, double Ftol, double dt0) {
  // local variables
  int i;
  double rho0;

  // check to see if cell linked-list has been initialized
  if (NBX == -1) {
    cerr << "	** ERROR: In dpm::fire, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
    exit(1);
  }

  // FIRE variables
  double P, fnorm, fcheck, vnorm, alpha, dtmax, dtmin;
  int npPos, npNeg, fireit;

  // set dt based on geometric parameters
  // NOTE: if calling this in a different program, know that it will modify dt! make sure to reset dt if not intended.
  setdt(dt0);

  // Initialize FIRE variables
  P = 0;
  fnorm = 0;
  vnorm = 0;
  alpha = alpha0;

  dtmax = 10.0 * dt;
  dtmin = 1e-2 * dt;

  npPos = 0;
  npNeg = 0;

  fireit = 0;
  fcheck = 10 * Ftol;

  // reset forces and velocities
  resetForcesAndEnergy();
  fill(v.begin(), v.end(), 0.0);

  // length scale
  rho0 = sqrt(a0.at(0));

  // relax forces using FIRE
  while (fcheck > Ftol && fireit < itmax) {
    // compute P
    P = 0.0;
    for (i = 0; i < vertDOF; i++) {
      P += v[i] * F[i];
    }

    // print to console
    if (fireit % NSKIP == 0) {
      cout << endl
           << endl;
      cout << "===========================================" << endl;
      cout << " 	F I R E 						" << endl;
      cout << "		M I N I M I Z A T I O N 	" << endl;
      cout << "===========================================" << endl;
      cout << endl;
      cout << "	** fireit 	= " << fireit << " / " << itmax << endl;
      cout << "	** fcheck 	= " << fcheck << endl;
      cout << " ** Ftol     = " << Ftol << endl;
      cout << "	** U 		= " << U << endl;
      cout << "	** dt 		= " << dt << endl;
      cout << "	** P 		= " << P << endl;
      cout << " ** phi (square boundaries)  = " << vertexPreferredPackingFraction2D() << endl;
      cout << " ** phi (polygonal boundaries) = " << vertexPreferredPackingFraction2D_polygon() << endl;
      cout << "	** alpha 	= " << alpha << endl;
      cout << "	** npPos 	= " << npPos << endl;
      cout << "	** npNeg 	= " << npNeg << endl;
      cout << "	** sxx  	= " << stress[0] << endl;
      cout << "	** syy 		= " << stress[1] << endl;
      cout << "	** sxy 		= " << stress[2] << endl;
    }

    // Adjust simulation based on net motion of degrees of freedom
    if (P > 0) {
      // increase positive counter
      npPos++;

      // reset negative counter
      npNeg = 0;

      // alter simulation if enough positive steps have been taken
      if (npPos > NDELAY) {
        // change time step
        if (dt * finc < dtmax)
          dt *= finc;

        // decrease alpha
        alpha *= falpha;
      }
    } else {
      // reset positive counter
      npPos = 0;

      // increase negative counter
      npNeg++;

      // check if simulation is stuck
      if (npNeg > NNEGMAX) {
        cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
        exit(1);
      }

      // take half step backwards, reset velocities
      for (i = 0; i < vertDOF; i++) {
        // take half step backwards
        x[i] -= 0.5 * dt * v[i];

        // reset vertex velocities
        v[i] = 0.0;
      }

      // decrease time step if past initial delay
      if (fireit > NDELAY) {
        // decrease time step
        if (dt * fdec > dtmin)
          dt *= fdec;

        // reset alpha
        alpha = alpha0;
      }
    }
    // VV VELOCITY UPDATE #1
    for (i = 0; i < vertDOF; i++)
      v[i] += 0.5 * dt * F[i];

    // compute fnorm, vnorm and P
    fnorm = 0.0;
    vnorm = 0.0;
    for (i = 0; i < vertDOF; i++) {
      fnorm += F[i] * F[i];
      vnorm += v[i] * v[i];
    }
    fnorm = sqrt(fnorm);
    vnorm = sqrt(vnorm);

    // update velocities (s.d. vs inertial dynamics) only if forces are acting
    if (fnorm > 0) {
      for (i = 0; i < vertDOF; i++)
        v[i] = (1 - alpha) * v[i] + alpha * (F[i] / fnorm) * vnorm;
    }
    // VV POSITION UPDATE
    for (i = 0; i < vertDOF; i++) {
      // update position
      x[i] += dt * v[i];

      // recenter in box
      if (x[i] > L[i % NDIM] && pbc[i % NDIM])
        x[i] -= L[i % NDIM];
      else if (x[i] < 0 && pbc[i % NDIM])
        x[i] += L[i % NDIM];
    }

    // update forces (function passed as argument)
    CALL_MEMBER_FN(*this, forceCall)
    ();

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++)
      v[i] += 0.5 * F[i] * dt;

    // update fcheck based on fnorm (= force per degree of freedom)
    fcheck = 0.0;
    for (i = 0; i < vertDOF; i++)
      fcheck += F[i] * F[i];
    fcheck = sqrt(fcheck / vertDOF);

    // update iterator
    fireit++;
  }
  // check if FIRE converged
  if (fireit == itmax) {
    cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
    cout << "fcheck = " << fcheck << '\n';
    if (fcheck < 1e-2) {
      cout << "fcheck < 1e-2 but fcheck > ftol, exceeded number of iterations, but small enough, proceeding \n";
    } else
      cout << "fcheck > 1e-2 and fcheck > ftol, exceeded number of iterations, deciding to proceed anyway, proceed with caution\n";
    // exit(1);
  } else {
    cout << endl;
    cout << "===========================================" << endl;
    cout << " 	F I R E 						" << endl;
    cout << "		M I N I M I Z A T I O N 	" << endl;
    cout << "	C O N V E R G E D! 				" << endl;
    cout << "===========================================" << endl;
    cout << endl;
    cout << "	** fireit 	= " << fireit << endl;
    cout << "	** fcheck 	= " << fcheck << endl;
    cout << "	** U 		= " << U << endl;

    cout << "	** fnorm	= " << fnorm << endl;
    cout << "	** vnorm 	= " << vnorm << endl;
    cout << "	** dt 		= " << dt << endl;
    cout << "	** P 		= " << P << endl;
    cout << "	** alpha 	= " << alpha << endl;
    cout << "	** sxx  	= " << stress[0] << endl;
    cout << "	** syy 		= " << stress[1] << endl;
    cout << "	** sxy 		= " << stress[2] << endl;
    cout << endl
         << endl;
  }
}

void dpm::vertexNVE2D(ofstream& enout, dpmMemFn forceCall, double T, double dt0, int NT, int NPRINTSKIP) {
  // local variables
  int t, i;
  double K, simclock;

  // set time step magnitude
  setdt(dt0);

  // initialize time keeper
  simclock = 0.0;

  // initialize velocities
  drawVelocities2D(T);

  // loop over time, print energy
  for (t = 0; t < NT; t++) {
    // VV VELOCITY UPDATE #1
    for (i = 0; i < vertDOF; i++)
      v[i] += 0.5 * dt * F[i];

    // VV POSITION UPDATE
    for (i = 0; i < vertDOF; i++) {
      // update position
      x[i] += dt * v[i];

      // recenter in box
      if (x[i] > L[i % NDIM] && pbc[i % NDIM])
        x[i] -= L[i % NDIM];
      else if (x[i] < 0 && pbc[i % NDIM])
        x[i] += L[i % NDIM];
    }

    // FORCE UPDATE
    CALL_MEMBER_FN(*this, forceCall)
    ();

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++)
      v[i] += 0.5 * F[i] * dt;

    // update sim clock
    simclock += dt;

    // print to console and file
    if (t % NPRINTSKIP == 0) {
      // compute kinetic energy
      K = vertexKineticEnergy();

      // print to console
      cout << endl
           << endl;
      cout << "===============================" << endl;
      cout << "	D P M  						" << endl;
      cout << " 			 					" << endl;
      cout << "		N V E 					" << endl;
      cout << "===============================" << endl;
      cout << endl;
      cout << "	** t / NT	= " << t << " / " << NT << endl;
      cout << "	** U 		= " << setprecision(12) << U << endl;
      cout << "	** K 		= " << setprecision(12) << K << endl;
      cout << "	** E 		= " << setprecision(12) << U + K << endl;

      // print to energy file
      cout << "** printing energy" << endl;
      enout << setw(w) << left << t;
      enout << setw(wnum) << left << simclock;
      enout << setw(wnum) << setprecision(12) << U;
      enout << setw(wnum) << setprecision(12) << K;
      enout << setw(wnum) << setprecision(12) << U + K;
      enout << endl;

      // print to configuration only if position file is open
      if (posout.is_open())
        printConfiguration2D();
    }
  }
}

/******************************

        D P M

                P R O T O C O L S

*******************************/

void dpm::vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0) {
  // local variables
  int it = 0, itmax = 1e3;
  double phi0 = vertexPreferredPackingFraction2D();
  double scaleFactor = 1.0, P, Sxy;

  // loop while phi0 < phi0Target
  while (phi0 < phi0Target && it < itmax) {
    // scale particle sizes
    scaleParticleSizes2D(scaleFactor);

    // update phi0
    phi0 = vertexPreferredPackingFraction2D();

    // relax configuration (pass member function force update)
    vertexFIRE2D(forceCall, Ftol, dt0);

    // get scale factor
    scaleFactor = sqrt((phi0 + dphi0) / phi0);

    // get updated pressure
    P = 0.5 * (stress[0] + stress[1]);
    Sxy = stress[2];

    // print to console
    if (it % 50 == 0) {
      cout << endl
           << endl;
      cout << "===============================" << endl;
      cout << "								" << endl;
      cout << " 	C O M P R E S S I O N 		" << endl;
      cout << "								" << endl;
      cout << "	P R O T O C O L 	  		" << endl;
      cout << "								" << endl;
      cout << "===============================" << endl;
      cout << endl;
      cout << "	** it 			= " << it << endl;
      cout << "	** phi0 curr	= " << phi0 << endl;
      if (phi0 + dphi0 < phi0Target)
        cout << "	** phi0 next 	= " << phi0 + dphi0 << endl;
      cout << "	** P 			= " << P << endl;
      cout << "	** Sxy 			= " << Sxy << endl;
      cout << "	** U 			= " << U << endl;
      // printConfiguration2D();
      cout << endl
           << endl;

      // update iterate
      it++;
    }
  }
}

void dpm::vertexCompress2Target2D_polygon(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0) {
  // same as vertexCompress2Target2D, but with polygonal boundaries (affects packing fraction calculation, and expects forceCall to
  //  account for polygonal boundary forces
  // local variables
  int it = 0, itmax = 1e4;
  double phi0 = vertexPreferredPackingFraction2D_polygon();
  double scaleFactor = 1.0, P, Sxy;

  // loop while phi0 < phi0Target
  while (phi0 < phi0Target && it < itmax) {
    // scale particle sizes
    scaleParticleSizes2D(scaleFactor);

    // update phi0
    phi0 = vertexPreferredPackingFraction2D_polygon();
    // relax configuration (pass member function force update)
    // make sure that forceCall is a force routine that includes a call to evaluatePolygonalWallForces
    vertexFIRE2D(forceCall, Ftol, dt0);

    // get scale factor
    scaleFactor = sqrt((phi0 + dphi0) / phi0);

    // get updated pressure
    P = 0.5 * (stress[0] + stress[1]);
    Sxy = stress[2];

    // print to console
    if (it % 50 == 0) {
      cout << " 	C O M P R E S S I O N 		" << endl;
      cout << "	** it 			= " << it << endl;
      cout << "	** phi0 curr	= " << phi0 << endl;
      if (phi0 + dphi0 < phi0Target)
        cout << "	** phi0 next 	= " << phi0 + dphi0 << endl;
      cout << "	** P 			= " << P << endl;
      cout << "	** Sxy 			= " << Sxy << endl;
      cout << "	** U 			= " << U << endl;
      // printConfiguration2D();
      cout << endl
           << endl;
    }

    // update iterate
    it++;
  }
}

void dpm::shrinkPolyWall(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0) {
  // shrink polygonal boundary walls to achieve desired packing fraction. preserves initial DP area, which is important because I scale my simulations using a0=1 and a0~25 micron^2
  //  account for polygonal boundary forces
  // local variables
  int it = 0, itmax = 1e4;
  double phi0 = vertexPreferredPackingFraction2D_polygon();
  double scaleFactor = 1.0, P, Sxy;

  string shrinkPositionsFile = "shrinkPositions.txt";
  ofstream shrinkStream(shrinkPositionsFile);

  // loop while phi0 < phi0Target
  while (phi0 < phi0Target && it < itmax) {
    for (int i = 0; i < NVTOT; i++) {
      shrinkStream << x[i * NDIM] << '\t' << x[i * NDIM + 1] << '\n';
    }
    shrinkStream << 0 << '\t' << 0 << '\n';

    if (it >= itmax)
      cout << "inside shrinkPolyWall, reached maxit. exiting compression steps\n";

    // scale poly boundary size to preserve a0=1
    scalePolyWallSize(1 / scaleFactor);
    cout << "scalePolyWallSize by factor of " << 1.0 / scaleFactor << "\n";
    cout << "vertex size = " << r[0] << '\n';
    cout << "cell size = " << a0[0] << '\n';

    // update phi0
    phi0 = vertexPreferredPackingFraction2D_polygon();
    // relax configuration (pass member function force update)
    // make sure that forceCall is a force routine that includes a call to evaluatePolygonalWallForces
    vertexFIRE2D(forceCall, Ftol, dt0);
    // vertexDampedMD(forceCall, dt0, 10.0, 0);

    // get scale factor
    scaleFactor = sqrt((phi0 + dphi0) / phi0);

    // get updated pressure
    P = 0.5 * (stress[0] + stress[1]);
    Sxy = stress[2];

    // print to console
    if (it % NSKIP == 0) {
      cout << " 	C O M P R E S S I O N 		" << endl;
      cout << "	** it 			= " << it << endl;
      cout << "	** phi0 curr	= " << phi0 << endl;
      if (phi0 + dphi0 < phi0Target)
        cout << "	** phi0 next 	= " << phi0 + dphi0 << endl;
      cout << "	** P 			= " << P << endl;
      cout << "	** Sxy 			= " << Sxy << endl;
      cout << "	** U 			= " << U << endl;
      cout << endl
           << endl;
    }
    // printConfiguration2D();

    // update iterate
    it++;
  }
}

void dpm::scalePolyWallSize(double scaleFactor) {
  double distanceBeforeScaling, changeInDistance;
  for (int tissueIt = 0; tissueIt < poly_bd_x.size(); tissueIt++) {
    double cx = 0, cy = 0, tempx = 0, tempy = 0;
    vector<double> poly_x = poly_bd_x[tissueIt];
    vector<double> poly_y = poly_bd_y[tissueIt];
    // calculate the center of the polygon
    for (int i = 0; i < poly_x.size(); i++) {
      cx += poly_x[i];
      cy += poly_y[i];
    }
    cx /= poly_x.size();
    cy /= poly_x.size();

    // record vertex where distance has max change, in order to scale the scaleFactor to less than a vertex size
    // this will prevent compression of boundaries from phasing through vertices
    double maxChangeInDistance = 0, maxDistanceBeforeScaling = 0;
    int vertexOfMaxChange = -1;
    for (int i = 0; i < poly_bd_x[tissueIt].size(); i++) {
      tempx = poly_bd_x[tissueIt][i] - cx;
      tempy = poly_bd_y[tissueIt][i] - cy;
      distanceBeforeScaling = sqrt(pow(tempx, 2) + pow(tempy, 2));
      changeInDistance = fabs(sqrt(pow(tempx * scaleFactor, 2) + pow(tempy * scaleFactor, 2)) - distanceBeforeScaling);
      if (changeInDistance > maxChangeInDistance) {
        vertexOfMaxChange = i;
        maxDistanceBeforeScaling = distanceBeforeScaling;
        maxChangeInDistance = changeInDistance;
      }
    }
    if (maxChangeInDistance > r[0] / 2) {
      // must have scaleFactor move boundary less than a vertex radius.
      scaleFactor = 1 - r[0] / (2 * maxDistanceBeforeScaling);
      cout << "modifying scaleFactor to be " << scaleFactor << '\n';
    }

    // new vertex is given by scaling the separation vector between center and vertex,
    //    and then adding it to the center
    for (int i = 0; i < poly_x.size(); i++) {
      poly_bd_x[tissueIt][i] = (poly_x[i] - cx) * scaleFactor + cx;
      poly_bd_y[tissueIt][i] = (poly_y[i] - cy) * scaleFactor + cy;
    }
  }
}

void dpm::vertexJamming2D(dpmMemFn forceCall, double Ftol, double Ptol, double dt0, double dphi0, bool plotCompression) {
  // local variables
  int k = 0, nr;
  bool jammed = 0, overcompressed = 0, undercompressed = 0;
  double pcheck, phi0, rH, r0, rL, rho0, scaleFactor = 1.0;
  // double pcheck, phi0, rH, r0, rL, rho0, scaleFactor;

  // initialize binary root search parameters
  r0 = sqrt(a0.at(0));
  rH = -1;
  rL = -1;

  // initialize preferred packing fraction
  phi0 = vertexPreferredPackingFraction2D();

  // save initial state
  vector<double> xsave(vertDOF, 0.0);
  vector<double> rsave(vertDOF, 0.0);
  vector<double> l0save(vertDOF, 0.0);
  vector<double> t0save(vertDOF, 0.0);
  vector<double> a0save(vertDOF, 0.0);

  xsave = x;
  rsave = r;
  l0save = l0;
  t0save = t0;
  a0save = a0;

  // loop until jamming is found
  while (!jammed && k < itmax) {
    // set length scale by 1st particle preferred area
    rho0 = sqrt(a0.at(0));

    // relax configuration (pass member function force update)
    vertexFIRE2D(forceCall, Ftol, dt0);

    // update pressure
    pcheck = 0.5 * (stress[0] + stress[1]);

    // remove rattlers
    nr = removeRattlers();

    // boolean checks for jamming
    undercompressed = ((pcheck < 2.0 * Ptol && rH < 0) || (pcheck < Ptol && rH > 0));
    overcompressed = (pcheck > 2.0 * Ptol);
    jammed = (pcheck < 2.0 * Ptol && pcheck > Ptol && rH > 0 && rL > 0);

    // output to console
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
    cout << "===============================================" << endl
         << endl;
    cout << " 	Q U A S I S T A T I C  						" << endl;
    cout << " 	  	I S O T R O P I C 						" << endl;
    cout << "			C O M P R E S S I O N 				" << endl
         << endl;
    cout << "===============================================" << endl;
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
    cout << endl;
    cout << "	* k 			= " << k << endl;
    cout << "	* phi0 			= " << phi0 << endl;
    cout << "	* phi 			= " << vertexPackingFraction2D() << endl;
    cout << "	* scaleFactor 	= " << scaleFactor << endl;
    cout << "	* r0 			= " << r0 << endl;
    cout << "	* rH 			= " << rH << endl;
    cout << "	* rL 			= " << rL << endl;
    cout << "	* pcheck 		= " << pcheck << endl;
    cout << "	* U 		 	= " << U << endl;
    cout << "	* Nvv  			= " << vvContacts() << endl;
    cout << "	* Ncc 			= " << ccContacts() << endl;
    cout << "	* # of rattlers = " << nr << endl
         << endl;
    cout << "	* undercompressed = " << undercompressed << endl;
    cout << "	* overcompressed = " << overcompressed << endl;
    cout << "	* jammed = " << jammed << endl
         << endl;
    if (plotCompression)
      printConfiguration2D();
    cout << endl
         << endl;

    // update particle scaleFactor based on target check
    if (rH < 0) {
      // if still undercompressed, then grow until overcompressed found
      if (undercompressed) {
        r0 = rho0;
        scaleFactor = sqrt((phi0 + dphi0) / phi0);
      }
      // if first overcompressed, decompress by dphi/2 until unjamming
      else if (overcompressed) {
        // current = upper bound length scale r
        rH = rho0;

        // save first overcompressed state
        r0 = rH;
        xsave = x;
        rsave = r;
        l0save = l0;
        t0save = t0;
        a0save = a0;

        // shrink particle sizes
        scaleFactor = sqrt((phi0 - 0.5 * dphi0) / phi0);

        // print to console
        cout << "	-- -- overcompressed for the first time, scaleFactor = " << scaleFactor << endl;
      }
    } else {
      if (rL < 0) {
        // if first undercompressed, save last overcompressed state, begin root search
        if (undercompressed) {
          // current = new lower bound length scale r
          rL = rho0;

          // load state
          x = xsave;
          r = rsave;
          l0 = l0save;
          t0 = t0save;
          a0 = a0save;

          // compute new scale factor by root search
          scaleFactor = 0.5 * (rH + rL) / r0;

          // print to console
          cout << "	-- -- undercompressed for the first time, scaleFactor = " << scaleFactor << endl;
          cout << "	-- -- BEGINNING ROOT SEARCH IN ENTHALPY MIN PROTOCOL..." << endl;
        }
        // if still overcompressed, decrement again
        else if (overcompressed) {
          // current = upper bound length scale r
          rH = rho0;

          // save overcompressed state
          r0 = rH;
          xsave = x;
          rsave = r;
          l0save = l0;
          t0save = t0;
          a0save = a0;

          // keep shrinking at same rate until unjamming
          scaleFactor = sqrt((phi0 - 0.5 * dphi0) / phi0);

          // print to console
          cout << "	-- -- overcompressed, still no unjamming, scaleFactor = " << scaleFactor << endl;
        }
      } else {
        // if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
        if (undercompressed) {
          // current = new lower bound length scale r
          rL = rho0;

          // load state
          x = xsave;
          r = rsave;
          l0 = l0save;
          t0 = t0save;
          a0 = a0save;

          // compute new scale factor by root search
          scaleFactor = 0.5 * (rH + rL) / r0;

          // print to console
          cout << "	-- -- undercompressed, scaleFactor = " << scaleFactor << endl;
        } else if (overcompressed) {
          // current = upper bound length scale r
          rH = rho0;

          // load state
          x = xsave;
          r = rsave;
          l0 = l0save;
          t0 = t0save;
          a0 = a0save;

          // compute new scale factor
          scaleFactor = 0.5 * (rH + rL) / r0;

          // print to console
          cout << "	-- -- overcompressed, scaleFactor = " << scaleFactor << endl;
        } else if (jammed) {
          cout << "	** At k = " << k << ", target pressure found!" << endl;
          cout << " WRITING ENTHALPY-MINIMIZED CONFIG TO FILE" << endl;
          cout << " ENDING COMPRESSION SIMULATION" << endl;
          scaleFactor = 1.0;
          if (!plotCompression)
            printConfiguration2D();
          break;
        }
      }
    }

    // scale particle sizes
    scaleParticleSizes2D(scaleFactor);

    // update packing fraction
    phi0 = vertexPreferredPackingFraction2D();

    // update iterate
    k++;
  }
}

void dpm::saveConfiguration(std::vector<double>& positionVector) {
  // save configuration data (for now, just vertex positions x)
  // if anything else changes like r, l0, NCELLS, etc. then this will be bugged.
  positionVector = x;
}

void dpm::loadConfiguration(std::vector<double>& positionVector) {
  // load position vector, with some checks to make sure crucial data hasn't been edited that will blow up the simulation
  // assumes we have a force balanced state with v = F = 0
  assert(positionVector.size() == x.size());
  x = positionVector;
  fill(v.begin(), v.end(), 0.0);
  fill(F.begin(), F.end(), 0.0);
}

/******************************

        D P M

                H E S S I A N

*******************************/

// wrapper function for total hessian
// note: dynamical matrix M = H - S
void dpm::dpmHessian2D(Eigen::MatrixXd& H, Eigen::MatrixXd& S) {
  // local variables
  int k, l;

  // print something to the console
  cout << "** Computing Hessian for configuration in dpmHessian2D ..." << endl;

  // initialize all possible matrices
  Eigen::MatrixXd Ha(vertDOF, vertDOF);   // stiffness matrix for area term
  Eigen::MatrixXd Sa(vertDOF, vertDOF);   // stress matrix for area term
  Eigen::MatrixXd Hl(vertDOF, vertDOF);   // stiffness matrix for perimeter term
  Eigen::MatrixXd Sl(vertDOF, vertDOF);   // stress matrix for perimeter term
  Eigen::MatrixXd Hb(vertDOF, vertDOF);   // stiffness matrix for bending energy
  Eigen::MatrixXd Sb(vertDOF, vertDOF);   // stress matrix for bending term
  Eigen::MatrixXd Hvv(vertDOF, vertDOF);  // stiffness matrix for interaction terms
  Eigen::MatrixXd Svv(vertDOF, vertDOF);  // stress matrix for interaction terms

  // initialize all matrices to be 0 initially
  for (k = 0; k < vertDOF; k++) {
    for (l = 0; l < vertDOF; l++) {
      Ha(k, l) = 0.0;
      Sa(k, l) = 0.0;
      Hl(k, l) = 0.0;
      Sl(k, l) = 0.0;
      Hb(k, l) = 0.0;
      Sb(k, l) = 0.0;
      Hvv(k, l) = 0.0;
      Svv(k, l) = 0.0;
      S(k, l) = 0.0;
      H(k, l) = 0.0;
    }
  }

  // find matrix elements for each term
  if (ka > 0)
    dpmAreaHessian2D(Ha, Sa);

  if (kl > 0)
    dpmPerimeterHessian2D(Hl, Sl);

  // if (kb > 0)
  // dpmBendingHessian2D(Hb,Sb);

  if (kc > 0)
    dpmRepulsiveHarmonicSprings2D(Hvv, Svv);

  // construct matrices
  for (k = 0; k < vertDOF; k++) {
    for (l = 0; l < vertDOF; l++) {
      H(k, l) = Ha(k, l) + Hl(k, l) + Hb(k, l) + Hvv(k, l);
      S(k, l) = -Sa(k, l) - Sl(k, l) - Sb(k, l) - Svv(k, l);
    }
  }
}

// construct hessian for area term
void dpm::dpmAreaHessian2D(Eigen::MatrixXd& Ha, Eigen::MatrixXd& Sa) {
  // local variables
  int nvtmp, ci, vim1, vi, vip1, vjm1, vj, vjp1;
  int kxm1, kx, kxp1, kym1, ky, kyp1, lxm1, lym1, lx, ly, lxp1, lyp1;
  double rho0, a0tmp, a02tmp, da, da_dxi, da_dyi, da_dxj, da_dyj;
  double lim1x, lix, liy, lim1y, ljm1x, ljm1y, ljx, ljy;

  // loop over cells
  rho0 = sqrt(a0[0]);
  for (ci = 0; ci < NCELLS; ci++) {
    // shape parameters for ci
    nvtmp = nv[ci];
    a0tmp = a0[ci];
    a02tmp = a0tmp * a0tmp;

    // fractional area strain
    da = (area(ci) / a0tmp) - 1.0;

    // loop over vertices
    for (vi = 0; vi < nvtmp; vi++) {
      // wrap vertices
      vim1 = (vi - 1 + nvtmp) % nvtmp;
      vip1 = (vi + 1) % nvtmp;

      // matrix indices
      kxm1 = NDIM * (gindex(ci, vim1));
      kym1 = NDIM * (gindex(ci, vim1)) + 1;

      kx = NDIM * (gindex(ci, vi));
      ky = NDIM * (gindex(ci, vi)) + 1;

      kxp1 = NDIM * (gindex(ci, vip1));
      kyp1 = NDIM * (gindex(ci, vip1)) + 1;

      // segment elements
      lim1x = x[kx] - x[kxm1];
      lim1y = x[ky] - x[kym1];

      lix = x[kxp1] - x[kx];
      liy = x[kyp1] - x[ky];

      if (pbc[0]) {
        lim1x -= L[0] * round(lim1x / L[0]);
        lix -= L[0] * round(lix / L[0]);
      }
      if (pbc[1]) {
        lim1y -= L[1] * round(lim1y / L[1]);
        liy -= L[1] * round(liy / L[1]);
      }

      // stress matrix
      Sa(kx, kyp1) = 0.5 * da * ((rho0 * rho0) / a0tmp);
      Sa(ky, kxp1) = -0.5 * da * ((rho0 * rho0) / a0tmp);

      Sa(kyp1, kx) = Sa(kx, kyp1);
      Sa(kxp1, ky) = Sa(ky, kxp1);

      // area derivatives (for stiffness matrix)
      da_dxi = 0.5 * (liy + lim1y);
      da_dyi = -0.5 * (lim1x + lix);

      // loop over other vertices, for area elasticity stiffness matrix
      for (vj = vi; vj < nvtmp; vj++) {
        // wrap jp1 and jm1
        vjp1 = (vj + 1) % nvtmp;
        vjm1 = (vj - 1 + nvtmp) % nvtmp;

        // dof elements
        lxm1 = NDIM * (gindex(ci, vjm1));
        lym1 = lxm1 + 1;

        lx = NDIM * (gindex(ci, vj));
        ly = lx + 1;

        lxp1 = NDIM * (gindex(ci, vjp1));
        lyp1 = lxp1 + 1;

        // j segments
        ljm1x = x[lx] - x[lxm1];
        if (pbc[0])
          ljm1x -= L[0] * round(ljm1x / L[0]);

        ljm1y = x[ly] - x[lym1];
        if (pbc[1])
          ljm1y -= L[1] * round(ljm1y / L[1]);

        ljx = x[lxp1] - x[lx];
        if (pbc[0])
          ljx -= L[0] * round(ljx / L[0]);

        ljy = x[lyp1] - x[ly];
        if (pbc[1])
          ljy -= L[1] * round(ljy / L[1]);

        // area derivatives
        da_dxj = 0.5 * (ljy + ljm1y);
        da_dyj = -0.5 * (ljm1x + ljx);

        // stiffness matrix
        Ha(kx, lx) = da_dxi * da_dxj * ((rho0 * rho0) / a02tmp);
        Ha(kx, ly) = da_dxi * da_dyj * ((rho0 * rho0) / a02tmp);

        Ha(ky, lx) = da_dyi * da_dxj * ((rho0 * rho0) / a02tmp);
        Ha(ky, ly) = da_dyi * da_dyj * ((rho0 * rho0) / a02tmp);

        Ha(lx, kx) = Ha(kx, lx);
        Ha(ly, kx) = Ha(kx, ly);

        Ha(lx, ky) = Ha(ky, lx);
        Ha(ly, ky) = Ha(ky, ly);
      }
    }
  }
}

// construct hessian for perimeter term
void dpm::dpmPerimeterHessian2D(Eigen::MatrixXd& Hl, Eigen::MatrixXd& Sl) {
  // local variables
  int nvtmp, ci, vim1, vi, vip1;
  int kxm1, kx, kxp1, kym1, ky, kyp1;
  double lim1x, lim1y, lix, liy, lim1, li, dlim1, dli, ulim1x, ulim1y, ulix, uliy;
  double l0im1, l0im1_sq, l0i, l0i_sq;
  double rho0, Kl;

  // loop over cells
  rho0 = sqrt(a0[0]);
  for (ci = 0; ci < NCELLS; ci++) {
    // number of vertices
    nvtmp = nv[ci];

    // prefactor scaled by length, will come out as dimensionless
    Kl = kl * (rho0 * rho0);

    for (vi = 0; vi < nvtmp; vi++) {
      // wrap vertices
      vim1 = (vi - 1 + nvtmp) % nvtmp;
      vip1 = (vi + 1) % nvtmp;

      // matrix indices
      kxm1 = NDIM * (gindex(ci, vim1));
      kym1 = NDIM * (gindex(ci, vim1)) + 1;

      kx = NDIM * (gindex(ci, vi));
      ky = NDIM * (gindex(ci, vi)) + 1;

      kxp1 = NDIM * (gindex(ci, vip1));
      kyp1 = NDIM * (gindex(ci, vip1)) + 1;

      // segment elements
      lim1x = x[kx] - x[kxm1];
      lim1y = x[ky] - x[kym1];

      lix = x[kxp1] - x[kx];
      liy = x[kyp1] - x[ky];

      if (pbc[0]) {
        lim1x -= L[0] * round(lim1x / L[0]);
        lix -= L[0] * round(lix / L[0]);
      }
      if (pbc[1]) {
        lim1y -= L[1] * round(lim1y / L[1]);
        liy -= L[1] * round(liy / L[1]);
      }

      // segment lengths
      lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);
      li = sqrt(lix * lix + liy * liy);

      // segment strains
      l0im1 = l0[gindex(ci, vim1)];
      l0i = l0[gindex(ci, vi)];

      dlim1 = (lim1 / l0im1) - 1.0;
      dli = (li / l0i) - 1.0;

      l0im1_sq = l0im1 * l0im1;
      l0i_sq = l0i * l0i;

      // -- PERIMETER SPRINGS

      // unit vectors
      ulim1x = lim1x / lim1;
      ulim1y = lim1y / lim1;

      ulix = lix / li;
      uliy = liy / li;

      // 	STIFFNESS MATRIX

      // main diagonal
      Hl(kx, kx) = Kl * ((ulix * ulix) / l0i_sq + (ulim1x * ulim1x) / l0im1_sq);
      Hl(ky, ky) = Kl * ((uliy * uliy) / l0i_sq + (ulim1y * ulim1y) / l0im1_sq);

      Hl(kx, ky) = Kl * ((ulix * uliy) / l0i_sq + (ulim1x * ulim1y) / l0im1_sq);
      Hl(ky, kx) = Hl(kx, ky);

      // 1off diagonal
      Hl(kx, kxp1) = -Kl * (ulix * ulix) / l0i_sq;
      Hl(ky, kyp1) = -Kl * (uliy * uliy) / l0i_sq;

      Hl(kx, kyp1) = -Kl * (ulix * uliy) / l0i_sq;
      Hl(ky, kxp1) = Hl(kx, kyp1);

      // enforce symmetry in lower triangle
      Hl(kxp1, kx) = Hl(kx, kxp1);
      Hl(kyp1, ky) = Hl(ky, kyp1);

      Hl(kyp1, kx) = Hl(kx, kyp1);
      Hl(kxp1, ky) = Hl(ky, kxp1);

      // 	STRESS MATRIX

      // main diagonal
      Sl(kx, kx) = Kl * (dlim1 * ((ulim1y * ulim1y) / (l0im1 * lim1)) + dli * ((uliy * uliy) / (l0i * li)));
      Sl(ky, ky) = Kl * (dlim1 * ((ulim1x * ulim1x) / (l0im1 * lim1)) + dli * ((ulix * ulix) / (l0i * li)));

      Sl(kx, ky) = -Kl * (dlim1 * ((ulim1x * ulim1y) / (l0im1 * lim1)) + dli * ((ulix * uliy) / (l0i * li)));
      Sl(ky, kx) = Sl(kx, ky);

      // 1off diagonal
      Sl(kx, kxp1) = -Kl * dli * ((uliy * uliy) / (l0i * li));
      Sl(ky, kyp1) = -Kl * dli * ((ulix * ulix) / (l0i * li));

      Sl(kx, kyp1) = Kl * dli * ((ulix * uliy) / (l0i * li));
      Sl(ky, kxp1) = Sl(kx, kyp1);

      // enforce symmetry in lower triangle
      Sl(kxp1, kx) = Sl(kx, kxp1);
      Sl(kyp1, ky) = Sl(ky, kyp1);

      Sl(kyp1, kx) = Sl(kx, kyp1);
      Sl(kxp1, ky) = Sl(ky, kxp1);
    }
  }
}

// TO-DO: need to make hessian function for bending term (th - th0)^2
// void dpm::dpmBendingHessian2D(Eigen::MatrixXd& Hb, Eigen::MatrixXd& Sb){
// 	// local variables

// 	//
// }

// construct hessian for interaction term
void dpm::dpmRepulsiveHarmonicSprings2D(Eigen::MatrixXd& Hvv, Eigen::MatrixXd& Svv) {
  // local variables
  int ci, cj, vi, vj, gi, gj;
  int mxi, myi, mxj, myj;
  double rho0, sij, dx, dy, rij, kij, h, uxij, uyij;

  // loop over cell pairs
  rho0 = sqrt(a0[0]);
  for (ci = 0; ci < NCELLS; ci++) {
    for (cj = ci + 1; cj < NCELLS; cj++) {
      // check if pair of cells is contact, only proceed if true
      if (cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2] > 0) {
        // loop over pairs of vertices on both cells, check for overlap, compute matrix elements
        for (vi = 0; vi < nv[ci]; vi++) {
          // matrix element indices (cell ci, vertex vi)
          gi = gindex(ci, vi);
          mxi = NDIM * gi;
          myi = mxi + 1;

          for (vj = 0; vj < nv[cj]; vj++) {
            // matrix element indices (cell cj, vertex vj)
            gj = gindex(cj, vj);
            mxj = NDIM * gj;
            myj = mxj + 1;

            // contact distance
            sij = r[gi] + r[gj];

            // get distance between vertices
            // particle distance
            dx = x[mxj] - x[mxi];
            if (pbc[0])
              dx -= L[0] * round(dx / L[0]);
            if (dx < sij) {
              dy = x[myj] - x[myi];
              if (pbc[1])
                dy -= L[1] * round(dy / L[1]);
              if (dy < sij) {
                rij = sqrt(dx * dx + dy * dy);
                if (rij < sij) {
                  // spring constant
                  kij = (kc * rho0 * rho0) / (sij * rij);

                  // dimensionless overlap
                  h = rij / sij;

                  // derivatives of distance w.r.t. coordinates
                  uxij = dx / rij;
                  uyij = dy / rij;

                  // compute stiffness and stress matrices (off diagonal, enforce symmetry in lower triangles)

                  // -- stiffness matrix
                  Hvv(mxi, mxj) = -((kc * rho0 * rho0) / (sij * sij)) * (uxij * uxij);
                  Hvv(myi, myj) = -((kc * rho0 * rho0) / (sij * sij)) * (uyij * uyij);
                  Hvv(mxi, myj) = -((kc * rho0 * rho0) / (sij * sij)) * (uxij * uyij);
                  Hvv(myi, mxj) = -((kc * rho0 * rho0) / (sij * sij)) * (uyij * uxij);

                  Hvv(mxj, mxi) = Hvv(mxi, mxj);
                  Hvv(myj, myi) = Hvv(myi, myj);
                  Hvv(mxj, myi) = Hvv(myi, mxj);
                  Hvv(myj, mxi) = Hvv(mxi, myj);

                  // -- stress matrix
                  Svv(mxi, mxj) = kij * (1.0 - h) * (uyij * uyij);
                  Svv(myi, myj) = kij * (1.0 - h) * (uxij * uxij);
                  Svv(mxi, myj) = -kij * (1.0 - h) * (uxij * uyij);
                  Svv(myi, mxj) = -kij * (1.0 - h) * (uxij * uyij);

                  Svv(mxj, mxi) = Svv(mxi, mxj);
                  Svv(myj, myi) = Svv(myi, myj);
                  Svv(mxj, myi) = Svv(myi, mxj);
                  Svv(myj, mxi) = Svv(mxi, myj);

                  // add to diagonal, using off diagonals and reciprocity

                  // -- stiffness matrix
                  Hvv(mxi, mxi) -= Hvv(mxi, mxj);
                  Hvv(myi, myi) -= Hvv(myi, myj);
                  Hvv(mxi, myi) -= Hvv(mxi, myj);
                  Hvv(myi, mxi) -= Hvv(myi, mxj);

                  Hvv(mxj, mxj) -= Hvv(mxi, mxj);
                  Hvv(myj, myj) -= Hvv(myi, myj);
                  Hvv(mxj, myj) -= Hvv(mxi, myj);
                  Hvv(myj, mxj) -= Hvv(myi, mxj);

                  // -- stress matrix
                  Svv(mxi, mxi) -= Svv(mxi, mxj);
                  Svv(myi, myi) -= Svv(myi, myj);
                  Svv(mxi, myi) -= Svv(mxi, myj);
                  Svv(myi, mxi) -= Svv(myi, mxj);

                  Svv(mxj, mxj) -= Svv(mxi, mxj);
                  Svv(myj, myj) -= Svv(myi, myj);
                  Svv(mxj, myj) -= Svv(mxi, myj);
                  Svv(myj, mxj) -= Svv(myi, mxj);
                }
              }
            }
          }
        }
      }
    }
  }
}

/******************************

        P R I N T   T O

        C O N S O L E  &  F I L E

*******************************/

void dpm::printContactMatrix() {
  int ci, cj;

  for (ci = 0; ci < NCELLS; ci++) {
    for (cj = 0; cj < NCELLS; cj++) {
      if (ci > cj)
        cout << setw(5) << cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];
      else if (ci < cj)
        cout << setw(5) << cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
      else
        cout << setw(5) << 0;
    }
    cout << endl;
  }
}

void dpm::printConfiguration2D() {
  // local variables
  int ci, cj, vi, gi, ctmp, zc, zv;
  double xi, yi, dx, dy, Lx, Ly;

  // check if pos object is open
  if (!posout.is_open()) {
    cerr << "** ERROR: in printConfiguration2D, posout is not open, but function call will try to use. Ending here." << endl;
    exit(1);
  } else
    cout << "** In printConfiguration2D, printing particle positions to file..." << endl;

  // save box sizes
  Lx = L.at(0);
  Ly = L.at(1);

  // print information starting information
  posout << setw(w) << left << "NEWFR"
         << " " << endl;
  posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << endl;
  posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction2D() << endl;

  // print box sizes
  posout << setw(w) << left << "BOXSZ";
  posout << setw(wnum) << setprecision(pnum) << left << Lx;
  posout << setw(wnum) << setprecision(pnum) << left << Ly;
  posout << endl;

  // print stress info
  posout << setw(w) << left << "STRSS";
  posout << setw(wnum) << setprecision(pnum) << left << stress.at(0);
  posout << setw(wnum) << setprecision(pnum) << left << stress.at(1);
  posout << setw(wnum) << setprecision(pnum) << left << stress.at(2);
  posout << endl;

  // print coordinate for rest of the cells
  for (ci = 0; ci < NCELLS; ci++) {
    // get cell contact data
    zc = 0;
    zv = 0;
    for (cj = 0; cj < NCELLS; cj++) {
      if (ci != cj) {
        // contact info from entry ci, cj
        if (ci < cj)
          ctmp = cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
        else
          ctmp = cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];

        // add to contact information
        zv += ctmp;
        if (ctmp > 0)
          zc++;
      }
    }

    // cell information
    posout << setw(w) << left << "CINFO";
    posout << setw(w) << left << nv.at(ci);
    posout << setw(w) << left << zc;
    posout << setw(w) << left << zv;
    posout << setw(wnum) << left << a0.at(ci);
    posout << setw(wnum) << left << area(ci);
    posout << setw(wnum) << left << perimeter(ci);
    posout << endl;

    // get initial vertex positions
    gi = gindex(ci, 0);
    xi = x.at(NDIM * gi);
    yi = x.at(NDIM * gi + 1);

    // place back in box center
    if (pbc[0])
      xi = fmod(xi, Lx);
    if (pbc[1])
      yi = fmod(yi, Ly);

    posout << setw(w) << left << "VINFO";
    posout << setw(w) << left << ci;
    posout << setw(w) << left << 0;

    // output initial vertex information
    posout << setw(wnum) << setprecision(pnum) << right << xi;
    posout << setw(wnum) << setprecision(pnum) << right << yi;
    posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
    posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
    posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
    posout << endl;

    // vertex information for next vertices
    for (vi = 1; vi < nv.at(ci); vi++) {
      // get global vertex index for next vertex
      gi++;

      // get next vertex positions
      dx = x.at(NDIM * gi) - xi;
      if (pbc[0])
        dx -= Lx * round(dx / Lx);
      xi += dx;

      dy = x.at(NDIM * gi + 1) - yi;
      if (pbc[1])
        dy -= Ly * round(dy / Ly);
      yi += dy;

      // Print indexing information
      posout << setw(w) << left << "VINFO";
      posout << setw(w) << left << ci;
      posout << setw(w) << left << vi;

      // output vertex information
      posout << setw(wnum) << setprecision(pnum) << right << xi;
      posout << setw(wnum) << setprecision(pnum) << right << yi;
      posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
      posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
      posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
      posout << endl;
    }
  }

  // print end frame
  posout << setw(w) << left << "ENDFR"
         << " " << endl;
}

void dpm::printHessianEigenvalues2D(ofstream& hessout, Eigen::MatrixXd& M) {
  // check if pos object is open
  if (!hessout.is_open()) {
    cerr << "** ERROR: in printMatrixEigenvalues2D, hessout is not open, but function call will try to use. Ending here." << endl;
    exit(1);
  } else
    cout << "** In printMatrixEigenvalues2D, printing particle positions to file..." << endl;

  // compute eigenvalues from matrix, plot
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> dynamicalMatrixEigenmodes(M);

  // print to file
  hessout << vertDOF << endl;
  hessout << dynamicalMatrixEigenmodes.eigenvalues() << endl;
}
