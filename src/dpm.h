#ifndef DPM_H
#define DPM_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

// pointer-to-member function call macro
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

class dpm;
typedef void (dpm::*dpmMemFn)(void);

// global constants
const double PI = 4.0 * atan(1.0);
const int nvmin = 12;

// printing constants
const int w = 10;
const int wnum = 25;
const int pnum = 14;

// FIRE constants
const double alpha0 = 0.2;
const double finc = 1.1;
const double fdec = 0.5;
const double falpha = 0.99;

const int NSKIP = 20000;
const int NMIN = 10;
const int NNEGMAX = 1000;
const int NDELAY = 20;
// int itmax = 1e7;

class dpm {
 protected:
  int itmax = 1e5;

  // int scalars
  int NCELLS;
  int NDIM;
  int NNN;
  int NVTOT;
  int vertDOF;

  // time step size
  double dt;

  // potential energy
  double U;
  std::vector<double> cellU;

  // particle spring constants
  double ka;
  double kl;
  double kb;
  double kc;

  // rheological parameters
  double maxwellRelaxationTime, taus;
  std::vector<double> vl0, Fl0, l00;

  // particle attraction constants
  double l1, l2;

  // boundary parameters
  std::vector<double> L;
  std::vector<bool> pbc;

  // alternative boundary parameters: non-rectangular boundaries
  std::vector<std::vector<double>> poly_bd_x;  // x coordinates of polygonal boundary condition (could be a triangle, square, n-0gon, star, etc.)
  std::vector<std::vector<double>> poly_bd_y;  // set poly_x and poly_y by writing a function like generateCircularBoundary

  // particle shape parameters
  std::vector<double> a0;
  std::vector<double> l0;
  std::vector<double> t0;
  std::vector<double> r;

  // indexing variables
  std::vector<int> nv;
  std::vector<int> szList;
  std::vector<int> im1;
  std::vector<int> ip1;

  // dynamical variables
  std::vector<double> x;
  std::vector<double> v;
  std::vector<double> F;
  double B;

  // macroscopic stress vector
  std::vector<double> stress;

  // local stress vector
  std::vector<std::vector<double>> fieldStress;
  std::vector<std::vector<double>> fieldStressCells;

  std::vector<std::vector<double>> fieldShapeStress;
  std::vector<std::vector<double>> fieldShapeStressCells;

  // contact network (vector, size N(N-1)/2), stores # vertex contacts between i-j (i,j are cells)
  //  aka flattened triangular matrix
  // cij is structured as follows: (0-1, 0-2, 0-3, ... ,0- (N-1), 1-2, 1-3, ..., 1- (N-1), 2-3,...)
  std::vector<int> cij;
  std::vector<std::vector<int>> numVertexContacts;

  // Box linked-list variables
  int NBX;
  std::vector<int> sb;
  std::vector<double> lb;
  std::vector<std::vector<int>> nn;
  std::vector<int> head;
  std::vector<int> last;
  std::vector<int> list;

  // output objects
  std::ofstream posout;

  void openFile(std::ofstream& file, const std::string& filename) {
    file.open(filename.c_str());
    if (!file.is_open()) {
      std::cerr << "ERROR: Could not open file " << filename << std::endl;
      std::cerr << "ERROR: Could not open file " << filename.c_str() << std::endl;
      exit(1);
    }
    std::cout << "** Opening file " << filename << " ..." << std::endl;
  }

 public:
  // Constructors and Destructors
  dpm(int n, int ndim, int seed);
  dpm(int n, int seed)
      : dpm(n, 2, seed) {}
  //~dpm();

  // -- G E T T E R S

  // main ints
  int getNCELLS() { return NCELLS; };
  int getNDIM() { return NDIM; };
  int getNNN() { return NNN; };
  int getNVTOT() { return NVTOT; };
  int getvertDOF() { return vertDOF; };
  int getNV(int ci) { return nv.at(ci); };

  // force parameters
  double getdt() { return dt; };
  double getka() { return ka; };
  double getkl() { return kl; };
  double getkb() { return kb; };
  double getkc() { return kc; };

  // static cell info
  double geta0(int ci) { return a0[ci]; };
  double getl0(int gi) { return l0[gi]; };
  double gett0(int gi) { return t0[gi]; };
  double getr(int gi) { return r[gi]; };

  // dynamic cell info
  double getx(int gi, int d) { return x[NDIM * gi + d]; };
  double getv(int gi, int d) { return v[NDIM * gi + d]; };
  double getF(int gi, int d) { return F[NDIM * gi + d]; };
  double getU() { return U; };

  // boundary variables
  double getL(int d) { return L.at(d); };
  bool getpbc(int d) { return pbc.at(d); };

  // cell shape indexing + information
  int gindex(int ci, int vi);
  void cindices(int& ci, int& vi, int gi);
  double area(int ci);
  double perimeter(int ci);
  void com2D(int ci, double& cx, double& cy);
  double vertexPackingFraction2D();
  double vertexPreferredPackingFraction2D();
  double vertexPreferredPackingFraction2D_polygon();
  double vertexKineticEnergy();
  int vvContacts();
  int ccContacts();
  void initializeFieldStress();

  // Setters
  void setitmax(int val) { itmax = val; };
  void setpbc(int d, bool val) { pbc.at(d) = val; };
  void setNCELLS(int val) { NCELLS = val; };
  void setdt(double val);
  void setka(double val) { ka = val; };
  void setkl(double val) { kl = val; };
  void setl00() { l00 = l0; };
  void setkb(double val) { kb = val; };
  void setkc(double val) { kc = val; };
  void setl1(double val) { l1 = val; };
  void setl2(double val) { l2 = val; };
  void setB(double val) { B = val; };
  void setMaxwellRelaxationTime(double val) { maxwellRelaxationTime = val; };
  void setTaus(double val) { taus = val; };
  void scaleL(int d, double val) { L.at(d) *= val; };
  void scaleRadius(double scalefactor) {
    for (int i = 0; i < NVTOT; i++) {
      r[i] *= scalefactor;
    }
  }
  void scaleVelocities(double scalefactor) {
    for (int i = 0; i < vertDOF; i++)
      v[i] *= scalefactor;
  }
  void displaceCell(int ci, double displaceX, double displaceY) {
    int firstIndex = szList[ci];
    for (int gi = firstIndex; gi < firstIndex + nv[ci]; gi++) {
      x[NDIM * gi] += displaceX;
      x[NDIM * gi + 1] += displaceY;
    }
  }
  void displaceVertex(int ci, int vi, double displaceX, double displaceY) {
    int firstIndex = szList[ci];
    int gi = firstIndex + vi;
    std::cout << "before displacement, vi gi = " << vi << '\t' << gi << '\t' << ", x y = " << x[NDIM * gi] << '\t' << x[NDIM * gi + 1] << '\n';
    x[NDIM * gi] += displaceX;
    x[NDIM * gi + 1] += displaceY;
    std::cout << "after displacement, vi gi = " << vi << '\t' << gi << '\t' << ", x y = " << x[NDIM * gi] << '\t' << x[NDIM * gi + 1] << '\n';
  }
  void setCellVelocity(int ci, double velocityX, double velocityY) {
    int firstIndex = szList[ci];
    for (int gi = firstIndex; gi < firstIndex + nv[ci]; gi++) {
      v[NDIM * gi] = velocityX;
      v[NDIM * gi + 1] = velocityY;
    }
  }

  void moveSimulationToPositiveCoordinates(double xshift = 0, double yshift = 0);

  // File openers
  void openPosObject(std::string& filename) { openFile(posout, filename); }

  // Initialize particles (two dimensions)
  void monodisperse2D(double calA0, int n);
  void bidisperse2D(double calA0, int nsmall, double smallfrac, double sizefrac);
  void gaussian2D(double dispersion, double calA0, int n1);
  void sinusoidalPreferredAngle(double thA, double thK);
  void initializeVertexShapeParameters(double calA0, int nref);
  void initializeVertexShapeParameters(std::vector<double> calA0, int nref);
  void initializeVertexIndexing2D();
  void initializePositions2D(double phi0, double Ftol, bool isFixedBoundary = false, double aspectRatio = 1.0, bool setUpCircularBoundary = false);
  void initializeAllPositions(std::string vertexPositionFile, int nref);
  void initializeFromConfigurationFile(std::string vertexPositionFile, double phi0);
  void initializeNeighborLinkedList2D(double boxLengthScale);
  void resizeNeighborLinkedList2D();

  // editing & updating
  void sortNeighborLinkedList2D();
  void scaleParticleSizes2D(double scaleFactor);
  int removeRattlers();
  void drawVelocities2D(double T);
  double distanceLineAndPoint(double x1, double y1, double x2, double y2, double x0, double y0);
  double distanceLinePointComponents(double x1, double y1, double x2, double y2, double x0, double y0, double& xcomp, double& ycomp);
  double linePointDistancesAndProjection(double x1, double y1, double x2, double y2, double x0, double y0, double& xcomp, double& ycomp, double& contactType);
  void generateCircularBoundary(int numEdges, double radius, double cx, double cy, std::vector<double>& poly_x, std::vector<double>& poly_y);
  void generateCircle(int numEdges, double cx, double cy, double r, std::vector<double>& poly_x, std::vector<double>& poly_y);
  void generateRectangularBoundary(double radius, double cx, double cy, std::vector<double>& poly_x, std::vector<double>& poly_y);
  void generateHorseshoeBoundary(double cx, double cy, std::vector<double>& poly_x, std::vector<double>& poly_y);
  void replaceCircularBoundary(int polyShapeID, double aspectRatio);
  std::vector<double> resample_polygon(std::vector<double> px, std::vector<double> py, double perimeter, int numPoints);
  bool isInsidePolygon(double x, double y, const std::vector<double>& poly_x, const std::vector<double>& poly_y);

  void scalePolyWallSize(double scaleFactor);
  void shrinkPolyWall(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);

  // force definitions
  void resetForcesAndEnergy();
  void shapeForces2D();
  void maxwellRelaxationRestLengths(std::vector<double>& l);
  void vertexRepulsiveForces2D();
  void vertexAttractiveForces2D();
  void evaluatePolygonalWallForces(const std::vector<double>& poly_x, const std::vector<double>& poly_y, bool attractionOn = false);

  // force updates
  void repulsiveForceUpdate();
  void attractiveForceUpdate();

  // simple integrators
  void vertexFIRE2D(dpmMemFn forceCall, double Ftol, double dt0);
  void vertexNVE2D(std::ofstream& enout, dpmMemFn forceCall, double T, double dt0, int NT, int NPRINTSKIP);

  // protocols
  void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void vertexCompress2Target2D_polygon(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void vertexJamming2D(dpmMemFn forceCall, double Ftol, double Ptol, double dt0, double dphi0, bool plotCompression);
  void saveConfiguration(std::vector<double>& positionVector);
  void loadConfiguration(std::vector<double>& positionVector);

  // hessian methods
  // note: dynamical matrix contribution is always M = H - S
  void dpmHessian2D(Eigen::MatrixXd& H, Eigen::MatrixXd& S);
  void dpmAreaHessian2D(Eigen::MatrixXd& Ha, Eigen::MatrixXd& Sa);
  void dpmPerimeterHessian2D(Eigen::MatrixXd& Hl, Eigen::MatrixXd& Sl);
  void dpmRepulsiveHarmonicSprings2D(Eigen::MatrixXd& Hvv, Eigen::MatrixXd& Svv);

  // print vertex information to file
  void printContactMatrix();
  void printConfiguration2D();
  void printHessianEigenvalues2D(std::ofstream& hessout, Eigen::MatrixXd& M);
};

#endif
