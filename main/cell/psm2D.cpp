// simulate transverse section, focusing on tissue structure
// 2D justification: cells don't move much in the Anterior-Posterior plane of the PSM
//   the measurements are taken in the 2D transverse section, so a 2D simulation would help
//   our understanding of tissue structure in those same 2D sections.
// Scientific question: Why do cells gain extracellular space in cadherin AND/OR fibronectin/fibrillin mutants?
//  cadherin mutants lose intercellular adhesion. fN/fbn mutants lose cell-ECM adhesion.
// My understanding - cell-ECM adhesion acts as a compression. The ECM proteins link up with the actin skeleton
//   of cells on the boundary of the PSM. This causes action like a purse-string at the boundary of the tissue
//   which might explain why PSM is rounded up like a cylinder.
/*
Compilation command:
g++ -O3 --std=c++11 -g -I src main/cell/psm2D.cpp src/dpm.cpp src/cell.cpp -o main/cell/psm2D.o
run command:

att_arr=(0.05)
att2_arr=(0.05)
t_stress_arr=(10000.0)
v0=0.1
phi_arr=(0.8)
tau_abp=1.0
gamma_arr=(0)
kon_arr=(1.0)
koff_arr=(0.5)
kl=0.5
ka=(5.0)
kb=0.01
calcMinPos=1
for att in ${att_arr[@]}; do
  for att2 in ${att2_arr[@]}; do
    for phi in ${phi_arr[@]}; do
      for t_stress in ${t_stress_arr[@]}; do
        for gamma in ${gamma_arr[@]}; do
          for k_on in ${kon_arr[@]}; do
            for k_off in ${koff_arr[@]}; do
              k_ecm=$att2
              #for k_ecm in ${kecm_arr[@]}; #do
              echo "./main/cell/psm2D.o   12  16 1.0 $phi $kl $ka $kb $att $att2 $t_stress    $v0    $tau_abp  $gamma $k_on $k_off $k_ecm $calcMinPos 1    50 testa_"$att"_a2_"$att2"_tm_"$t_stress"_p_"$phi"_t_"$tau_abp"_gamma_"$gamma"_k_on_"$k_on"_k_off_"$k_off"_k_ecm_"$k_ecm
              #done
            done
          done
        done
      done
    done
  done
done
*/
//                NCELLS NV  A0  phi att att2 t_maxwell_bd v0  tau_abp gamma k_on k_off k_ecm calcMinPos seed duration outFileStem

#include <sstream>
#include "cell.h"

#define NDIM 2
using namespace std;

// helper functions to parse arguments
template <typename T>
T parseArg(const std::string& arg) {
  T value;
  std::stringstream(arg) >> value;
  return value;
}

// global constants
const bool plotCompression = 0;     // whether or not to plot configuration during compression protocol (0 saves memory)
const double dphi0 = 0.005;         // packing fraction increment
const double kc = 1.0;              // interaction force spring constant (should be unit)
const double boxLengthScale = 2.5;  // neighbor list box size in units of initial l0
const double dt0 = 0.01;  // initial magnitude of time step in units of MD time
const double Ptol = 1e-5;
const double Ftol = 1e-3;
const double att_range = 0.3;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  double B = 1.0, phi0 = 0.7;
  //  Read command-line arguments into corresponding variables
  int NCELLS = parseArg<int>(argv[1]);
  int nv = parseArg<int>(argv[2]);
  double calA0 = parseArg<double>(argv[3]);
  double phi = parseArg<double>(argv[4]);
  phi0 = min(phi0, phi);
  double kl = parseArg<double>(argv[5]);
  double ka = parseArg<double>(argv[6]);
  double kb = parseArg<double>(argv[7]);
  double att = parseArg<double>(argv[8]);
  double att2 = parseArg<double>(argv[9]);
  double t_stress = parseArg<double>(argv[10]);
  double v0_abp = parseArg<double>(argv[11]);
  double tau_abp = parseArg<double>(argv[12]);
  double gamma = parseArg<double>(argv[13]);
  double k_on = parseArg<double>(argv[14]);
  double k_off = parseArg<double>(argv[15]);
  double k_ecm = parseArg<double>(argv[16]);
  bool calcMinPos = parseArg<int>(argv[17]);
  int seed = parseArg<int>(argv[18]);
  double runTime = parseArg<double>(argv[19]);
  std::string outFileStem = argv[20];

  int numCellTypes = 2;  // 0 = interior cell type (PSM) and 1 = exterior cell type (boundary)
  cell cell2D(NCELLS, 0.0, 0.0, seed, numCellTypes);
  cell2D.openFileStreams(outFileStem, calcMinPos);

  // set spring constants
  cell2D.setka(1.0);  // set to 1 for initialization and FIRE, then set to larger for dynamic
  cell2D.setkl(kl);
  cell2D.setkb(kb);
  cell2D.setkc(kc);
  cell2D.setkon(k_on);
  cell2D.setkoff(k_off);
  cell2D.setkecm(10 * k_ecm);
  cell2D.setB(B);
  if (t_stress > 0.0)
    cell2D.setMaxwellRelaxationTime(t_stress);  // t_stress is infinity unless this is uncommented
  cell2D.setpbc(0, false);                      //  specify non-periodic boundaries
  cell2D.setpbc(1, false);
  cell2D.setl1(att);  // set adhesion scales
  cell2D.setl2(att_range);
  assert(att < att_range && att2 < att_range);  // required to have a differentiable, finite adhesive potential
  assert(att < att_range / 2);                  // existing bug: l1 > l2/2 leads to purely repulsive, not sure why

  cout
      << "att, att2, att_range = " << att << '\t' << att2 << '\t' << att_range << '\n';

  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdate);
  dpmMemFn attractionWithActiveBrownianUpdate = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdateWithCrawling);
  dpmMemFn attractionSmoothActive = static_cast<void (dpm::*)()>(&cell::attractiveSmoothActive);
  dpmMemFn attractionSmoothActiveBrownianCatchBondsUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothActiveCatchBonds);
  dpmMemFn attractionSmoothActiveBrownianECMBondsUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothActiveECMBonds);
  dpmMemFn attractiveSmoothForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothForceUpdate);
  dpmMemFn attractiveSmoothWithPolyWalls = static_cast<void (dpm::*)()>(&cell::attractiveSmoothForceUpdateWithPolyWall);
  dpmMemFn attractivePolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveWithPolarityForceUpdate);
  dpmMemFn repulsiveForceUpdateWithPolyWalls = static_cast<void (dpm::*)()>(&cell::repulsiveForceUpdateWithPolyWall);
  dpmMemFn attractiveForceUpdateWithPolyWalls = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdateWithPolyWall);

  // cellTypeIntMat = ones to begin with
  for (int cellIDi = 0; cellIDi < numCellTypes; cellIDi++) {
    for (int cellIDj = 0; cellIDj < numCellTypes; cellIDj++) {
      // other than boundaries, off-diagonals are zero (only cells of same type interact)
      // if either cellID is a boundary, set it to 0 (no attraction with boundaries)
      // otherwise, one of the cellIDs must not be a boundary. If the cells are the same, they should attract (set to 1).
      double boundaryID = numCellTypes - 1;
      // when att > 0, use l1 and l1*modifier (att2/att * att)
      // when att = 0, use l1=0 and modifier (att2)
      double attractionRatio = (att > 0) ? att2 / att : att2;
      double attractionModifier =
          (cellIDi == boundaryID || cellIDj == boundaryID) ? attractionRatio : (cellIDi == cellIDj) ? 1.0
                                                                                                    : 0.0;
      cell2D.setCellTypeAttractionModifiers(cellIDi, cellIDj, attractionModifier);
    }
  }
  cell2D.printInteractionMatrix();

  // initialize particles with the same number of vertices and the same preferred shape parameter calA0
  cell2D.monodisperse2D(calA0, nv);
  int circleID = 0, rectangleID = 1, horseshoeID = 2;
  cell2D.initializeTransverseTissue(phi0, Ftol, circleID);  // initialize within a ring boundary
  cell2D.initializeNeighborLinkedList2D(boxLengthScale);

  // switch ring boundary for rectangular boundary
  // cell2D.replaceCircularBoundary(rectangleID, 2.0);

  // compress to desired density
  cell2D.resizeNeighborLinkedList2D();
  cell2D.shrinkPolyWall(attractiveSmoothWithPolyWalls, Ftol, dt0, phi, dphi0);
  // cell2D.shrinkPolyWall(repulsiveForceUpdateWithPolyWalls, Ftol, dt0, phi, dphi0);

  cell2D.replacePolyWallWithDP(numCellTypes);
  double restPackingFraction = 1.0;
  bool setPreferredShape = true;
  cell2D.resizeTissuePreferredArea(restPackingFraction, setPreferredShape);
  cell2D.resizeCatchBonds();
  cell2D.resizeNeighborLinkedList2D();

  double shortRelaxTime = 25.0;
  double relaxTime = 400.0;
  cell2D.setka(ka);
  cell2D.vertexDampedMD(attractiveSmoothForceUpdate, dt0, shortRelaxTime, 0.0);

  cell2D.shrinkCellVertices(attractiveSmoothForceUpdate, dt0, 0.0);

  cell2D.setl00();  // set l00 to be l0 before setting maxwell relaxation time
  cell2D.setMaxwellRelaxationTime(t_stress);
  cell2D.setActiveBrownianParameters(v0_abp, tau_abp);
  cell2D.resizeNeighborLinkedList2D();

  cell2D.setgamma(gamma);
  cell2D.vertexDampedMD(attractionSmoothActive, dt0, relaxTime, 0.0);

  // begin production run after all of the initialization and equilibration settles
  double v0_decay_rate = 0.0, v0_min = v0_abp;
  cout << "relaxing with final potential \n";
  cell2D.vertexDampedMD(attractionSmoothActiveBrownianECMBondsUpdate, dt0, relaxTime, 0.0, 0.0, v0_min);

  cout << "measuring with final potential!\n";
  cell2D.setitmax(1e4);
  cell2D.vertexDampedMD(attractionSmoothActiveBrownianECMBondsUpdate, dt0, runTime, 5.0, v0_decay_rate * v0_abp, v0_min);

  bool testingSteadyState = false;

  if (testingSteadyState) {
    cell2D.setl1(0);
    cell2D.vertexDampedMD(attractionSmoothActiveBrownianECMBondsUpdate, dt0, runTime, 5.0, v0_decay_rate * v0_abp, v0_min);

    cell2D.setl1(att);
    cell2D.vertexDampedMD(attractionSmoothActiveBrownianECMBondsUpdate, dt0, runTime, 5.0, v0_decay_rate * v0_abp, v0_min);
  }

  cout
      << "\n** Finished psm.cpp (2D transverse section of pre-somitic mesoderm), ending. " << endl;

  return 0;
}