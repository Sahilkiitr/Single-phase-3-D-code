#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace plb;
typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

class PressureGradient {
public:
PressureGradient(T deltaP_, plint nx_) : deltaP(deltaP_), nx(nx_)
{ }
void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity)
{
velocity.resetToZero();
density = (T)1 - deltaP*DESCRIPTOR<T>::invCs2 / (T)(nx-1) * (T)iX;
}
private:
T deltaP;
plint nx;
};

void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField3D<int>&˓→geometry)
{
const plint nx = geometry.getNx();
const plint ny = geometry.getNy();
const plint nz = geometry.getNz();
Box3D sliceBox(0,0, 0,ny-1, 0,nz-1);
std::unique_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>→(geometry, sliceBox);
plb_ifstream geometryFile(fNameIn.c_str());
for (plint iX=0; iX<nx-1; ++iX) {
if (!geometryFile.is_open()) {
pcout << "Error: could not open geometry file " << fNameIn << std::endl;
exit(EXIT_FAILURE);
}
geometryFile >> *slice;
copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
VtkImageOutput3D<T> vtkOut("porousMedium", 1.0);
vtkOut.writeData<float>(*copyConvert<int,T>(geometry, geometry.˓→getBoundingBox()), "tag", 1.0);
  }

{
std::unique_ptr<MultiScalarField3D<T> > floatTags = copyConvert<int,T>˓→(geometry, geometry.getBoundingBox());
std::vector<T> isoLevels;
isoLevels.push_back(0.5);
typedef TriangleSet<T>::Triangle Triangle;
std::vector<Triangle> triangles;
Box3D domain = floatTags->getBoundingBox().enlarge(-1);
domain.x0++;
domain.x1--;
isoSurfaceMarchingCube(triangles, *floatTags, isoLevels, domain);
TriangleSet<T> set(triangles);
std::string outDir = fNameOut + "/";
set.writeBinarySTL(outDir + "porousMedium.stl");
}
}

void porousMediaSetup(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition,
MultiScalarField3D<int>& geometry, T deltaP)
{
const plint nx = lattice.getNx();
const plint ny = lattice.getNy();
const plint nz = lattice.getNz();
pcout << "Definition of inlet/outlet." << std::endl;
Box3D inlet (0,0, 1,ny-2, 1,nz-2);
boundaryCondition->addPressureBoundary0N(inlet, lattice);
setBoundaryDensity(lattice, inlet, (T) 1.);

Box3D outlet(nx-1,nx-1, 1,ny-2, 1,nz-2);
boundaryCondition->addPressureBoundary0P(outlet, lattice);
setBoundaryDensity(lattice, outlet, (T) 1. - deltaP*DESCRIPTOR<T>::invCs2);

pcout << "Definition of the geometry." << std::endl;
Where "geometry" evaluates to 1, use bounce-back.
defineDynamics(lattice, geometry, new BounceBack<T,DESCRIPTOR>(), 1);
defineDynamics(lattice, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);
pcout << "Initilization of rho and u." << std::endl;
initializeAtEquilibrium( lattice, lattice.getBoundingBox(),˓→PressureGradient(deltaP, nx) );
lattice.initialize();
delete boundaryCondition;

void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter)
 {
const plint nx = lattice.getNx();
const plint ny = lattice.getNy();
const plint nz = lattice.getNz();

const plint imSize = 600;
ImageWriter<T> imageWriter("leeloo");
imageWriter.writeScaledGif(createFileName("ux_inlet", iter, 6),
computeVelocityNorm(lattice, Box3D(0,0, 0,ny-1, 0,nz-1)),imSize, imSize );

imageWriter.writeScaledGif(createFileName("ux_half", iter, 6),
*computeVelocityNorm(lattice, Box3D(nx/2,nx/2, 0,ny-1, 0,nz-1)),imSize, imSize );
