// ~~~~ OpenFoam's FV headers
#include "fvCFD.H"

// ~~~~ Utility headers
#include <vector>
#include <string>

// ~~~~ Just for mpi timer:
#include "mpi.h"

int main( int argc, char *argv[] ) {
	#define NO_CONTROL
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	
	// ~~~~ Create / Read the input dictionary

    IOdictionary inputDict(
            IOobject (
                    "input",                                // dictionary name
                    runTime.constant(),             // dict is found in "constant"
                    mesh,                                   // registry for the dict
                    IOobject::MUST_READ,    // must exist, otherwise failure
                    IOobject::NO_WRITE              // dict is only read by the solver
            )
    );


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GET MESH QUANTITIES:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// mesh quantities:
	const vectorField& C = mesh.cellCentres(); 
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    // Get number of cells 	
    int nCells = mesh.nCells(); 
	    
    // Loop through the boundary objects to get the total number of boundary faces: 
    std::vector<int> nPatchFaces;
    std::vector<std::string> patchNames; 
    int nBoundaryFaces = 0; 
    forAll(mesh.boundary(), patch)
    {
        const word& patchName = mesh.boundary()[patch].name(); // name of boundary patch
        nPatchFaces.push_back(mesh.boundary()[patch].size());
        patchNames.push_back(patchName); 
        Info << "boundary patch name: " << patchNames[patch] << "\t, faces: " << nPatchFaces[patch] << endl;
        nBoundaryFaces += nPatchFaces[patch];
    }
    std::cout << "Number of internal faces: " << mesh.nInternalFaces() << std::endl;
    std::cout << "Number of boundary faces: " << nBoundaryFaces << std::endl;
    std::cout << "Number of cell centers: " << mesh.nCells() << std::endl;

    // Get the coordinates of the boundary in one list: 
    double Cf_b[nBoundaryFaces][3]; 
    int boundary_count = 0;
    forAll(mesh.boundary(), patch)
    {
        forAll(mesh.boundary()[patch], facei)
        {
            const label& face = boundaryMesh[patch].start() + facei; // global faceID
            
            for (int j = 0; j < 3; j++)
            {
                Cf_b[boundary_count][j] = boundaryMesh[patch].faceCentres()[facei][j];
            }
            boundary_count = boundary_count + 1;
        }
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GET DISTANCE FIELD (BRUTE FORCE):
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ~~~~ Initialize a distance field on the mesh: 
    volScalarField distance ( 
        IOobject("distance", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar( "distance", dimensionSet(0,1,0,0,0,0,0), 0.0));      
  
    // ~~~~ Compute distance. 
    Info << "Computing distance field..." << endl;
    double diff_temp = 0.;
    double dist_temp = 0.;
    double time_loop = MPI_Wtime(); 
    for (int cell = 0; cell < nCells; cell++)
    {
        // Initialize distance to some high value:
        distance.ref()[cell] = 1.0e20;

        // Loop through boundary
        for (int bface = 0; bface < nBoundaryFaces; bface++)
        {
            // Compute the distance from cell to boundary face: 
            dist_temp = 0; 
            for (int j = 0; j < 3; j++)
            {
                diff_temp = C[cell][j] - Cf_b[bface][j]; 
                dist_temp = dist_temp +  diff_temp * diff_temp; 
            }
            dist_temp = std::sqrt(dist_temp);

            // Update distance: 
            if (dist_temp < distance()[cell])
            {
                distance.ref()[cell] = dist_temp; 
            } 
        }
        // Info << "Distance for cell " << cell << " : " << distance()[cell] << endl;
    }
    time_loop = MPI_Wtime() - time_loop;
    Info << "Done. Total loop time was " << time_loop << " s." << endl;

    Info << "\nWriting..." << endl;
    distance.write();  
    Info << "Done." << endl;
}
