// ~~~~ OpenFoam's FV headers
#include "fvCFD.H"

// ~~~~ Utility headers
#include <vector>
#include <string>
// #include "time.h"

int main( int argc, char *argv[] ) {
	#define NO_CONTROL
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	
	// ~~~~ Create / Read the input dictionary

    // IOdictionary inputDict(
    //         IOobject (
    //                 "input",                                // dictionary name
    //                 runTime.constant(),             // dict is found in "constant"
    //                 mesh,                                   // registry for the dict
    //                 IOobject::MUST_READ,    // must exist, otherwise failure
    //                 IOobject::NO_WRITE              // dict is only read by the solver
    //         )
    // );


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GET BOUNDARY FACES:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Info << "Getting mesh properties..." << endl;
    // ~~~~ mesh quantities:
    const vectorField& C = mesh.cellCentres(); 
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    // ~~~~ Get number of cells 	
    int nCells = mesh.nCells(); 
            
    // ~~~~ Loop through the boundary objects to get the total number of boundary faces: 
    int nBoundaryFaces = 0; 
    
    forAll(mesh.boundary(), patch)
    {
        const word& patchName = mesh.boundary()[patch].name(); // name of boundary patch
        // substring "proc" should not be found in boundary patch name 
        if ((patchName.find("proc") == std::string::npos) && 
            (mesh.boundary()[patch].size() > 0))
        {     
            nBoundaryFaces += mesh.boundary()[patch].size();
        }
    }
    
    // ~~~~ Check to see if any processors are not neighboring the geometry: 
    int no_boundary_sum = 0; 
    if (nBoundaryFaces == 0)
    {
        no_boundary_sum = 1; 
        nBoundaryFaces = 1; // add dummy boundary face 
    }
    
    // ~~~~ How many processors are disconnected?
    int no_boundary = no_boundary_sum;
    reduce(no_boundary_sum, sumOp<scalar>());
    Info << no_boundary_sum << " processor[s] not connected to geometry boundary." << endl;  

    // ~~~~ Do an allreduce on the total number of boundary faces.
    int nBoundaryFaces_ = nBoundaryFaces; 
    reduce(nBoundaryFaces, sumOp<scalar>());
    Info << "Total number of boundary faces on geometry: " << nBoundaryFaces - no_boundary_sum << endl; 

    // ~~~~ Get processor boundary coordinates: 
    List<scalar> Cf_b_x(nBoundaryFaces_);
    List<scalar> Cf_b_y(nBoundaryFaces_);
    List<scalar> Cf_b_z(nBoundaryFaces_);

    int boundary_count = 0;
    
    // ~~~~ Take care of dummy boundaries:
    if (no_boundary == 1)
    {
        Cf_b_x[boundary_count] = 1.0e20;
        Cf_b_y[boundary_count] = 1.0e20;
        Cf_b_z[boundary_count] = 1.0e20;
        boundary_count = boundary_count + 1;
    }
   
    // ~~~~ Loop through actual domain boundaries:  
    forAll(mesh.boundary(), patch)
    {   
        const word& patchName = mesh.boundary()[patch].name(); // name of boundary patch
        if ((patchName.find("proc") == std::string::npos) && 
            (mesh.boundary()[patch].size() > 0) && 
            (no_boundary == 0)) // if "proc" is not in patch name, and boundary mesh has faces
        {
            forAll(mesh.boundary()[patch], facei)
            {   
                // Store processor-local coordinates: 
                Cf_b_x[boundary_count] = boundaryMesh[patch].faceCentres()[facei][0];
                Cf_b_y[boundary_count] = boundaryMesh[patch].faceCentres()[facei][1];
                Cf_b_z[boundary_count] = boundaryMesh[patch].faceCentres()[facei][2];
                    
                boundary_count = boundary_count + 1;
            }
        }
        // else if ((no_boundary == 1) && (boundary_count < 1)) // if this is domain region disconnected from geometry
        // {
        //     // Update dummy boundary value with very high number. 
        //     for (int i = 0; i < 1; i++)
        //     {
        //         Cf_b_x[boundary_count] = 1.0e20; 
        //         Cf_b_y[boundary_count] = 1.0e20;
        //         Cf_b_z[boundary_count] = 1.0e20;
        //         boundary_count = boundary_count + 1;
        //     } 
        // }  
    }
    
    if (nBoundaryFaces_ != boundary_count)
    {
        std::cout << "Mismatch in available boundaries for " << Pstream::myProcNo() 
                    << " --- aborting." << endl; 
        exit(1);
    }
   
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GATHER LOCAL INFORMATION ON ALL PROCS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Info << "Gathering local info on all processors..." << endl; 
    // Gather the number of boundary faces: 
    List<scalar> nBoundaryFaces_list(Pstream::nProcs()); 
    nBoundaryFaces_list[Pstream::myProcNo()] = nBoundaryFaces_; 
    Pstream::gatherList(nBoundaryFaces_list);
    Pstream::scatterList(nBoundaryFaces_list);

    // Gather the boundary coordinates:  
    List < List<scalar> > CF_B_X(Pstream::nProcs());
    CF_B_X[Pstream::myProcNo()] = Cf_b_x; 
    Pstream::gatherList(CF_B_X); 
    Pstream::scatterList(CF_B_X);

    List < List<scalar> > CF_B_Y(Pstream::nProcs());
    CF_B_Y[Pstream::myProcNo()] = Cf_b_y; 
    Pstream::gatherList(CF_B_Y);
    Pstream::scatterList(CF_B_Y);

    List < List<scalar> > CF_B_Z(Pstream::nProcs());
    CF_B_Z[Pstream::myProcNo()] = Cf_b_z; 
    Pstream::gatherList(CF_B_Z);
    Pstream::scatterList(CF_B_Z);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GET DISTANCE FIELD (BRUTE FORCE):
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // ~~~~ Initialize a distance field on the mesh: 
    volScalarField distance ( 
        IOobject("distance", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar( "distance", dimensionSet(0,1,0,0,0,0,0), 0.0));  

    // ~~~~ Compute distance: 
    Info << "Computing distance field..." << endl; 
    double diff_temp_x = 0.; 
    double diff_temp_y = 0.;
    double diff_temp_z = 0.; 
    double dist_temp = 0.; 
    clock_t start, end;

    start = std::clock();
    for (int cell = 0; cell < nCells; cell++)
    {   
        // Initialize distance to some high value:
        distance.ref()[cell] = 1.0e20;

        // Boundary: Loop through processors: 
        for (int proc = 0; proc < Pstream::nProcs(); proc++)
        {
            // Boundary: Loop through local boundary:
            for (int bface = 0; bface < nBoundaryFaces_list[proc]; bface++)
            {
                // Compute the distance from cell to boundary face: 
                diff_temp_x = C[cell][0] - CF_B_X[proc][bface];
                diff_temp_y = C[cell][1] - CF_B_Y[proc][bface];
                diff_temp_z = C[cell][2] - CF_B_Z[proc][bface];
               
                dist_temp = diff_temp_x*diff_temp_x + 
                            diff_temp_y*diff_temp_y + 
                            diff_temp_z*diff_temp_z; 

                dist_temp = std::sqrt(dist_temp); 
                
                // Update distance:
                if (dist_temp < distance()[cell])
                {
                    distance.ref()[cell] = dist_temp;
                }
            }
        } 
    }
    end = std::clock();
    double time_taken = double(end - start)/double(CLOCKS_PER_SEC);
    reduce(time_taken, maxOp<scalar>());
    Info << "Done in " << time_taken << " sec" << endl;
    
    Info << "\nWriting..." << endl;
    distance.write();  
    Info << "Done." << endl; 
}
