// accepts cmdline parameters
// parameter 1 MTX file
// tolerance is parameter 2 optionally 


#include "AztecOO.h"
#include "AztecOO_Version.h"

/* The next two ifdefs should be removed when when the revision of
 az_aztec_defs.h is complete. */
#ifndef TRILINOS_NO_CONFIG_H
 
#ifdef HAVE_CSTDLIB
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#else
#include <cstdio>
#include <cstdlib>
#endif

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_BlockMapIn.h"
#include "Epetra_Export.h"

#define COUT std::cout << myPID << "> " <<
#define ENDL << std::endl

int myPID = -1;
int numProc = 1;
 
#define MAXELPR 100

int lclphase = -1, gblphase;
int readPID = 0;
bool distribMtx = false;

void syncnext( char *what, Epetra_MpiComm &Comm) {
	/*
	lclphase++;
	if (myPID == 0) COUT "Sync " << what ENDL;
	do {
		(void)Comm.MinAll(&lclphase, &gblphase, 1);
		if ( gblphase < lclphase )
			COUT "*" ENDL;
	} while (gblphase < lclphase);
	if (myPID == 0) COUT "Sync OK." ENDL;
	*/
}

char *MPIReadFileContents(const char *fn, Epetra_MpiComm &Comm) {
	MPI_File fh;		 
	if ( MPI_File_open(Comm.Comm(), fn, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS )
		return NULL;
	MPI_Offset size;		 
	MPI_File_get_size(fh, &size);
	char *buf = (char*)malloc(size+1);
	MPI_Status s;		 
	MPI_File_read(fh, (void*)buf, size, MPI_CHAR, &s);
	MPI_File_close(&fh);
	buf[size] = 0;
	return buf;
}

int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI
	// Initialize MPI
	MPI_Init(&argc, &argv);
	Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm Comm;
#endif
	int NumGlobalElements = 10000;
	double tol = 1E-5;
	 
	myPID = Comm.MyPID();
	numProc = Comm.NumProc();

	if (argc < 2) {
		COUT "Missing filename.";
		return 0;
	}
	const char *f = argv[1];

	if (argc > 2) {
		double t = strtof(argv[2], NULL);
		if (t > 0) tol = t;
	}


	std::cout << Comm << std::endl;
	if (myPID == 0)
		COUT AztecOO_Version() ENDL;

	syncnext("Init", Comm); // 1

	// Construct a Map that puts same number of equations on each processor
	Epetra_Map *PMap; // parallel map
	Epetra_Map *P0Map; // process 0 map, used when distributing matrix
	// Create a Epetra_Matrix
	Epetra_CrsMatrix *A=NULL,*M=NULL;
	Epetra_Vector *x, *b, *x0;

	////********************* CREATE MATRIX ****************************************

	// READ MTX FILE 
	if (!distribMtx || myPID == 0)
	{
/*		int phase0 = lclphase;
		for (int pn = 0; pn < numProc; pn++)
		{
			if (lclphase - phase0 == myPID)
			{
*/
				COUT "Read file (#rows)" ENDL;
				std::string fn = f;
				COUT "read row count " << fn ENDL;
				char *buf = MPIReadFileContents(fn.c_str(), Comm);
				if (buf == NULL || strlen(buf) == 0) {
					COUT "FAIL";
					return 0;
				}

				std::stringstream *buffer = new std::stringstream(buf);

				int n = 0;
				std::string s;
				std::getline(*buffer, s);
				*buffer >> NumGlobalElements;
				if (distribMtx)
					P0Map = new Epetra_Map(NumGlobalElements, NumGlobalElements, Comm);
				else
					PMap = new Epetra_Map(NumGlobalElements, 0, Comm);

				delete buffer;
				free(buf);
				COUT "OK read row count " << NumGlobalElements ENDL;

				//EpetraExt::MatrixMarketFileToMap(f,Comm,Map) doesntwork,expects strange data %numprocs etc.

				if (distribMtx)
					PMap = new Epetra_Map(NumGlobalElements, 0, Comm); // distributed map
				x = new Epetra_Vector(*PMap);
				b = new Epetra_Vector(*PMap);
				x0 = new Epetra_Vector(*PMap);

				Epetra_Map *Map = (distribMtx ? P0Map : PMap);
				COUT Map ENDL;
				{
					COUT "Read file (all)" ENDL;
					if (
						EpetraExt::MatrixMarketFileToCrsMatrix(f, *Map, A, false,
#ifdef _DEBUG 
							true
#else 
							false
#endif 
						) != 0)
					{
						COUT "Error reading " << f;
						return 0;
					}
					COUT "OK read file (all)" ENDL;
				}

				{
					std::string fn = f;
					fn = fn.substr(0, fn.find_last_of('.')) + ".rhs";
					COUT "read RHS" << fn ENDL;
					char *buf = MPIReadFileContents(fn.c_str(), Comm);
					if (buf == NULL || strlen(buf) == 0) {
						COUT "FAIL";
						return 0;
					}

					std::stringstream *buffer = new std::stringstream(buf);

					int n = 0;
					for (int i = 0; i < NumGlobalElements; i++) 
					{
						double f;
						std::string s;
						std::getline(*buffer, s);
						s.erase(s.find_last_not_of(" \n\r\t") + 1);
						std::replace(s.begin(), s.end(), ',', '.');
						f = stof(s);
						(*b)[i] = f;
						n++;
						if (buffer->eof())
							break;
					}
					delete buffer;
					free(buf);
					COUT "OK read RHS " << n ENDL;
				}

				{
					std::string fn = f;
					fn = fn.substr(0, fn.find_last_of('.')) + ".ans";
					COUT "read solution " << fn ENDL;
					char *buf = MPIReadFileContents(fn.c_str(), Comm);
					if (buf == NULL || strlen(buf) == 0) {
						COUT "FAIL";
						return 0;
					}
						
					std::stringstream *buffer = new std::stringstream(buf);

					int n = 0;
					for (int i = 0; i < NumGlobalElements; i++)
					{
						double f;
						std::string s;
						std::getline(*buffer, s);
						s.erase(s.find_last_not_of(" \n\r\t") + 1);
						std::replace(s.begin(), s.end(), ',', '.');
						f = stof(s);
						(*x0)[i] = f;
						n++;
						if (buffer->eof())
							break;
					}
					delete buffer;
					free(buf);
					COUT "OK read solution " << n ENDL;
				}

/*			}
			else
				COUT "Read file skipped now, trying int next turn." ENDL;
			
			if (distribMtx) // process 0 will send matrix and vectors to others...
				break;
*/			
						
//		}
	}
	syncnext("Read file", Comm);

	if ( distribMtx && myPID != 0) {
		Comm.MaxAll(&NumGlobalElements, &NumGlobalElements, 1);
		PMap = new Epetra_Map(NumGlobalElements, 0, Comm);
		x = new Epetra_Vector(*PMap);
		b = new Epetra_Vector(*PMap);
		x0 = new Epetra_Vector(*PMap);
	}
	int NumMyElements = PMap->NumMyElements();

	if (!distribMtx) {
		M = A;
	} else
	if (distribMtx) {
		syncnext("Before distribution", Comm);

		M = new Epetra_CrsMatrix(Copy, *PMap, 0);
		Epetra_Export exporter(*P0Map, *PMap);
		M->Export(*A, exporter, Insert);
	}

	// some info
	if (myPID == 0) {
		COUT "#elements " << NumGlobalElements ENDL;
		COUT "#myelements " << NumMyElements ENDL;
		COUT " tolerance " << tol ENDL;
	}
	 

	COUT "START SETUP MATRIX AND SOLVE" ENDL;
	std::cout.flush();
	//_sleep(500);


	M->FillComplete();

	syncnext("Fill complete", Comm); 
							
	// Create Linear Problem
	Epetra_LinearProblem *P;
	P = new Epetra_LinearProblem(M, x, b);
	// Create AztecOO instance 
	AztecOO *solver;
	solver = new AztecOO(*P); 

	//int n = solver.GetAztecOption(AZ_precond); // AZ_dom_decomp
	solver->SetAztecOption(AZ_precond, AZ_Jacobi);
	//n = solver.GetAztecOption(AZ_solver); // AZ_gmres
	solver->SetAztecOption(AZ_solver, AZ_gmres); //! AZ_cg only for symmetric!

	syncnext("SolverInit", Comm); 

	solver->Iterate(10000, tol); 

	syncnext("Solve", Comm); // 4

	if (myPID == 0) {
		COUT "Solver performed " << solver->NumIters() << " iterations." ENDL
			<< "Norm of true residual = " << solver->TrueResidual() ENDL
			<< "Time " << solver->SolveTime() ENDL;

		// compare x with x0
		double diff = 0;
		for (int i = 0; i < NumGlobalElements; i++)
		{
			if (i < 10)
				COUT "x: " << (*x)[i] << " x0: " << (*x0)[i] ENDL;
			if (i == 10)
				COUT "..." ENDL;
			double x1 = (*x0)[i];
			double x2 = (*x)[i];
			diff += abs(x2 - x1);
		}
		COUT "Sum diff X to X0 vect: " << diff << " avg: " << diff / NumGlobalElements ENDL;
	}

	delete A;
	if ( distribMtx )
		delete M;
	delete x;
	delete b;
	delete x0;

//END:
#ifdef EPETRA_MPI
	MPI_Finalize();
#endif

	// getch();

	return 0;
}

