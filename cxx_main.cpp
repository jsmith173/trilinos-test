// reads parameter 1 as MTX file, RHS file with same name, but RHS extension, and ANS as expected answer if they are there
// tolerance is parameter 2 optionally then

//#include <conio.h>

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
//#include "EpetraExt_CrsMatrixIn.h"
//#include "EpetraExt_BlockMapIn.h"
#include "Epetra_Time.h"
#include "MTXReader.h"

#define COUT std::cout << myPID << "> " <<
#define ENDL << std::endl; std::cout.flush()

#define EXIT MPI_Finalize(); return 1;

// I tried putting sync. points (BARRIER) to the program but that does not make any difference...
//#define BARRIER MPI_Barrier(Comm.Comm())

char *MPIReadFileContents(const char *fn, Epetra_MpiComm &Comm) {
	MPI_File fh;
	if (MPI_File_open(Comm.Comm(), fn, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
		return NULL;
	MPI_Offset size;
	MPI_File_get_size(fh, &size);
	char *buf = (char*)malloc(size + 1);
	MPI_Status s;
	MPI_File_read(fh, (void*)buf, size, MPI_CHAR, &s);
	MPI_File_close(&fh);
	buf[size] = 0;
	return buf;
}

int main(int argc, char *argv[])
{
	int myPID = -1;

#ifdef EPETRA_MPI
	// Initialize MPI
	MPI_Init(&argc, &argv);
	Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm Comm;
#endif
	int NumGlobalElements;

	double tol = 1E-5;

	if (argc < 2) {
		std::cout << "Missing filename.";
		EXIT;
	}

	if (argc > 2) {
		double t = strtof(argv[2], NULL);
		if (t > 0) tol = t;
	}

	myPID = Comm.MyPID();
	if (myPID == 0)
		std::cout << AztecOO_Version() << std::endl << std::endl;

	std::cout << Comm << std::endl;
	
	// Construct a Map that puts same number of equations on each processor
	Epetra_Map *Map;
	// Create a Epetra_Matrix
	Epetra_CrsMatrix *A;
	Epetra_Vector *x, *b, *x0;

////********************* CREATE MATRIX ****************************************

	const char *f = argv[1];

	COUT "Read file (#rows)" ENDL;
	std::string fn = f;
	COUT "read row count " << fn ENDL;
	char *buf = MPIReadFileContents(fn.c_str(), Comm);
	if (buf == NULL || strlen(buf) == 0) {
		COUT "FAIL";
		EXIT;
	}

	std::stringstream *buffer = new std::stringstream(buf);

	std::string s;
	std::getline(*buffer, s);
	*buffer >> NumGlobalElements;
	COUT "NumGlobalElements: " << NumGlobalElements ENDL;
	// THE MAP:
	Map = new Epetra_Map(NumGlobalElements, 0, Comm);

	delete buffer;
	free(buf);

	//EpetraExt::MatrixMarketFileToMap(f,Comm,Map) doesntwork,expects strange data %numprocs etc.

	x = new Epetra_Vector(*Map);
	x->PutScalar(0);
	b = new Epetra_Vector(*Map);
	b->PutScalar(0);
	x0 = new Epetra_Vector(*Map);
	x0->PutScalar(0);

	COUT Map ENDL;
	{
		COUT "Read file (all)" ENDL;
		if (
			//EpetraExt::MatrixMarketFileToCrsMatrix(f, *Map, A, false,
			MatrixMarketFileToCrsMatrix(f,  *Map, A, false,
#ifdef _DEBUG 
				true
#else 
				false
#endif 
			) != 0)
		{
			COUT "Error reading " << f;
			EXIT;
		}
		COUT "OK read file (all)" ENDL;
	}

	int NumMyElements = Map->NumMyElements();
	int MyMapIdx = Map->MinMyGID();
	COUT "My MAPGID " << MyMapIdx ENDL;
	{
		std::string fn = f;
		fn = fn.substr(0, fn.find_last_of('.')) + ".rhs";
		COUT "read RHS" << fn ENDL;
		char *buf = MPIReadFileContents(fn.c_str(), Comm);
		if (buf == NULL || strlen(buf) == 0) {
			COUT "FAIL";
			EXIT;
		} 
		COUT "file read OK " << fn ENDL;

		std::stringstream *buffer = new std::stringstream(buf);

		int n = NumMyElements, cnt=0;
		for (int i = 0; i < NumGlobalElements; i++)
			if ( i >= MyMapIdx) {
				double f;
				std::string s;
				std::getline(*buffer, s);
				s.erase(s.find_last_not_of(" \n\r\t") + 1);
				std::replace(s.begin(), s.end(), ',', '.');
				f = stof(s);
				(*b)[cnt++] = f;
				n--;
				if (n==0 || buffer->eof())
					break;
			}
		delete buffer;
		free(buf);
		COUT "OK read RHS " << cnt << " +" << MyMapIdx ENDL;
	}

	bool hasAns = false;
	{
		std::string fn = f;
		fn = fn.substr(0, fn.find_last_of('.')) + ".ans";
		COUT "read solution " << fn ENDL;
		char *buf = MPIReadFileContents(fn.c_str(), Comm);
		if (buf == NULL || strlen(buf) == 0) {
			COUT "Solution file missing" ENDL;
		}
		else {
			COUT "file read OK " << fn ENDL;

			std::stringstream *buffer = new std::stringstream(buf);

			int n = NumMyElements, cnt=0;
			for (int i = 0; i < NumGlobalElements; i++)
				if (i >= MyMapIdx) {
					double f;
					std::string s;
					std::getline(*buffer, s);
					s.erase(s.find_last_not_of(" \n\r\t") + 1);
					std::replace(s.begin(), s.end(), ',', '.');
					f = stof(s);
					(*x0)[cnt++] = f;
					n--;
					if (n==0 || buffer->eof())
						break;
				}
			delete buffer;
			free(buf);
			hasAns = true;
			COUT "OK read solution " << cnt << " +" << MyMapIdx ENDL;
		}
	}

	// some info
	if (myPID == 0) {
		std::cout << "#elements " << NumGlobalElements << std::endl;
		std::cout << "#myelements " << NumMyElements << std::endl;
		std::cout << " tolerance " << tol << std::endl;
	}
	 
	COUT "Finish up matrix" ENDL;
	A->FillComplete();

	// Create Linear Problem
	Epetra_LinearProblem problem(A, x, b);
	// Create AztecOO instance 
	AztecOO solver(problem);

	//int n = solver.GetAztecOption(AZ_precond); // AZ_dom_decomp
	solver.SetAztecOption(AZ_precond, AZ_Jacobi);
	//n = solver.GetAztecOption(AZ_solver); // AZ_gmres
	solver.SetAztecOption(AZ_solver, AZ_gmres); //! AZ_cg only for symmetric!

	COUT "Solve..." ENDL;
	solver.Iterate(100, tol);
	
	if (myPID == 0) 
		std::cout << "Solver performed " << solver.NumIters() << " iterations." << std::endl
			<< "Norm of true residual = " << solver.TrueResidual() << std::endl
			<< "Time " << solver.SolveTime() << std::endl;

	// compare with pre-read solution
	if (!hasAns) {
		COUT "*** No solution file found to compare with." ENDL;
	} else {
		double diff = 0;
		// compare only our own part of the vector
		for (int i = 0; i < NumMyElements; i++)
		{
			int idx = i + MyMapIdx;
			if (i < 10)
				COUT "x[" << idx << "]: " << (*x)[i] << " x0[" << idx << "]: " << (*x0)[i] ENDL;
			if (i == 10) {
				COUT "..." ENDL; break;
			}

			double x1 = (*x0)[i];
			double x2 = (*x)[i];
			diff += abs(x2 - x1);
		}
		COUT "Sum diff X to X0 vect (indexes " << MyMapIdx << ".." << MyMapIdx+NumMyElements-1 << "): " << diff << " avg: " << diff / NumGlobalElements << std::endl;
	}

#ifdef _WINDOWS
	_beep(1000, 200);
#endif

#ifdef EPETRA_MPI
	MPI_Finalize();
#endif

	// getch();

	return 0;
}

