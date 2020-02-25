#pragma once

template<typename int_type>
static void quickpart_list_inc_int(
	int_type *list, int_type *parlista, double *parlistb,
	int start, int end, int *equal, int *larger)
{
	int i;
	int_type key, itmp;
	double dtmp;

	key = list ? list[(end + start) / 2] : 1;

	*equal = *larger = start;
	for (i = start; i <= end; i++)
		if (list[i] < key) {
			itmp = parlista[i];
			parlista[i] = parlista[*larger];
			parlista[(*larger)] = parlista[*equal];
			parlista[(*equal)] = itmp;
			dtmp = parlistb[i];
			parlistb[i] = parlistb[*larger];
			parlistb[(*larger)] = parlistb[*equal];
			parlistb[(*equal)] = dtmp;
			itmp = list[i];
			list[i] = list[*larger];
			list[(*larger)++] = list[*equal];
			list[(*equal)++] = itmp;
		}
		else if (list[i] == key) {
			itmp = parlista[i];
			parlista[i] = parlista[*larger];
			parlista[(*larger)] = itmp;
			dtmp = parlistb[i];
			parlistb[i] = parlistb[*larger];
			parlistb[(*larger)] = dtmp;
			list[i] = list[*larger];
			list[(*larger)++] = key;
		}
}

template<typename int_type>
static void sort_three(int_type* list, int_type *parlista, double *parlistb,
	int start, int end)
{
	int  equal, larger;

	if (start < end) {
		quickpart_list_inc_int(list, parlista, parlistb, start, end,
			&equal, &larger);
		sort_three(list, parlista, parlistb, start, equal - 1);
		sort_three(list, parlista, parlistb, larger, end);
	}
}

int MatrixMarketFileToCrsMatrixHandle(const char *filename,
	const Epetra_Comm & comm,
	Epetra_CrsMatrix * & A,
	const Epetra_Map * rowMap,
	const Epetra_Map * colMap,
	const Epetra_Map * rangeMap,
	const Epetra_Map * domainMap,
	const bool transpose,
	const bool verbose
)
{
	const int chunk_read = 500000;  //  Modify this variable to change the
									//  size of the chunks read from the file.
	const int headerlineLength = 257;
	const int lineLength = 81;
	const int tokenLength = 35;
	char line[headerlineLength];
	char token1[tokenLength];
	char token2[tokenLength];
	char token3[tokenLength];
	char token4[tokenLength];
	char token5[tokenLength];
	int M, N, NZ;      // Matrix dimensions
	int me = comm.MyPID();

	Epetra_Time timer(comm);

	// Make sure domain and range maps are either both defined or undefined
	if ((domainMap != 0 && rangeMap == 0) || (domainMap == 0 && rangeMap != 0)) {
		EPETRA_CHK_ERR(-3);
	}

	// check maps to see if domain and range are 1-to-1

	if (domainMap != 0) {
		if (!domainMap->UniqueGIDs()) { EPETRA_CHK_ERR(-2); }
		if (!rangeMap->UniqueGIDs()) { EPETRA_CHK_ERR(-2); }
	}
	else {
		// If domain and range are not specified,
		// rowMap becomes both and must be 1-to-1 if specified
		if (rowMap != 0) {
			if (!rowMap->UniqueGIDs()) { EPETRA_CHK_ERR(-2); }
		}
	}

	FILE * handle = 0;
	if (me == 0) {
		if (verbose) std::cout << "Reading MatrixMarket file " << filename << std::endl;
		handle = fopen(filename, "r");  // Open file
		if (handle == 0)
			EPETRA_CHK_ERR(-1); // file not found

								// Check first line, which should be
								// %%MatrixMarket matrix coordinate real general
		if (fgets(line, headerlineLength, handle) == 0) {
			if (handle != 0) fclose(handle);
			EPETRA_CHK_ERR(-1);
		}
		if (sscanf(line, "%s %s %s %s %s", token1, token2, token3, token4, token5) == 0) {
			if (handle != 0) fclose(handle);
			EPETRA_CHK_ERR(-1);
		}
		if (strcmp(token1, "%%MatrixMarket") ||
			strcmp(token2, "matrix") ||
			strcmp(token3, "coordinate") ||
			strcmp(token4, "real") ||
			strcmp(token5, "general")) {
			if (handle != 0) fclose(handle);
			EPETRA_CHK_ERR(-1);
		}

		// Next, strip off header lines (which start with "%")
		do {
			if (fgets(line, headerlineLength, handle) == 0) {
				if (handle != 0) fclose(handle);
				EPETRA_CHK_ERR(-1);
			}
		} while (line[0] == '%');

		// Next get problem dimensions: M, N, NZ
		if (sscanf(line, "%d %d %d", &M, &N, &NZ) == 0) {
			if (handle != 0) fclose(handle);
			EPETRA_CHK_ERR(-1);
		}
	}
	comm.Broadcast(&M, 1, 0);
	comm.Broadcast(&N, 1, 0);
	comm.Broadcast(&NZ, 1, 0);

	// Now create matrix using user maps if provided.


	// Now read in chunks of triplets and broadcast them to other processors.
	char *buffer = new char[chunk_read*lineLength];
	int nchunk;
	int nmillion = 0;
	int nread = 0;
	int rlen;

	// Storage for this processor's nonzeros.
	const int localblock = 100000;
	int localsize = NZ / comm.NumProc() + localblock;
	int *iv = (int *)malloc(localsize * sizeof(int));
	int *jv = (int *)malloc(localsize * sizeof(int));
	double *vv = (double *)malloc(localsize * sizeof(double));
	int lnz = 0;   //  Number of non-zeros on this processor.

	if (!iv || !jv || !vv)
		EPETRA_CHK_ERR(-1);

	Epetra_Map *rowMap1;
	bool allocatedHere = false;
	if (rowMap != 0)
		rowMap1 = const_cast<Epetra_Map *>(rowMap);
	else {
		rowMap1 = new Epetra_Map(M, 0, comm);
		allocatedHere = true;
	}
	int ioffset = rowMap1->IndexBase() - 1;
	int joffset = (colMap != 0 ? colMap->IndexBase() - 1 : ioffset);

	int rowmajor = 1;  // Assume non-zeros are listed in row-major order,
	int prevrow = -1;  // but test to detect otherwise.  If non-zeros
					   // are row-major, we can avoid the sort.

	while (nread < NZ) {
		if (NZ - nread > chunk_read) nchunk = chunk_read;
		else nchunk = NZ - nread;

		if (me == 0) {
			char *eof;
			rlen = 0;
			for (int i = 0; i < nchunk; i++) {
				eof = fgets(&buffer[rlen], lineLength, handle);
				if (eof == NULL) {
					fprintf(stderr, "%s", "Unexpected end of matrix file.");
					EPETRA_CHK_ERR(-1);
				}
				rlen += strlen(&buffer[rlen]);
			}
			buffer[rlen++] = '\n';
		}
		comm.Broadcast(&rlen, 1, 0);
		comm.Broadcast(buffer, rlen, 0);

		buffer[rlen++] = '\0';
		nread += nchunk;

		// Process the received data, saving non-zeros belonging on this proc.
		char *lineptr = buffer;
		for (rlen = 0; rlen < nchunk; rlen++) {
			char *next = strchr(lineptr, '\n');
			int I = atoi(strtok(lineptr, " \t\n")) + ioffset;
			int J = atoi(strtok(NULL, " \t\n")) + joffset;
			double V = atof(strtok(NULL, " \t\n"));
			lineptr = next + 1;
			if (transpose) {
				// swap I and J indices.
				int tmp = I;
				I = J;
				J = tmp;
			}
			if (rowMap1->MyGID(I)) {
				//  This processor keeps this non-zero.
				if (lnz >= localsize) {
					// Need more memory to store nonzeros.
					localsize += localblock;
					iv = (int *)realloc(iv, localsize * sizeof(int));
					jv = (int *)realloc(jv, localsize * sizeof(int));
					vv = (double *)realloc(vv, localsize * sizeof(double));
				}
				iv[lnz] = I;
				jv[lnz] = J;
				vv[lnz] = V;
				lnz++;
				if (I < prevrow) rowmajor = 0;
				prevrow = I;
			}
		}

		// Status check
		if (nread / 1000000 > nmillion) {
			nmillion++;
			if (verbose && me == 0) std::cout << nmillion << "M ";
		}
	}

	delete[] buffer;

	if (!rowmajor) {
		// Sort into row-major order (by iv) so can insert entire rows at once.
		// Reorder jv and vv to parallel iv.
		if (verbose && me == 0) std::cout << std::endl << "   Sorting local nonzeros" << std::endl;
		sort_three(iv, jv, vv, 0, lnz - 1);
	}

	//  Compute number of entries per local row for use in constructor;
	//  saves reallocs in FillComplete.

	//  Now construct the matrix.
	//  Count number of entries in each row so can use StaticProfile=2.
	if (verbose && me == 0) std::cout << std::endl << "   Constructing the matrix" << std::endl;
	int numRows = rowMap1->NumMyElements();
	int *numNonzerosPerRow = new int[numRows];
	for (int i = 0; i < numRows; i++) numNonzerosPerRow[i] = 0;
	for (int i = 0; i < lnz; i++)
		numNonzerosPerRow[rowMap1->LID(iv[i])]++;

	if (rowMap != 0 && colMap != 0)
		A = new Epetra_CrsMatrix(Copy, *rowMap, *colMap, numNonzerosPerRow);
	else if (rowMap != 0) {
		// Construct with StaticProfile=true since we know numNonzerosPerRow.
		// Less memory will be needed in FillComplete.
		A = new Epetra_CrsMatrix(Copy, *rowMap, numNonzerosPerRow, true);
	}
	else {
		// Construct with StaticProfile=true since we know numNonzerosPerRow.
		// Less memory will be needed in FillComplete.
		A = new Epetra_CrsMatrix(Copy, *rowMap1, numNonzerosPerRow, true);
	}
	A->SetTracebackMode(1);

	// Rows are inserted in ascending global number, and the insertion uses numNonzerosPerRow.
	// Therefore numNonzerosPerRow must be permuted, as it was created in ascending local order.
	int *gidList = new int[numRows];
	for (int i = 0; i < numRows; i++) gidList[i] = rowMap1->GID(i);
	Epetra_Util::Sort(true, numRows, gidList, 0, 0, 1, &numNonzerosPerRow);
	delete[] gidList;

	//  Insert the global values into the matrix row-by-row.
	if (verbose && me == 0) std::cout << "   Inserting global values" << std::endl;
	{
		int i = 0;
		for (int sum = 0; i < numRows; i++) {
			if (numNonzerosPerRow[i]) {
				int ierr = A->InsertGlobalValues(iv[sum], numNonzerosPerRow[i],
					&vv[sum], &jv[sum]);
				if (ierr<0) EPETRA_CHK_ERR(ierr);
				sum += numNonzerosPerRow[i];
			}
		}
	}

	delete[] numNonzerosPerRow;
	free(iv);
	free(jv);
	free(vv);

	if (verbose && me == 0) std::cout << "   Completing matrix fill" << std::endl;
	if (rangeMap != 0 && domainMap != 0) {
		EPETRA_CHK_ERR(A->FillComplete(*domainMap, *rangeMap));
	}
	else if (M != N) {
		Epetra_Map newDomainMap(N, rowMap1->IndexBase(), comm);
		EPETRA_CHK_ERR(A->FillComplete(newDomainMap, *rowMap1));
	}
	else {
		EPETRA_CHK_ERR(A->FillComplete());
	}

	if (allocatedHere) delete rowMap1;

	if (handle != 0) fclose(handle);
	double dt = timer.ElapsedTime();
	if (verbose && me == 0) std::cout << "File Read time (secs):  " << dt << std::endl;
	return(0);
}

int MatrixMarketFileToCrsMatrix(const char *filename,
	const Epetra_Map & rowMap, Epetra_CrsMatrix * & A, const bool transpose,
	const bool verbose)
{
	EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename,
		rowMap.Comm(), A,
		&rowMap, 0, 0, 0,
		transpose, verbose));
	return(0);
}