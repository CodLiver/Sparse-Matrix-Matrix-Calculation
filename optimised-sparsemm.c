#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
	const COO, const COO, const COO,
	COO *);
void CSRmaker(struct sortCOO *list1, double *mtxCSRnz, int *mtxCSRn, int *mtxCSRm, int *Anz);

void matrixAdder(int *curr, double *C2Vec, int *mtxAn, double *mtxAnz, int *mtxAm);

void mtxGiantAdder(COO A, COO B, COO C, COO *spC);

void dataFiller(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC);
void dataFillerv2(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC);

int qsorter(const void *i, const void *j);


void CSRmaker(struct sortCOO *list1, double *mtxCSRnz, int *mtxCSRn, int *mtxCSRm, int *Anz) {
	for (int i = 0; i < *Anz; i++) {
		mtxCSRnz[i] = list1[i].d;
		mtxCSRn[i] = list1[i].j;
		mtxCSRm[list1[i].i+1]=i+1;
	}

}

void matrixAdder(int *curr, double *C2Vec, int *mtxAn, double *mtxAnz, int *mtxAm) {

	int getARowSize = mtxAm[*curr] - mtxAm[*curr - 1];
	int curPos = mtxAm[*curr - 1];
	int tmpPos = 0;


	for (int row = 0; row < getARowSize; ++row)
	{
		tmpPos = curPos + row;
		C2Vec[mtxAn[tmpPos]] += mtxAnz[tmpPos];
	}

}

/*void mtxGiantAdder(COO A, COO B, COO C, COO *spC) {

	int Anz, Bnz, Cnz;

	Anz = A->NZ, Bnz = B->NZ, Cnz = C->NZ;

	double *mtxAnz = malloc(Anz * sizeof(double));//has nz vals
	int *mtxAm = malloc((A->m + 1) * sizeof(int));//has cumulative els
	mtxAm[0] = 0;
	mtxAm[A->m] = Anz;
	int *mtxAn = malloc(Anz * sizeof(int));//col index
	//
	double *mtxBnz = malloc(Bnz * sizeof(double));
	int *mtxBm = malloc((B->m + 1) * sizeof(int));
	mtxBm[0] = 0;
	mtxBm[B->m] = Bnz;
	int *mtxBn = malloc(Bnz * sizeof(int));
	//
	double *mtxCnz = malloc(Cnz * sizeof(double));
	int *mtxCm = malloc((C->m + 1) * sizeof(int));
	mtxCm[0] = 0;
	mtxCm[C->m] = Cnz;
	int *mtxCn = malloc(Cnz * sizeof(int));

	double *mtxResAnz = malloc((Cnz + Bnz + Anz) * sizeof(double));
	int *mtxResAm = malloc((Cnz + Bnz + Anz) * sizeof(int));
	int *mtxResAn = malloc((Cnz + Bnz + Anz) * sizeof(int));

	//fill CSR matrix
#pragma acc parallel
	{
		CSRmaker(A, &Anz, mtxAnz, mtxAn, mtxAm);//NZ here
		CSRmaker(B, &Bnz, mtxBnz, mtxBn, mtxBm);
		CSRmaker(C, &Cnz, mtxCnz, mtxCn, mtxCm);
	}
	printf("MM starts \n");

	double *C2Vec = calloc(C->n, sizeof(double));//vectorized format
	int posC = 0;

	//#pragma acc kernels loop
	for (int cur = 1; cur < C->m + 1; ++cur, ++cur)
	{
#pragma acc parallel
		{
			matrixAdder(&cur, C2Vec, mtxAn, mtxAnz, mtxAm);
			matrixAdder(&cur, C2Vec, mtxBn, mtxBnz, mtxBm);
			matrixAdder(&cur, C2Vec, mtxCn, mtxCnz, mtxCm);
		}
		--cur;


		for (int adder = 0; adder < C->n; ++adder)//put indice here.
		{
			if (C2Vec[adder])
			{
				mtxResAnz[posC] = C2Vec[adder];//data indice here
				mtxResAn[posC] = adder;//col indice
				mtxResAm[posC] = cur;
				++posC;
			}
		}


		memset(C2Vec, 0, C->n * sizeof(double));
	}

	alloc_sparse(C->m, C->n, posC, spC);//??

	dataFiller(*spC, mtxResAn, mtxResAnz, mtxResAm, &posC);

	//print_sparse(*spC);

	free(C2Vec);
	free(mtxResAm);
	free(mtxResAn);
	free(mtxResAnz);
	free(mtxAm);
	free(mtxAn);
	free(mtxAnz);
	free(mtxBm);
	free(mtxBn);
	free(mtxBnz);
	free(mtxCm);
	free(mtxCn);
	free(mtxCnz);

}*/

void dataFiller(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC) {

	//#pragma acc parallel loop
	for (int NZ = 0; NZ < *posC; NZ++) {
		spC->coords[NZ].i = mtxResAm[NZ];
		spC->coords[NZ].j = mtxResAn[NZ];
		spC->data[NZ] = mtxResAnz[NZ];
	}
}

void dataFillerv2(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC) {

	//#pragma acc parallel loop
	for (int NZ = 0; NZ < *posC; NZ++) {
		spC->coords[NZ].i = mtxResAm[NZ];
		spC->coords[NZ].j = mtxResAn[NZ];
		spC->data[NZ] = mtxResAnz[NZ];
	}
}


int qsorter(const void *i, const void *j) {
	int ii;
	int jj;

	ii = ((struct sortCOO *)i)->ij;
	jj = ((struct sortCOO *)j)->ij;

	//printf("test %d %d \n",ii,jj);
	return ii - jj;
}

/* Computes C = A*B.
* C should be allocated by this routine.
*/
void optimised_sparsemm(const COO A, const COO B, COO *C)
{
	struct sortCOO *list1 = malloc(A->NZ * sizeof(struct sortCOO));

	struct sortCOO *list2 = malloc(B->NZ * sizeof(struct sortCOO));

	/*COO newA;
	alloc_sparse(A->m, A->n, A->NZ, &newA);
	COO newB;
	alloc_sparse(B->m, B->n, B->NZ, &newB);*/

		#pragma acc parallel loop
		for (int i = 0; i < A->NZ; i++) {
			list1[i].i = A->coords[i].i;
			list1[i].j = A->coords[i].j;
			list1[i].ij = A->coords[i].i * A->n + A->coords[i].j;
			list1[i].d = A->data[i];
		}
		qsort(list1, A->NZ, sizeof(struct sortCOO), qsorter);

		#pragma acc parallel loop
			for (int i = 0; i < A->NZ; i++) {
				list2[i].i = B->coords[i].j;
				list2[i].j = B->coords[i].i;
				list2[i].ij = B->coords[i].j * B->m + B->coords[i].i;
				list2[i].d = B->data[i];
			}
			qsort(list2, B->NZ, sizeof(struct sortCOO), qsorter);


			//CSR for A
			double *mtxCSRnz = malloc(A->NZ * sizeof(double));//data of NZ
			int *mtxCSRm = malloc((A->m + 1) * sizeof(int));// first zero last cumulative per row
			mtxCSRm[0] = 0;
			int *mtxCSRn = malloc(A->NZ * sizeof(int));//all coloumn indice

			//CSR for B
			double *mtxCSCnz = malloc(B->NZ * sizeof(double));
			int *mtxCSCm = malloc((B->n + 1) * sizeof(int));
			mtxCSRm[0] = 0;
			int *mtxCSCn = malloc(B->NZ * sizeof(int));

			CSRmaker(list1, mtxCSRnz, mtxCSRn, mtxCSRm, &A->NZ);
			CSRmaker(list2, mtxCSCnz, mtxCSCn, mtxCSCm, &B->NZ);

			for (int ii = 0; ii < A->m + 1; ii++)
			{
			printf("%d, ", mtxCSRm[ii]);

		}

		printf("\n CSRm ");

		for (int ii = 0; ii < A->m + 1; ii++)
		{
		printf("%d, ", mtxCSCm[ii]);

	}
		printf("hey\n" );

			/*int counter = 1;
			int curCol = 0;
			//fill CSR matrix
			for (int i = 0; i < *Anz; i++)
			{
				mtxCSRnz[i] = A->data[i];
				mtxCSRn[i] = A->coords[i].j;

				if (A->coords[i].i > curCol) {
					while (A->coords[i].i > curCol)
					{
						mtxCSRm[counter] = i;
						curCol++;
						counter++;
					}
				}
			}*/


	/*structSort(list1, A, newA);
	structSort(list2, B, newB);*/

/*	print_sparse(A);
	printf("next \n");
	print_sparse(B);
	printf("next \n");

	for (int i = 0; i < A->NZ; i++) {
		printf("i %d j %d ij %d d %f \n", list1[i].i,list1[i].j,list1[i].ij, list1[i].d );
	}
printf("\n");
	for (int i = 0; i < B->NZ; i++) {

		printf("i %d j %d ij %d d %f \n", list2[i].i,list2[i].j,list2[i].ij, list2[i].d );
	}*/

  /*
	printf("A CSR init: \n");
	clock_t begin = clock();

	int Am, An, Anz;
	Am = A->m;
	An = A->n;
	Anz = A->NZ;

	double *mtxCSRnz = malloc(Anz * sizeof(double));
	int *mtxCSRm = malloc((Am + 1) * sizeof(int));
	mtxCSRm[0] = 0;
	mtxCSRm[Am] = Anz;
	int *mtxCSRn = malloc(Anz * sizeof(int));

	//fill CSR matrix
	CSRmaker(newA, &Anz, mtxCSRnz, mtxCSRn, mtxCSRm);

	printf("B CSC init: \n");

	int Bm, Bn, Bnz;
	Bm = newB->m;
	Bn = newB->n;
	Bnz = newB->NZ;

	double *mtxCSCnz = malloc(Bnz * sizeof(double));
	int *mtxCSCm = malloc(Bnz * sizeof(int));//seg fault
	int *mtxCSCn = malloc((Bn + 1) * sizeof(int));
	mtxCSCn[0] = 0;

	int counter = 0;
	for (int j = 0; j < Bn; ++j)
	{
		for (int i = 0; i < Bnz; ++i)
		{
			if (j == newB->coords[i].j)
			{
				mtxCSCnz[counter] = newB->data[i];
				mtxCSCm[counter] = newB->coords[i].i;
				++counter;
			}

		}
		mtxCSCn[j + 1] = counter;
	}

//CSC part
	int point = 1;
	int getColSize = 0;
	int pointCSC = 1;
	double *CSC2Vec = calloc(Bm, sizeof(double));

	//CSR part
	int csrCount = 1;
	int getRowSize = 0;
	int pointCSR = 1;

	//results
	double *c = calloc(Am*Bn, sizeof(double));//maybe turn to

	printf("sizes: %d,%d \n", Am, Bn);
	double ress = 0;




	printf("MM starts \n");
	for (int lister = 1; lister < Bnz + 1;lister++)
	{

		getColSize = mtxCSCn[point] - mtxCSCn[point - 1];
		CSC2Vec[mtxCSCm[lister - 1]] = mtxCSCnz[lister - 1];

		if (pointCSC == getColSize)
		{
			getRowSize = mtxCSRm[pointCSR] - mtxCSRm[pointCSR - 1];

			for (int csrInd = 1; csrInd < Anz + 1; csrInd++)
			{
				ress += CSC2Vec[mtxCSRn[csrInd - 1]] * mtxCSRnz[csrInd - 1];


				if (csrCount == getRowSize)
				{
					c[(pointCSR - 1) * Bn + (point - 1)] = ress;//still change the structure.

					pointCSR++;
					csrCount = 1;
					ress = 0;

					if (csrInd != Anz)
					{
						getRowSize = mtxCSRm[pointCSR] - mtxCSRm[pointCSR - 1];
					}

				}
				else {
					++csrCount;
				}
			}

			pointCSR = 1;
			memset(CSC2Vec, 0, Bn * sizeof(double));
			++point;
			pointCSC = 1;
		}
		else {
			++pointCSC;
		}
	}
	printf("MM done \n");

	free(mtxCSCnz);
	free(mtxCSCn);
	free(mtxCSCm);
	free(mtxCSRnz);
	free(mtxCSRn);
	free(mtxCSRm);
	free(CSC2Vec);


	int NZ = 0;
	/*IMPORTANT FOR NZ-
	for (int ii = 0; ii < Am*Bn; ii++)
	{
		//printf("%f,", c[ii]);
		if (fabs(c[ii]) > 1e-15) {
			NZ++;
		}
		if ((ii + 1) % Bn == 0 && ii > 0) {
			//printf("\n");
		}

	}
	//giant matrix booms

	double *mtxResAnz = malloc((NZ) * sizeof(double));
	int *mtxResCol = malloc((NZ) * sizeof(int));
	int *mtxResRow = malloc((NZ) * sizeof(int));//not csr form, just row indice!!

	int tmpNZ = 0;
	int colChange = 0;
	for (int ii = 0; ii < Am*Bn; ii++)
	{

		if (fabs(c[ii]) > 1e-15) {
			mtxResAnz[tmpNZ] = c[ii];
			mtxResCol[tmpNZ] = ii % Bn;
			mtxResRow[tmpNZ] = colChange;
			tmpNZ++;
		}

		if ((ii + 1) % Bn == 0 && ii > 0) {
			//printf("tmp %d,", tmpNZ);
			//printf("\n");
			colChange++;

		}

	}

	alloc_sparse(Am, Bn, NZ, C);

	dataFillerv2(*C, mtxResCol, mtxResAnz, mtxResRow, &NZ);

	free(mtxResRow);
	free(mtxResCol);
	free(mtxResAnz);


	free(c);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("1st Loop took %g seconds\n", time_spent);*/
}




/* Computes O = (A + B + C) (D + E + F).
* O should be allocated by this routine.
*/
void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
	const COO D, const COO E, const COO F,
	COO *O)
{

/*
		struct sortCOO *list1 = malloc(A->NZ * sizeof(struct sortCOO));
		struct sortCOO *list2 = malloc(B->NZ * sizeof(struct sortCOO));
		struct sortCOO *list3 = malloc(C->NZ * sizeof(struct sortCOO));
		struct sortCOO *list4 = malloc(D->NZ * sizeof(struct sortCOO));
		struct sortCOO *list5 = malloc(E->NZ * sizeof(struct sortCOO));
		struct sortCOO *list6 = malloc(F->NZ * sizeof(struct sortCOO));


		COO newA;
		alloc_sparse(A->m, A->n, A->NZ, &newA);
		COO newB;
		alloc_sparse(B->m, B->n, B->NZ, &newB);
		COO newC;
		alloc_sparse(C->m, C->n, C->NZ, &newC);
		COO newD;
		alloc_sparse(D->m, D->n, D->NZ, &newD);
		COO newE;
		alloc_sparse(E->m, E->n, E->NZ, &newE);
		COO newF;
		alloc_sparse(F->m, F->n, F->NZ, &newF);


		structSort(list1, A, newA);
		structSort(list2, B, newB);
		structSort(list3, C, newC);
		structSort(list4, D, newD);
		structSort(list5, E, newE);
		structSort(list6, F, newF);

		COO spC;
		COO spD;
		mtxGiantAdder(newA, newB, newC, &spC);
		mtxGiantAdder(newD, newE, newF, &spD);

		optimised_sparsemm(spC, spD, O);

		printf("MM done!! \n");
		free(newA);
		free(newB);
		free(newC);
		free(newD);
		free(newE);
		free(newF);
		free_sparse(&spC);
		free_sparse(&spD);*/

}
