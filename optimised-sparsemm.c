#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

//http://shodhganga.inflibnet.ac.in/bitstream/10603/151755/11/11_chapter%203.pdf 44

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
	const COO, const COO, const COO,
	COO *);
	void CSRmaker(struct sortCOO *list1, double *mtxCSRnz, int *mtxCSRn, int *mtxCSRm, COO A);
	
	void matrixAdder(int *curr, double *C2Vec, int *mtxAn, double *mtxAnz, int *mtxAm);
	
	void mtxGiantAdder(COO A, COO B, COO C, COO *spC);
	
	void dataFiller(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC);
	void dataFillerv2(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC);
	void dataFillerV3(struct sortCOO *listRes, int *finalNZ, COO C );
	
	int qsorter(const void *i, const void *j);
	
	struct sortCOO * sortedNormal(COO A);
	struct sortCOO * sortedTranspose(COO B);
	
	
	void CSRmaker(struct sortCOO *list1, double *mtxCSRnz, int *mtxCSRn, int *mtxCSRm, COO A) {
		
		//maybe dont do cumulative but normal 1,3,1 and eradicate the need of back-dependency
		//#pragma ivdep
		for (int i = 0; i < A->NZ; i++) {
			mtxCSRnz[i] = list1[i].d;
			mtxCSRn[i] = list1[i].j;
			mtxCSRm[list1[i].i+1]=i+1;
			//or get for at a time and then the remained. make it in a separate branch in github.
			
			//++mtxCSRm[list1[i].i]; //and make the rest calloc with NO +1
		}
		
		int sum=0;
		for (int fixer = 1; fixer < A->n +1; fixer++) {
			
			if (!mtxCSRm[fixer]) {
				mtxCSRm[fixer]=sum;
			}else{
				sum=mtxCSRm[fixer];
			}
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

void dataFillerV3(struct sortCOO *listRes, int *finalNZ, COO C ){
	
	
	#pragma acc parallel loop
	#pragma ivdep
	for (int i = 0; i < *finalNZ; i++) {
		C->coords[i].i = listRes[i].i;
		C->coords[i].j = listRes[i].j;
		C->data[i] = listRes[i].d;
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

struct sortCOO * sortedNormal(COO A){
	
	struct sortCOO *list1 = malloc(A->NZ * sizeof(struct sortCOO));
	
	#pragma acc parallel loop
	//#pragma ivdep //check if there is overlaps
	for (int i = 0; i < A->NZ; i++) {
		list1[i].i = A->coords[i].i;
		list1[i].j = A->coords[i].j;
		list1[i].ij = A->coords[i].i * A->n + A->coords[i].j;
		list1[i].d = A->data[i];
	}
	
	//or get for at a time and then the remained. make it in a separate branch in github.
	qsort(list1, A->NZ, sizeof(struct sortCOO), qsorter);
	
	return list1;
}

struct sortCOO * sortedTranspose(COO B){
	
	struct sortCOO *list2 = malloc(B->NZ * sizeof(struct sortCOO));
	
	
	#pragma acc parallel loop
	//#pragma ivdep //check if there is overlaps
	for (int i = 0; i < B->NZ; i++) {
		list2[i].i = B->coords[i].j;
		list2[i].j = B->coords[i].i;
		list2[i].ij = B->coords[i].j * B->m + B->coords[i].i;
		list2[i].d = B->data[i];
	}
	//or get for at a time and then the remained. make it in a separate branch in github.
	qsort(list2, B->NZ, sizeof(struct sortCOO), qsorter);
	
	return list2;
}

/* Computes C = A*B.
* C should be allocated by this routine.
*/
void optimised_sparsemm(const COO A, const COO B, COO *C)
{
	
	printf("create lists\n");
	struct sortCOO *list1 = sortedNormal(A);//malloc(A->NZ * sizeof(struct sortCOO));
	
	struct sortCOO *list2 = sortedTranspose(B);
	
	
	
	
	//CSR for A
	double *mtxCSRnz = malloc(A->NZ * sizeof(double));//data of NZ
	int *mtxCSRm = calloc((A->m + 1)  , sizeof(int));// first zero last cumulative per row //wo +1 calloc to get rid of get col size
	mtxCSRm[0] = 0;
	int *mtxCSRn = malloc(A->NZ * sizeof(int));//all coloumn indice
	
	//CSR for B
	double *mtxCSCnz = malloc(B->NZ * sizeof(double));
	int *mtxCSCm = calloc((B->n + 1) , sizeof(int));//wo +1 calloc to get rid of get col size
	mtxCSRm[0] = 0;
	int *mtxCSCn = malloc(B->NZ * sizeof(int));
	
	printf("create CSR\n");
	CSRmaker(list1, mtxCSRnz, mtxCSRn, mtxCSRm, A);
	
	CSRmaker(list2, mtxCSCnz, mtxCSCn, mtxCSCm, B);
	printf("CSR DONE\n");
	
	/*
	for (int ii = 0; ii < A->m + 1; ii++)
	{
	printf("%d, ", mtxCSRm[ii]);
}
printf("hey1\n");

printf("\n CSRm ");

for (int ii = 0; ii < A->m + 1; ii++)
{
printf("%d, ", mtxCSCm[ii]);

}
printf("hey2\n" );*/


	printf("MM start\n");



	double *vector = calloc(A->n , sizeof(double));


	double res=0;


	int finalNZ=0;


	struct sortCOO *listRes=NULL;malloc(0* sizeof(struct sortCOO));//(struct sortCOO *)malloc(0* sizeof(struct sortCOO));
	//struct sortCOO *tmp=NULL;

	for (int i = 0; i < A->m; i++) {//for all cols or NZ.
		for (int rowMaker = mtxCSRm[i]; rowMaker < mtxCSRm[i+1]; rowMaker++) {
			//printf("indice[%d]= value: %f \n",mtxCSRn[rowMaker] , mtxCSRnz[rowMaker]);
			vector[mtxCSRn[rowMaker]] = mtxCSRnz[rowMaker];
			//printf("res: %f \n", Vector[mtxCSRn[rowMaker]] );
		}
		
		/*for (int ij	 = 0; ij < A->n; ij++) {
		printf("%f,",	Vector[ij] );
	}
	printf("\n");*/

	//#pragma acc parallel loop . ress can be overlapped. if you really want, use copy.
	for (int j = 0; j < B->n; j++) {//for all cols
		
		//#pragma acc parallel loop //do vectorization %%
		for (int colMaker = mtxCSCm[j]; colMaker < mtxCSCm[j+1]; colMaker++) {
			res+=vector[mtxCSCn[colMaker]]*mtxCSCnz[colMaker];
		}
		
		if (res) {//size_t n = sizeof(a)/sizeof(a[0]); if 0, then all 0.
			//printf("hey %d\n",res);
			//printf("ls: %d \n",&listRes);
			
			listRes=realloc(listRes, (finalNZ+1) * sizeof(struct sortCOO));
			listRes[finalNZ].i=i;
			listRes[finalNZ].j=j;
			listRes[finalNZ].ij=i* A->n +j;
			listRes[finalNZ].d=res;
			//printf("i %d, j %d, ij= %d, d %f \n",listRes[finalNZ].i,listRes[finalNZ].j,listRes[finalNZ].ij,listRes[finalNZ].d);
			++finalNZ;
			res=0;
		}
	}
	memset(vector, 0, A->n * sizeof(double));
	}

	printf("MM DONE %d %d %d \n", A->m, B->n,finalNZ);
	alloc_sparse(A->m, B->n, finalNZ, C);
	dataFillerV3(listRes, &finalNZ, *C);
	printf("EOF\n");

	free(listRes);
	free(mtxCSRnz);
	free(mtxCSRn);
	free(mtxCSRm);
	free(mtxCSCnz);
	free(mtxCSCn);
	free(mtxCSCm);
	free(vector);

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
