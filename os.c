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
	void CSRmakerA(struct sortCOO *list1, double *mtxCSRnz, int *mtxCSRn, int *mtxCSRm, COO A);
	void CSRmakerB(struct sortCOO *list1, double *mtxCSRnz, int *mtxCSRn, int *mtxCSRm, COO A);
	void matrixAdder(int *curr, double *C2Vec, int *mtxAn, double *mtxAnz, int *mtxAm);
	void mtxGiantAdder(struct sortCOO *list1,struct sortCOO *list2,struct sortCOO *list3,COO A, COO B, COO C, COO *spC);
	//void dataFiller(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC);
	// void dataFillerv2(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC);
	void dataFiller(struct sortCOO *listRes, int *finalNZ, COO C );
	int qsorter(const void *i, const void *j);
	struct sortCOO * sortedNormal(COO A);
	struct sortCOO * sortedTranspose(COO B);

	void CSRmakerA(struct sortCOO *list1, double *mtxCSRnz, int *mtxCSRn, int *mtxCSRm, COO A) {

		//maybe dont do cumulative but normal 1,3,1 and eradicate the need of back-dependency
		////#pragma ivdep
		for (int i = 0; i < A->NZ; i++) {
			mtxCSRnz[i] = list1[i].d;
			mtxCSRn[i] = list1[i].j;
			//printf("mtxCSRm[%d]=%d \n",list1[i].i+1,i+1);
			mtxCSRm[list1[i].i+1]=i+1;
			//printf("mtxCSRm[%d]=> %d\n",list1[i].i+1,mtxCSRm[list1[i].i+1]);
			//or get for at a time and then the remained. make it in a separate branch in github.

			//++mtxCSRm[list1[i].i]; //and make the rest calloc with NO +1
		}

		/*printf("A:\n");
		for (int i = 0; i < A->m +1; i++) {
			printf("%d,",mtxCSRm[i]);
		}
		printf("new \n");*/

		int sum=0;
		for (int fixer = 1; fixer < A->m +1; fixer++) {

			if (!mtxCSRm[fixer]) {
				mtxCSRm[fixer]=sum;
			}else{
				sum=mtxCSRm[fixer];
			}
		}
		/*for (int i = 0; i < A->m +1; i++) {
			printf("%d,",mtxCSRm[i]);
		}
		printf("\n");*/
	}

	void CSRmakerB(struct sortCOO *list1, double *mtxCSRnz, int *mtxCSRn, int *mtxCSRm, COO A) {

		//maybe dont do cumulative but normal 1,3,1 and eradicate the need of back-dependency

		////#pragma ivdep
		for (int i = 0; i < A->NZ; i++) {
			mtxCSRnz[i] = list1[i].d;
			mtxCSRn[i] = list1[i].j;
			mtxCSRm[list1[i].i+1]=i+1;
			//or get for at a time and then the remained. make it in a separate branch in github.

			//++mtxCSRm[list1[i].i]; //and make the rest calloc with NO +1
		}

		/*printf("B:\n");
		for (int i = 0; i < A->n +1; i++) {
			printf("%d,",mtxCSRm[i]);
		}
		printf("new \n");*/

		int sum=0;
		for (int fixer = 1; fixer < A->n +1; fixer++) {

			if (!mtxCSRm[fixer]) {
				mtxCSRm[fixer]=sum;
			}else{
				sum=mtxCSRm[fixer];//
			}
		}
		/*for (int i = 0; i < A->n +1; i++) {
			printf("%d,",mtxCSRm[i]);
		}
		printf("\n");*/
	}

	void matrixAdder(int *curr, double *C2Vec, int *mtxAn, double *mtxAnz, int *mtxAm) {

		int getARowSize = mtxAm[*curr] - mtxAm[*curr - 1];
		int curPos = mtxAm[*curr - 1];
		int tmpPos = 0;

		#pragma acc parallel loop//??????
		for (int row = 0; row < getARowSize; ++row)
		{
			tmpPos = curPos + row;
			C2Vec[mtxAn[tmpPos]] += mtxAnz[tmpPos];
		}

	}

	void mtxGiantAdder(struct sortCOO *list1,struct sortCOO *list2,struct sortCOO *list3,COO A, COO B, COO C, COO *spC) {

	int Anz=A->NZ;
	int Bnz=B->NZ;
	int Cnz=C->NZ;


	double *mtxAnz = malloc(Anz * sizeof(double));//has nz vals
	int *mtxAm = calloc((A->m + 1) , sizeof(int));//has cumulative els
	//mtxAm[0] = 0;
	// mtxAm[A->m] = Anz;
	int *mtxAn = malloc(Anz * sizeof(int));//col index
	//
	double *mtxBnz = malloc(Bnz * sizeof(double));
	int *mtxBm = calloc((B->m + 1) , sizeof(int));
	//mtxBm[0] = 0;
	// mtxBm[B->m] = Bnz;
	int *mtxBn = malloc(Bnz * sizeof(int));
	//
	double *mtxCnz = malloc(Cnz * sizeof(double));
	int *mtxCm = calloc((C->m + 1) , sizeof(int));
	// mtxCm[C->m] = Cnz;
	int *mtxCn = malloc(Cnz * sizeof(int));

	struct sortCOO *listRes=malloc((Cnz + Bnz + Anz) * sizeof(struct sortCOO));

	//fill CSR matrix

		CSRmakerA(list1, mtxAnz, mtxAn, mtxAm, A);
		CSRmakerA(list2, mtxBnz, mtxBn, mtxBm, B);
		CSRmakerA(list3, mtxCnz, mtxCn, mtxCm, C);

	printf("M+M starts \n");
	printf("%d\n",C->n);
	double *C2Vec = calloc(C->n, sizeof(double));//vectorized format  /////////////problemmmmmmmmmmmmmmmmmmmmm

	int posC = 0;

	////#pragma acc kernels loop
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
				listRes[posC].i=cur;//row indice
				listRes[posC].j=adder;//col indice
				listRes[posC].d=C2Vec[adder];//value of cur vector.
				++posC;
			}
		}
		memset(C2Vec, 0, C->n * sizeof(double));
	}

	alloc_sparse(C->m, C->n, posC, spC);//??

	dataFiller(listRes, &posC,*spC);
	//dataFiller(*spC, mtxResAn, mtxResAnz, mtxResAm, &posC);

	//print_sparse(*spC);

	free(C2Vec);
	free(listRes);
	free(mtxAm);
	free(mtxAn);
	free(mtxAnz);
	free(mtxBm);
	free(mtxBn);
	free(mtxBnz);
	free(mtxCm);
	free(mtxCn);
	free(mtxCnz);


}
/*void dataFillerv2(COO spC, int *mtxResAn, double *mtxResAnz, int *mtxResAm, int *posC) {

	////#pragma acc parallel loop
	for (int NZ = 0; NZ < *posC; NZ++) {
		spC->coords[NZ].i = mtxResAm[NZ];
		spC->coords[NZ].j = mtxResAn[NZ];
		spC->data[NZ] = mtxResAnz[NZ];
	}
}*/

void dataFiller(struct sortCOO *listRes, int *finalNZ, COO C ){


	#pragma acc parallel loop
	//#pragma ivdep
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
	////#pragma ivdep //check if there is overlaps
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
	////#pragma ivdep //check if there is overlaps
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
	int *mtxCSCm = calloc((B->n + 1) , sizeof(int));//wo +1 calloc to get rid of get col size ///////to n
	mtxCSRm[0] = 0;
	int *mtxCSCn = malloc(B->NZ * sizeof(int));

	printf("create CSR\n");
	CSRmakerA(list1, mtxCSRnz, mtxCSRn, mtxCSRm, A);
	CSRmakerB(list2, mtxCSCnz, mtxCSCn, mtxCSCm, B);
	printf("CSR DONE\n");

	printf("MM start\n");
	double *vector = calloc(A->n , sizeof(double));
	double res=0;
	int finalNZ=0;
	struct sortCOO *listRes=malloc(0* sizeof(struct sortCOO));//(struct sortCOO *)malloc(0* sizeof(struct sortCOO));

// #pragma acc parallel{}
	for (int i = 0; i < A->m; ++i) {//for all cols or NZ.
		for (int rowMaker = mtxCSRm[i]; rowMaker < mtxCSRm[i+1]; ++rowMaker) {
			vector[mtxCSRn[rowMaker]] = mtxCSRnz[rowMaker];
		}
//vector print
		/*for (int ij	 = 0; ij < A->n; ij++) {
				printf("%f,",	vector[ij] );
			}
			printf("\n");*/

	for (int j = 0; j < B->n; ++j) {//for all cols /////////7to n

		//printf()
			for (int colMaker = mtxCSCm[j]; colMaker < mtxCSCm[j+1]; ++colMaker) {
				res+=vector[mtxCSCn[colMaker]]*mtxCSCnz[colMaker];
			}
				if (res) {
					listRes=realloc(listRes, (finalNZ+1) * sizeof(struct sortCOO));
					listRes[finalNZ].i=i;
					listRes[finalNZ].j=j;
					//listRes[finalNZ].ij=i* A->n +j;
					listRes[finalNZ].d=res;
					++finalNZ;
					res=0;
				}
		}
	memset(vector, 0, A->n * sizeof(double));
	}

	printf("MM DONE %d %d %d \n", A->m, B->n,finalNZ);
	alloc_sparse(A->m, B->n, finalNZ, C);
	dataFiller(listRes, &finalNZ, *C);


	//print_sparse(*C);


	free(list1);
	free(list2);
	free(listRes);
	//printf("lr EOF\n");

	free(mtxCSRnz);

	//printf("rnz EOF\n");

	free(mtxCSRn);//
	//printf("r n EOF\n");

	free(mtxCSRm);//
	//printf("r m EOF\n");

	free(mtxCSCnz);

	free(mtxCSCn);
	//printf("c n EOF\n");

	free(mtxCSCm);//
	//printf("c m EOF\n");

	free(vector);
	printf("EOF\n");

}

/* Computes O = (A + B + C) (D + E + F).
* O should be allocated by this routine.
*/
void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
	const COO D, const COO E, const COO F,
	COO *O)
	{


		struct sortCOO *list1 = sortedNormal(A);
		struct sortCOO *list2 = sortedNormal(B);
		struct sortCOO *list3 = sortedNormal(C);
		struct sortCOO *list4 = sortedNormal(D);
		struct sortCOO *list5 = sortedNormal(E);
		struct sortCOO *list6 = sortedNormal(F);

		COO spC;
		COO spD;
		mtxGiantAdder(list1,list2,list3,A, B, C, &spC);
		mtxGiantAdder(list4,list5,list6,D, E, F, &spD);

		optimised_sparsemm(spC, spD, O);

		printf("M+M done!! \n");
		free_sparse(&spC);
		free_sparse(&spD);

	}
