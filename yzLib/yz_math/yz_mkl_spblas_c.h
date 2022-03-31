/***********************************************************/
/**	\file
	\brief		MKL spblas c interface
	\details	Although mkl spblas provide c interface, but 
				it is very hard to use because it is fortran 
				style. In this file, all the interface are 
				changed to c style interface. with c_ as prefix
				for each function. \n
				Currently, I just keep real value functions 
				with lower case.
	\author		Yizhong Zhang
	\date		9/11/2012
*/
/***********************************************************/
#ifndef __YZ_MKL_SPBLAS_C_H__
#define __YZ_MKL_SPBLAS_C_H__

#include "yzLib/yz_setting.h"

#ifndef YZ_mkl_spblas_h
#	error yz_mkl_spblas_c.h must be included after mkl_spblas.h
#endif

namespace yz{	namespace mkl{

//	========================================
///@{
/**	@name Float, Sparse BLAS Level2 lower case c interface
*/
//	========================================
inline void c_mkl_scsrmv(char transa, MKL_INT m, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float beta, float *y){
	mkl_scsrmv(&transa, &m, &k, &alpha, matdescra, val, indx,  pntrb, pntre, x, &beta, y);
}

inline void c_mkl_scsrsv(char transa, MKL_INT m, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y){
	mkl_scsrsv(&transa, &m, &alpha, matdescra, val, indx,  pntrb, pntre, x, y);
}

inline void c_mkl_scsrgemv(char transa, MKL_INT m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_scsrgemv(&transa, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_scsrgemv(char transa, MKL_INT m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_cspblas_scsrgemv(&transa, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_scsrsymv(char uplo, MKL_INT m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_scsrsymv(&uplo, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_scsrsymv(char uplo, MKL_INT m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_cspblas_scsrsymv(&uplo, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_scsrtrsv(char uplo, char transa, char diag, MKL_INT m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_scsrtrsv(&uplo, &transa, &diag, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_scsrtrsv(char uplo, char transa, char diag, MKL_INT m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_cspblas_scsrtrsv(&uplo, &transa, &diag, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_scscmv(char transa, MKL_INT m, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float beta, float *y){
	mkl_scscmv(&transa, &m, &k, &alpha, matdescra, val, indx,  pntrb, pntre, x, &beta, y);
}

inline void c_mkl_scscsv(char transa, MKL_INT m, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y){
	mkl_scscsv(&transa, &m, &alpha, matdescra, val, indx,  pntrb, pntre, x, y);
}

inline void c_mkl_scoomv(char transa, MKL_INT m, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, float *x, float beta, float *y){
	mkl_scoomv(&transa, &m, &k, &alpha, matdescra, val, rowind,  colind, &nnz, x, &beta, y);
}

inline void c_mkl_scoosv(char transa, MKL_INT m, float alpha, char *matdescra, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, float *x, float *y){
	mkl_scoosv(&transa, &m, &alpha, matdescra, val, rowind,  colind, &nnz, x, y);
}

inline void c_mkl_scoogemv(char transa, MKL_INT m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, float *x,  float *y){
	mkl_scoogemv(&transa, &m, val, rowind,  colind, &nnz, x,  y);
}

inline void c_mkl_cspblas_scoogemv(char transa, MKL_INT m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, float *x,  float *y){
	mkl_cspblas_scoogemv(&transa, &m, val, rowind,  colind, &nnz, x,  y);
}

inline void c_mkl_scoosymv(char uplo, MKL_INT m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, float *x,  float *y){
	mkl_scoosymv(&uplo, &m, val, rowind,  colind, &nnz, x,  y);
}

inline void c_mkl_cspblas_scoosymv(char uplo, MKL_INT m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, float *x,  float *y){
	mkl_cspblas_scoosymv(&uplo, &m, val, rowind,  colind, &nnz, x,  y);
}

inline void c_mkl_scootrsv(char uplo, char transa, char diag, MKL_INT m, float *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT nnz, float *x,  float *y){
	mkl_scootrsv(&uplo, &transa, &diag, &m, val, rowind, colind, &nnz, x,  y);
}

inline void c_mkl_cspblas_scootrsv(char uplo, char transa, char diag, MKL_INT m, float *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT nnz, float *x,  float *y){
	mkl_cspblas_scootrsv(&uplo, &transa, &diag, &m, val, rowind, colind, &nnz, x,  y);
}

inline void c_mkl_sdiamv (char transa, MKL_INT m, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT lval, MKL_INT *idiag,  MKL_INT ndiag, float *x, float beta, float *y){
	mkl_sdiamv (&transa, &m, &k, &alpha, matdescra, val, &lval,	idiag,  &ndiag, x, &beta, y);
}

inline void c_mkl_sdiasv (char transa, MKL_INT m, float alpha, char *matdescra, float  *val, MKL_INT lval, MKL_INT *idiag,  MKL_INT ndiag, float *x, float *y){
	mkl_sdiasv (&transa, &m, &alpha, matdescra, val, &lval, idiag,  &ndiag, x, y);
}

inline void c_mkl_sdiagemv(char transa, MKL_INT m, float *val, MKL_INT lval,  MKL_INT *idiag, MKL_INT ndiag, float *x,  float *y){
	mkl_sdiagemv(&transa, &m, val, &lval,  idiag, &ndiag, x,  y);
}

inline void c_mkl_sdiasymv(char uplo, MKL_INT m, float *val, MKL_INT lval,  MKL_INT *idiag, MKL_INT ndiag, float *x,  float *y){
	mkl_sdiasymv(&uplo, &m, val, &lval,  idiag, &ndiag, x,  y);
}

inline void c_mkl_sdiatrsv(char uplo, char transa, char diag, MKL_INT m, float *val, MKL_INT lval,  MKL_INT  *idiag, MKL_INT ndiag, float *x,  float *y){
	mkl_sdiatrsv(&uplo, &transa, &diag, &m, val, &lval,  idiag, &ndiag, x,  y);
}

inline void c_mkl_sskymv (char transa, MKL_INT m, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT *pntr, float *x, float beta, float *y){
	mkl_sskymv (&transa, &m, &k, &alpha, matdescra, val, pntr, x, &beta, y);
}

inline void c_mkl_sskysv(char transa, MKL_INT m, float alpha, char *matdescra, float  *val, MKL_INT *pntr,  float *x, float *y){
	mkl_sskysv(&transa, &m, &alpha, matdescra, val, pntr,  x, y);
}

inline void c_mkl_sbsrmv (char transa, MKL_INT m, MKL_INT k, MKL_INT lb, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float beta, float *y){
	mkl_sbsrmv (&transa, &m, &k, &lb, &alpha, matdescra, val, indx,  pntrb, pntre, x, &beta, y);
}

inline void c_mkl_sbsrsv(char transa, MKL_INT m, MKL_INT lb, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y){
	mkl_sbsrsv(&transa, &m, &lb, &alpha, matdescra, val, indx,  pntrb, pntre, x, y);
}

inline void c_mkl_sbsrgemv(char transa, MKL_INT m, MKL_INT lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_sbsrgemv(&transa, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_sbsrgemv(char transa, MKL_INT m, MKL_INT lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_cspblas_sbsrgemv(&transa, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_sbsrsymv(char uplo, MKL_INT m, MKL_INT lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_sbsrsymv(&uplo, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_sbsrsymv(char uplo, MKL_INT m, MKL_INT lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_cspblas_sbsrsymv(&uplo, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_sbsrtrsv(char uplo, char transa, char diag, MKL_INT m, MKL_INT lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_sbsrtrsv(&uplo, &transa, &diag, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_sbsrtrsv(char uplo, char transa, char diag, MKL_INT m, MKL_INT lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y){
	mkl_cspblas_sbsrtrsv(&uplo, &transa, &diag, &m, &lb, a, ia,  ja, x,  y);
}

///@}

//	========================================
///@{
/**	@name Float, Sparse BLAS Level3 lower case c interface
*/
//	========================================
inline void c_mkl_scsrmm(char transa, MKL_INT m, MKL_INT n, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT ldb, float beta, float *c, MKL_INT ldc){
	mkl_scsrmm(&transa, &m, &n, &k, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_scsrsm(char transa, MKL_INT m, MKL_INT n, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT ldb,  float *c, MKL_INT ldc){
	mkl_scsrsm(&transa, &m, &n, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb,  c, &ldc);
}

inline void c_mkl_scscmm(char transa, MKL_INT m, MKL_INT n, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT ldb, float beta, float *c, MKL_INT ldc){
	mkl_scscmm(&transa, &m, &n, &k, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_scscsm(char transa, MKL_INT m, MKL_INT n, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT ldb,  float *c, MKL_INT ldc){
	mkl_scscsm(&transa, &m, &n, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb,  c, &ldc);
}

inline void c_mkl_scoomm(char transa, MKL_INT m, MKL_INT n, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, float *b, MKL_INT ldb, float beta, float *c, MKL_INT ldc){
	mkl_scoomm(&transa, &m, &n, &k, &alpha, matdescra, val, rowind,  colind, &nnz, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_scoosm(char transa, MKL_INT m, MKL_INT n, float alpha, char *matdescra, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, float *b, MKL_INT ldb,  float *c, MKL_INT ldc){
	mkl_scoosm(&transa, &m, &n, &alpha, matdescra, val, rowind,  colind, &nnz, b, &ldb,  c, &ldc);
}

inline void c_mkl_sdiamm (char transa, MKL_INT m, MKL_INT n, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT lval, MKL_INT *idiag,  MKL_INT ndiag, float *b, MKL_INT ldb, float beta, float *c, MKL_INT ldc){
	mkl_sdiamm (&transa, &m, &n, &k, &alpha, matdescra, val, &lval, idiag,  &ndiag, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_sdiasm (char transa, MKL_INT m, MKL_INT n, float alpha, char *matdescra, float  *val, MKL_INT lval, MKL_INT *idiag,  MKL_INT ndiag, float *b, MKL_INT ldb, float *c, MKL_INT ldc){
	mkl_sdiasm (&transa, &m, &n, &alpha, matdescra, val, &lval, idiag,  &ndiag, b, &ldb, c, &ldc);
}

inline void c_mkl_sskysm (char transa, MKL_INT m, MKL_INT n, float alpha, char *matdescra, float  *val, MKL_INT *pntr,  float *b, MKL_INT ldb, float *c, MKL_INT ldc){
	mkl_sskysm (&transa, &m, &n, &alpha, matdescra, val, pntr,  b, &ldb, c, &ldc);
}

inline void c_mkl_sskymm (char transa, MKL_INT m, MKL_INT n, MKL_INT k, float alpha, char *matdescra, float  *val, MKL_INT *pntr, float *b, MKL_INT ldb, float beta, float *c, MKL_INT ldc){
	mkl_sskymm (&transa, &m, &n, &k, &alpha, matdescra, val, pntr, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_sbsrmm(char transa, MKL_INT m, MKL_INT n, MKL_INT k, MKL_INT lb, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT ldb, float beta, float *c, MKL_INT ldc){
	mkl_sbsrmm(&transa, &m, &n, &k, &lb, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_sbsrsm(char transa, MKL_INT m, MKL_INT n, MKL_INT lb, float alpha, char *matdescra, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT ldb,  float *c, MKL_INT ldc){
	mkl_sbsrsm(&transa, &m, &n, &lb, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb,  c, &ldc);
}

///@}

//	========================================
///@{
/**	@name Double, Sparse BLAS Level2 lower case c interface
*/
//	========================================
inline void c_mkl_dcsrmv(char transa, MKL_INT m, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double beta, double *y){
	mkl_dcsrmv(&transa, &m, &k, &alpha, matdescra, val, indx,  pntrb, pntre, x, &beta, y);
}

inline void c_mkl_dcsrsv(char transa, MKL_INT m, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y){
	mkl_dcsrsv(&transa, &m, &alpha, matdescra, val, indx,  pntrb, pntre, x, y);
}

inline void c_mkl_dcsrgemv(char transa, MKL_INT m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_dcsrgemv(&transa, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_dcsrgemv(char transa, MKL_INT m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_cspblas_dcsrgemv(&transa, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_dcsrsymv(char uplo, MKL_INT m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_dcsrsymv(&uplo, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_dcsrsymv(char uplo, MKL_INT m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_cspblas_dcsrsymv(&uplo, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_dcsrtrsv(char uplo, char transa, char diag, MKL_INT m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_dcsrtrsv(&uplo, &transa, &diag, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_dcsrtrsv(char uplo, char transa, char diag, MKL_INT m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_cspblas_dcsrtrsv(&uplo, &transa, &diag, &m, a, ia,  ja, x,  y);
}

inline void c_mkl_dcscmv(char transa, MKL_INT m, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double beta, double *y){
	mkl_dcscmv(&transa, &m, &k, &alpha, matdescra, val, indx,  pntrb, pntre, x, &beta, y);
}

inline void c_mkl_dcscsv(char transa, MKL_INT m, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y){
	mkl_dcscsv(&transa, &m, &alpha, matdescra, val, indx,  pntrb, pntre, x, y);
}

inline void c_mkl_dcoomv(char transa, MKL_INT m, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, double *x, double beta, double *y){
	mkl_dcoomv(&transa, &m, &k, &alpha, matdescra, val, rowind,  colind, &nnz, x, &beta, y);
}

inline void c_mkl_dcoosv(char transa, MKL_INT m, double alpha, char *matdescra, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, double *x, double *y){
	mkl_dcoosv(&transa, &m, &alpha, matdescra, val, rowind,  colind, &nnz, x, y);
}

inline void c_mkl_dcoogemv(char transa, MKL_INT m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, double *x,  double *y){
	mkl_dcoogemv(&transa, &m, val, rowind,  colind, &nnz, x,  y);
}

inline void c_mkl_cspblas_dcoogemv(char transa, MKL_INT m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, double *x,  double *y){
	mkl_cspblas_dcoogemv(&transa, &m, val, rowind,  colind, &nnz, x,  y);
}

inline void c_mkl_dcoosymv(char uplo, MKL_INT m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, double *x,  double *y){
	mkl_dcoosymv(&uplo, &m, val, rowind,  colind, &nnz, x,  y);
}

inline void c_mkl_cspblas_dcoosymv(char uplo, MKL_INT m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, double *x,  double *y){
	mkl_cspblas_dcoosymv(&uplo, &m, val, rowind,  colind, &nnz, x,  y);
}

inline void c_mkl_dcootrsv(char uplo, char transa, char diag, MKL_INT m, double *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT nnz, double *x,  double *y){
	mkl_dcootrsv(&uplo, &transa, &diag, &m, val, rowind, colind, &nnz, x,  y);
}

inline void c_mkl_cspblas_dcootrsv(char uplo, char transa, char diag, MKL_INT m, double *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT nnz, double *x,  double *y){
	mkl_cspblas_dcootrsv(&uplo, &transa, &diag, &m, val, rowind, colind, &nnz, x,  y);
}

inline void c_mkl_ddiamv (char transa, MKL_INT m, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT lval, MKL_INT *idiag,  MKL_INT ndiag, double *x, double beta, double *y){
	mkl_ddiamv (&transa, &m, &k, &alpha, matdescra, val, &lval,	idiag,  &ndiag, x, &beta, y);
}

inline void c_mkl_ddiasv (char transa, MKL_INT m, double alpha, char *matdescra, double  *val, MKL_INT lval, MKL_INT *idiag,  MKL_INT ndiag, double *x, double *y){
	mkl_ddiasv (&transa, &m, &alpha, matdescra, val, &lval, idiag,  &ndiag, x, y);
}

inline void c_mkl_ddiagemv(char transa, MKL_INT m, double *val, MKL_INT lval,  MKL_INT *idiag, MKL_INT ndiag, double *x,  double *y){
	mkl_ddiagemv(&transa, &m, val, &lval,  idiag, &ndiag, x,  y);
}

inline void c_mkl_ddiasymv(char uplo, MKL_INT m, double *val, MKL_INT lval,  MKL_INT *idiag, MKL_INT ndiag, double *x,  double *y){
	mkl_ddiasymv(&uplo, &m, val, &lval,  idiag, &ndiag, x,  y);
}

inline void c_mkl_ddiatrsv(char uplo, char transa, char diag, MKL_INT m, double *val, MKL_INT lval,  MKL_INT  *idiag, MKL_INT ndiag, double *x,  double *y){
	mkl_ddiatrsv(&uplo, &transa, &diag, &m, val, &lval,  idiag, &ndiag, x,  y);
}

inline void c_mkl_dskymv (char transa, MKL_INT m, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT *pntr, double *x, double beta, double *y){
	mkl_dskymv (&transa, &m, &k, &alpha, matdescra, val, pntr, x, &beta, y);
}

inline void c_mkl_dskysv(char transa, MKL_INT m, double alpha, char *matdescra, double  *val, MKL_INT *pntr,  double *x, double *y){
	mkl_dskysv(&transa, &m, &alpha, matdescra, val, pntr,  x, y);
}

inline void c_mkl_dbsrmv (char transa, MKL_INT m, MKL_INT k, MKL_INT lb, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double beta, double *y){
	mkl_dbsrmv (&transa, &m, &k, &lb, &alpha, matdescra, val, indx,  pntrb, pntre, x, &beta, y);
}

inline void c_mkl_dbsrsv(char transa, MKL_INT m, MKL_INT lb, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y){
	mkl_dbsrsv(&transa, &m, &lb, &alpha, matdescra, val, indx,  pntrb, pntre, x, y);
}

inline void c_mkl_dbsrgemv(char transa, MKL_INT m, MKL_INT lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_dbsrgemv(&transa, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_dbsrgemv(char transa, MKL_INT m, MKL_INT lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_cspblas_dbsrgemv(&transa, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_dbsrsymv(char uplo, MKL_INT m, MKL_INT lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_dbsrsymv(&uplo, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_dbsrsymv(char uplo, MKL_INT m, MKL_INT lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_cspblas_dbsrsymv(&uplo, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_dbsrtrsv(char uplo, char transa, char diag, MKL_INT m, MKL_INT lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_dbsrtrsv(&uplo, &transa, &diag, &m, &lb, a, ia,  ja, x,  y);
}

inline void c_mkl_cspblas_dbsrtrsv(char uplo, char transa, char diag, MKL_INT m, MKL_INT lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y){
	mkl_cspblas_dbsrtrsv(&uplo, &transa, &diag, &m, &lb, a, ia,  ja, x,  y);
}

///@}

//	========================================
///@{
/**	@name Float, Sparse BLAS Level3 lower case c interface
*/
//	========================================
inline void c_mkl_dcsrmm(char transa, MKL_INT m, MKL_INT n, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT ldb, double beta, double *c, MKL_INT ldc){
	mkl_dcsrmm(&transa, &m, &n, &k, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_dcsrsm(char transa, MKL_INT m, MKL_INT n, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT ldb,  double *c, MKL_INT ldc){
	mkl_dcsrsm(&transa, &m, &n, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb,  c, &ldc);
}

inline void c_mkl_dcscmm(char transa, MKL_INT m, MKL_INT n, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT ldb, double beta, double *c, MKL_INT ldc){
	mkl_dcscmm(&transa, &m, &n, &k, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_dcscsm(char transa, MKL_INT m, MKL_INT n, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT ldb,  double *c, MKL_INT ldc){
	mkl_dcscsm(&transa, &m, &n, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb,  c, &ldc);
}

inline void c_mkl_dcoomm(char transa, MKL_INT m, MKL_INT n, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, double *b, MKL_INT ldb, double beta, double *c, MKL_INT ldc){
	mkl_dcoomm(&transa, &m, &n, &k, &alpha, matdescra, val, rowind,  colind, &nnz, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_dcoosm(char transa, MKL_INT m, MKL_INT n, double alpha, char *matdescra, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT nnz, double *b, MKL_INT ldb,  double *c, MKL_INT ldc){
	mkl_dcoosm(&transa, &m, &n, &alpha, matdescra, val, rowind,  colind, &nnz, b, &ldb,  c, &ldc);
}

inline void c_mkl_ddiamm (char transa, MKL_INT m, MKL_INT n, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT lval, MKL_INT *idiag,  MKL_INT ndiag, double *b, MKL_INT ldb, double beta, double *c, MKL_INT ldc){
	mkl_ddiamm (&transa, &m, &n, &k, &alpha, matdescra, val, &lval, idiag,  &ndiag, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_ddiasm (char transa, MKL_INT m, MKL_INT n, double alpha, char *matdescra, double  *val, MKL_INT lval, MKL_INT *idiag,  MKL_INT ndiag, double *b, MKL_INT ldb, double *c, MKL_INT ldc){
	mkl_ddiasm (&transa, &m, &n, &alpha, matdescra, val, &lval, idiag,  &ndiag, b, &ldb, c, &ldc);
}

inline void c_mkl_dskysm (char transa, MKL_INT m, MKL_INT n, double alpha, char *matdescra, double  *val, MKL_INT *pntr,  double *b, MKL_INT ldb, double *c, MKL_INT ldc){
	mkl_dskysm (&transa, &m, &n, &alpha, matdescra, val, pntr,  b, &ldb, c, &ldc);
}

inline void c_mkl_dskymm (char transa, MKL_INT m, MKL_INT n, MKL_INT k, double alpha, char *matdescra, double  *val, MKL_INT *pntr, double *b, MKL_INT ldb, double beta, double *c, MKL_INT ldc){
	mkl_dskymm (&transa, &m, &n, &k, &alpha, matdescra, val, pntr, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_dbsrmm(char transa, MKL_INT m, MKL_INT n, MKL_INT k, MKL_INT lb, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT ldb, double beta, double *c, MKL_INT ldc){
	mkl_dbsrmm(&transa, &m, &n, &k, &lb, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb, &beta, c, &ldc);
}

inline void c_mkl_dbsrsm(char transa, MKL_INT m, MKL_INT n, MKL_INT lb, double alpha, char *matdescra, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT ldb,  double *c, MKL_INT ldc){
	mkl_dbsrsm(&transa, &m, &n, &lb, &alpha, matdescra, val, indx,  pntrb, pntre, b, &ldb,  c, &ldc);
}

///@}

//	========================================
///@{
/**	@name Converters lower case
*/
//	========================================
inline void c_mkl_dcsrbsr(MKL_INT * job,MKL_INT m,MKL_INT mblk,MKL_INT ldAbsr,double *Acsr,MKL_INT * AJ,MKL_INT * AI,double *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info){
	mkl_dcsrbsr(job,&m,&mblk,&ldAbsr,Acsr,AJ,AI,Absr, AJB, AIB, info);
}

inline void c_mkl_dcsrcoo(MKL_INT * job,MKL_INT n,double *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT nnz,double *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info){
	mkl_dcsrcoo(job,&n,Acsr,AJR,AIR,&nnz,Acoo, ir, jc, info);
}

inline void c_mkl_ddnscsr(MKL_INT *job,MKL_INT m,MKL_INT n,double *Adns,MKL_INT lda,double *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info){
	mkl_ddnscsr(job,&m,&n,Adns,&lda,Acsr,AJ,AI,info);
}

inline void c_mkl_dcsrcsc(MKL_INT * job,MKL_INT n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0,double *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info){
	mkl_dcsrcsc(job,&n,Acsr,AJ0,AI0,Acsc,AJ1,AI1,info);
}

inline void c_mkl_dcsrdia(MKL_INT * job,MKL_INT n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0,double *Adia,MKL_INT ndiag,MKL_INT * distance,MKL_INT idiag,double *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info){
	mkl_dcsrdia(job,&n,Acsr,AJ0,AI0,Adia,&ndiag,distance,&idiag,Acsr_rem,AJ0_rem,AI0_rem,info);
}

inline void c_mkl_dcsrsky(MKL_INT * job,MKL_INT n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0, double *Asky,MKL_INT * pointers,MKL_INT * info){
	mkl_dcsrsky(job,&n,Acsr,AJ0,AI0, Asky,pointers,info);
}

inline void c_mkl_scsrbsr(MKL_INT * job,MKL_INT m,MKL_INT mblk,MKL_INT ldAbsr,float *Acsr,MKL_INT * AJ,MKL_INT * AI,float *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info){
	mkl_scsrbsr(job,&m,&mblk,&ldAbsr,Acsr,AJ,AI,Absr, AJB, AIB, info);
}

inline void c_mkl_scsrcoo(MKL_INT * job,MKL_INT n,float *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT nnz,float *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info){
	mkl_scsrcoo(job,&n,Acsr,AJR,AIR,&nnz,Acoo, ir, jc, info);
}

inline void c_mkl_sdnscsr(MKL_INT *job,MKL_INT m,MKL_INT n,float *Adns,MKL_INT lda,float *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info){
	mkl_sdnscsr(job,&m,&n,Adns,&lda,Acsr,AJ,AI,info);
}

inline void c_mkl_scsrcsc(MKL_INT * job,MKL_INT n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0,float *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info){
	mkl_scsrcsc(job,&n,Acsr,AJ0,AI0,Acsc,AJ1,AI1,info);
}

inline void c_mkl_scsrdia(MKL_INT * job,MKL_INT n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0,float *Adia,MKL_INT ndiag,MKL_INT * distance,MKL_INT idiag,float *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info){
	mkl_scsrdia(job,&n,Acsr,AJ0,AI0,Adia,&ndiag,distance,&idiag,Acsr_rem,AJ0_rem,AI0_rem,info);
}

inline void c_mkl_scsrsky(MKL_INT * job,MKL_INT n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0, float *Asky,MKL_INT * pointers,MKL_INT * info){
	mkl_scsrsky(job,&n,Acsr,AJ0,AI0, Asky,pointers,info);
}



///@}

//	========================================
///@{
/**	@name Sparse BLAS Level2 (CSR-CSR) lower case
*/
//	========================================
inline void c_mkl_dcsrmultcsr(char transa, MKL_INT job, MKL_INT sort, MKL_INT n, MKL_INT k, MKL_INT m, double *a, MKL_INT *ja, MKL_INT *ia, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT nnzmax, MKL_INT *ierr){
	mkl_dcsrmultcsr(&transa, &job, &sort, &n, &k, &m, a, ja, ia, b, jb, ib, c, jc, ic, &nnzmax, ierr);
}

inline void c_mkl_dcsrmultd(char transa,  MKL_INT n, MKL_INT k, MKL_INT m, double *a, MKL_INT *ja, MKL_INT *ia, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT ldc){
	mkl_dcsrmultd(&transa,  &n, &k, &m, a, ja, ia, b, jb, ib, c, &ldc);
}

inline void c_mkl_dcsradd(char transa, MKL_INT job, MKL_INT sort, MKL_INT n, MKL_INT k, double *a, MKL_INT *ja, MKL_INT *ia, double beta, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT nnzmax, MKL_INT *ierr){
	mkl_dcsradd(&transa, &job, &sort, &n, &k, a, ja, ia, &beta, b, jb, ib, c, jc, ic, &nnzmax, ierr);
}

inline void c_mkl_scsrmultcsr(char transa, MKL_INT job, MKL_INT sort, MKL_INT n, MKL_INT k, MKL_INT m, float *a, MKL_INT *ja, MKL_INT *ia, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT *jc, MKL_INT *ic, MKL_INT nnzmax, MKL_INT *ierr){
	mkl_scsrmultcsr(&transa, &job, &sort, &n, &k, &m, a, ja, ia, b, jb, ib, c, jc, ic, &nnzmax, ierr);
}

inline void c_mkl_scsrmultd(char transa,  MKL_INT n, MKL_INT k, MKL_INT m, float *a, MKL_INT *ja, MKL_INT *ia, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT ldc){
	mkl_scsrmultd(&transa,  &n, &k, &m, a, ja, ia, b, jb, ib, c, &ldc);
}

inline void c_mkl_scsradd(char transa, MKL_INT job, MKL_INT sort, MKL_INT n, MKL_INT k, float *a, MKL_INT *ja, MKL_INT *ia, float beta, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT *jc, MKL_INT *ic, MKL_INT nnzmax, MKL_INT *ierr){
	mkl_scsradd(&transa, &job, &sort, &n, &k, a, ja, ia, &beta, b, jb, ib, c, jc, ic, &nnzmax, ierr);
}

///@}

}}	//	end namespace yz::mkl

#endif	//	__YZ_MKL_SPBLAS_C_H__