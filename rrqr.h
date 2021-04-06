/* Starting from version 7.8, MATLAB BLAS expects ptrdiff_t arguments for integers */
#if MATLAB_VERSION >= 0x0708
#include <stddef.h>
#include <stdlib.h>
#endif 
#include <string.h>
        
/* Define MX_HAS_INTERLEAVED_COMPLEX for version <9.4 */
#ifndef MX_HAS_INTERLEAVED_COMPLEX
#define MX_HAS_INTERLEAVED_COMPLEX 0
#endif

/* Starting from version 7.6, MATLAB BLAS is seperated */
#if MATLAB_VERSION >= 0x0705
#include <blas.h>
#else
#define dcabs1 FORTRAN_WRAPPER(dcabs1)
extern doublereal dcabs1(
        doublereal *z
        );
#endif
#include <lapack.h>
#include "f2c.h"

#define cgeqpb FORTRAN_WRAPPER(cgeqpb)
int cgeqpb(integer *, integer *, integer *, integer *,
        complex *, integer *, complex *, integer *, integer *,
        real *, real *, integer *, real *, complex *, integer *, 
        real *, integer *);

#define cgeqpc FORTRAN_WRAPPER(cgeqpc)
int cgeqpc(integer *, integer *, integer *, integer *, complex *,
        integer *, complex *, integer *, integer *, integer *, real *,
        integer *, integer *, complex *, complex *, real *, real *,
        complex *, integer *, real *);

#define cgeqpw FORTRAN_WRAPPER(cgeqpw)
int cgeqpw(integer *, integer *, integer *, integer *, integer *,
        complex *, integer *, integer *, real *, complex *, real *, real *,
        complex *, complex *, real *);

#define cgeqpx FORTRAN_WRAPPER(cgeqpx)
int cgeqpx(integer*, integer*, integer*, integer*, complex*,
        integer*, complex*, integer*, integer*, real*, real*, integer*,
        real*, complex*, integer*, real*, integer*);

#define cgeqpy FORTRAN_WRAPPER(cgeqpy)
int cgeqpy(integer*, integer*, integer*, integer*, complex*,
        integer*, complex*, integer*, integer*, real*, real*, integer*,
        real*, complex*, integer*, real*, integer*);

#define clasmx FORTRAN_WRAPPER(clasmx)
doublereal clasmx(integer *);

#define clauc1 FORTRAN_WRAPPER(clauc1)
logical clauc1(integer *, complex *, real *, complex *, complex *, real *);

#define ctrqpx FORTRAN_WRAPPER(ctrqpx)
int ctrqpx(integer *, integer *, integer *, integer *, complex *,
        integer *, complex *, integer *, integer *, real *, real *,
        integer *, real *, complex *, real *, integer *);

#define ctrqpy FORTRAN_WRAPPER(ctrqpy)
int ctrqpy(integer *, integer *, integer *, integer *,
        complex *, integer *, complex *c, integer *, integer *,
        real *, real *, integer *, real *, complex *, real *, integer *);

#define ctrqxc FORTRAN_WRAPPER(ctrqxc)
int ctrqxc(integer *, integer *, integer *, integer *, complex *,
        integer *, complex *, integer *, integer *, integer *, real *,
        real *, real *, complex *, real *, integer *);

#define cglbif FORTRAN_WRAPPER(cglbif)
int cglbif(integer *, integer *, integer *, integer *,
        complex *, integer *, complex *, integer *, integer *,
        real *, integer *, logical *, complex *, real *, integer *);

#define ccniif FORTRAN_WRAPPER(ccniif)
int ccniif(integer *, integer *, integer *, integer *,
        complex *, integer *, complex *, integer *, integer *,
        real *, real *, integer *, logical *, complex *, real *, integer *);

#define cgret FORTRAN_WRAPPER(cgret)
int cgret(integer *, integer *, integer *, integer *,
       complex *, integer *, complex *, integer *, complex *,
       real *, integer *);

#define shess FORTRAN_WRAPPER(shess)
int chess(integer *, integer *, integer *, integer *,
        complex *, integer *, complex *, integer *,
       complex *, real *, integer *);

#define ctrqyc FORTRAN_WRAPPER(ctrqyc)
int ctrqyc(integer *, integer *, integer *, integer *,
        complex *, integer *, complex *, integer *, integer *,
        integer *, real *, real *, real *, complex *, real *, integer *);

#define ctrrnk FORTRAN_WRAPPER(ctrrnk)
int ctrrnk(integer *, complex *, integer *, real *, integer *, complex *, integer *);


#define dgeqpb FORTRAN_WRAPPER(dgeqpb)
int dgeqpb(integer *, integer *, integer *,
        integer *, doublereal *, integer *, doublereal *, integer *,
        integer *, doublereal *, doublereal *, integer *, doublereal *,
        doublereal *, integer *, integer *);
  
#define dgeqpc FORTRAN_WRAPPER(dgeqpc)
int dgeqpc(integer *, integer *, integer *,
        integer *, doublereal *, integer *, doublereal *, integer *,
        integer *, integer *, doublereal *, integer *, integer *,
        doublereal *, doublereal *, doublereal *, doublereal *,
        doublereal *, integer *);

#define dgeqpw FORTRAN_WRAPPER(dgeqpw)
int dgeqpw(integer *, integer *, integer *,
        integer *, integer *, doublereal *, integer *, integer *,
        doublereal *, doublereal *, doublereal *, doublereal *,
        doublereal *, doublereal *);

#define dgeqpx FORTRAN_WRAPPER(dgeqpx)
int dgeqpx(integer*, integer*, integer*, integer*, doublereal*,
        integer*, doublereal*, integer*, integer*, doublereal*,
        doublereal*, integer*, doublereal*, doublereal*, integer*,
        integer*);

#define dgeqpy FORTRAN_WRAPPER(dgeqpy)
int dgeqpy(integer*, integer*, integer*, integer*, doublereal*,
        integer*, doublereal*, integer*, integer*, doublereal*,
        doublereal*, integer*, doublereal*, doublereal*, integer*,
        integer*);

#define dlasmx FORTRAN_WRAPPER(dlasmx)
doublereal dlasmx(integer *);

#define dlauc1 FORTRAN_WRAPPER(dlauc1)
logical dlauc1(integer *, doublereal *, doublereal *, doublereal *,
        doublereal *, doublereal *);

#define dtrqpx FORTRAN_WRAPPER(dtrqpx)
int dtrqpx(integer *, integer *, integer *, integer *,
        doublereal *, integer *, doublereal *, integer *, integer *,
        doublereal *, doublereal *, integer *, doublereal *, doublereal *,
        integer *, integer *);

#define dtrqpy FORTRAN_WRAPPER(dtrqpy)
int dtrqpy(integer *, integer *, integer *, integer *,
        doublereal *, integer *, doublereal *, integer *, integer *,
        doublereal *, doublereal *, integer *, doublereal *, doublereal *,
        integer *, integer *);

#define dtrqxc FORTRAN_WRAPPER(dtrqxc)
int dtrqxc(integer *, integer *, integer *, integer *, doublereal *,
        integer*, doublereal *, integer *, integer *, integer *,
        doublereal *, doublereal *, doublereal *, doublereal *, integer *);

#define dcniif FORTRAN_WRAPPER(dcniif)
int dcniif(integer *, integer *, integer *, integer *,
        doublereal *, integer *, doublereal *, integer *, integer *,
        doublereal *, doublereal *, integer *, logical *, doublereal *,
        integer *);

#define dglbif FORTRAN_WRAPPER(dglbif)
int dglbif(integer *, integer *, integer *, integer *,
        doublereal *, integer *, doublereal *, integer *, integer *,
        doublereal *, integer *, logical *, doublereal *, integer *);

#define dgret FORTRAN_WRAPPER(dgret)
int dgret(integer *, integer *, integer *, integer *,
        doublereal *, integer *, doublereal *, integer *,
        doublereal *, integer *);

#define dhess FORTRAN_WRAPPER(dhess)
int dhess(integer *, integer *, integer *, integer *,
       doublereal *, integer *, doublereal *, integer *,
       doublereal *, integer *);

#define dtrqyc FORTRAN_WRAPPER(dtrqyc)
int dtrqyc(integer *, integer *, integer *, integer *, doublereal *,
        integer*, doublereal *, integer *, integer *, integer *,
        doublereal *, doublereal *, doublereal *, doublereal *, integer *);

#define dtrrnk FORTRAN_WRAPPER(dtrrnk)
int dtrrnk(integer *, doublereal *, integer *,
        doublereal *, integer *, doublereal *, integer *);


#define sgeqpb FORTRAN_WRAPPER(sgeqpb)
int sgeqpb(integer *, integer *, integer *, integer *, real *,
        integer *, real *, integer *, integer *, real *, real *, integer *,
        real *, real *, integer *, integer *);

#define sgeqpc FORTRAN_WRAPPER(sgeqpc)
int sgeqpc(integer *, integer *, integer *, integer *, real *,
        integer *, real *, integer *, integer *, integer *, real *,
        integer *, integer *, real *, real *, real *, real *, real *, integer *);

#define sgeqpw FORTRAN_WRAPPER(sgeqpw)
int sgeqpw(integer *, integer *, integer *, integer *,
        integer *, real *, integer *, integer *, real *, real *, real *,
        real *, real *, real *);

#define sgeqpx FORTRAN_WRAPPER(sgeqpx)
int sgeqpx(integer*, integer*, integer*, integer*, real*, integer*,
        real*, integer*, integer*, real*, real*, integer*, real*, real*,
        integer*, integer*);

#define sgeqpy FORTRAN_WRAPPER(sgeqpy)
int sgeqpy(integer*, integer*, integer*, integer*, real*, integer*,
        real*, integer*, integer*, real*, real*, integer*, real*, real*,
        integer*, integer*);

#define slasmx FORTRAN_WRAPPER(slasmx)
doublereal slasmx(integer *);

#define slauc1 FORTRAN_WRAPPER(slauc1)
logical slauc1(integer *, real *, real *, real *, real *, real *);

#define strqpx FORTRAN_WRAPPER(strqpx)
int strqpx(integer *, integer *, integer *, integer *, real *,
        integer *, real *, integer *, integer *, real *, real *, integer *,
        real *, real *, integer *, integer *);

#define strqpy FORTRAN_WRAPPER(strqpy)
int strqpy(integer *, integer *, integer *, integer *, real *,
        integer *, real *, integer *, integer *, real *, real *, integer *,
        real *, real *, integer *, integer *);

#define strqxc FORTRAN_WRAPPER(strqxc)
int strqxc(integer *, integer *, integer *, integer *, real *,
        integer *, real *, integer *, integer *, integer *, real *,
        real *, real *, real *, integer *);

#define sglbif FORTRAN_WRAPPER(sglbif)
int sglbif(integer *, integer *, integer *, integer *,
        real *, integer *, real *, integer *, integer *,
        real *, integer *, logical *, real *, integer *);

#define scniif FORTRAN_WRAPPER(scniif)
int scniif(integer *, integer *, integer *, integer *,
        real *, integer *, real *, integer *, integer *,
        real *, real *, integer *, logical *, real *, integer *);

#define sgret FORTRAN_WRAPPER(sgret)
int sgret(integer *, integer *, integer *, integer *,
        real *, integer *, real *, integer *, real *, integer *);

#define shess FORTRAN_WRAPPER(shess)
int shess(integer *, integer *, integer *, integer *,
        real *, integer *, real *, integer *, real *, integer *);

#define strqyc FORTRAN_WRAPPER(strqyc)
int strqyc(integer *, integer *, integer *, integer *, real *,
        integer *, real *, integer *, integer *, integer *, real *,
        real *, real *, real *, integer *);

#define strrnk FORTRAN_WRAPPER(strrnk)
int strrnk(integer *, real *, integer *, real *,
        integer *, real *, integer *);


#define zgeqpb FORTRAN_WRAPPER(zgeqpb)
int zgeqpb(integer *, integer *, integer *, integer *, doublecomplex *,
        integer *, doublecomplex *, integer *, integer *, doublereal *,
        doublereal *, integer *, doublereal *, doublecomplex *, integer *,
        doublereal *, integer *);

#define zgeqpc FORTRAN_WRAPPER(zgeqpc)
int zgeqpc(integer *, integer *, integer *, integer *, doublecomplex *,
        integer *, doublecomplex *, integer *, integer *, integer *, doublereal *,
        integer *, integer *, doublecomplex *, doublecomplex *, doublereal *, doublereal *,
        doublecomplex *, integer *, doublereal *);

#define zgeqpw FORTRAN_WRAPPER(zgeqpw)
int zgeqpw(integer *, integer *, integer *, integer *, integer *,
        doublecomplex *, integer *, integer *, doublereal *, doublecomplex *,
        doublereal *, doublereal *, doublecomplex *, doublecomplex *, doublereal *);

#define zgeqpx FORTRAN_WRAPPER(zgeqpx)
int zgeqpx(integer*, integer*, integer*, integer*, doublecomplex*,
        integer*, doublecomplex*, integer*, integer*, doublereal*,
        doublereal*, integer*, doublereal*, doublecomplex*, integer*,
        doublereal*, integer*);

#define zgeqpy FORTRAN_WRAPPER(zgeqpy)
int zgeqpy(integer*, integer*, integer*, integer*, doublecomplex*,
        integer*, doublecomplex*, integer*, integer*, doublereal*,
        doublereal*, integer*, doublereal*, doublecomplex*, integer*,
        doublereal*, integer*);

#define zlasmx FORTRAN_WRAPPER(zlasmx)
doublereal zlasmx(integer *i);

#define zlauc1 FORTRAN_WRAPPER(zlauc1)
logical zlauc1(integer *, doublecomplex *, doublereal *,
        doublecomplex *, doublecomplex *, doublereal *);

#define ztrqpx FORTRAN_WRAPPER(ztrqpx)
int ztrqpx(integer *, integer *, integer *, integer *,
        doublecomplex *, integer *, doublecomplex *, integer*,
        integer *, doublereal *, doublereal *, integer *, doublereal *,
        doublecomplex *, doublereal *, integer *);

#define ztrqpy FORTRAN_WRAPPER(ztrqpy)
int ztrqpy(integer *, integer *, integer *, integer *,
        doublecomplex *, integer *, doublecomplex *, integer*,
        integer *, doublereal *, doublereal *, integer *, doublereal *,
        doublecomplex *, doublereal *, integer *);

#define ztrqxc FORTRAN_WRAPPER(ztrqxc)
int ztrqxc(integer *, integer *, integer *, integer *, doublecomplex *,
        integer *, doublecomplex *, integer *, integer *, integer *,
        doublereal *, doublereal *, doublereal *, doublecomplex *,
        doublereal *, integer *);

#define zcniif FORTRAN_WRAPPER(zcniif)
int zcniif(integer *, integer *, integer *, integer *,
        doublecomplex *, integer *, doublecomplex *, integer *, integer *,
        doublereal *, doublereal *, integer *, logical *, doublecomplex *,
        doublereal *, integer *);

#define zglbif FORTRAN_WRAPPER(zglbif)
int zglbif(integer *, integer *, integer *, integer *,
        doublecomplex *, integer *, doublecomplex *, integer *, integer *,
        doublereal *, integer *, logical *, doublecomplex *, doublereal *,
        integer *);
        
#define zgret FORTRAN_WRAPPER(zgret)
int zgret(integer *, integer *, integer *, integer *,
        doublecomplex *, integer *,doublecomplex *, integer *, 
        doublecomplex *, doublereal *, integer *);

#define dzhess FORTRAN_WRAPPER(dhess)
int zhess(integer *, integer *, integer *, integer *,
        doublecomplex *, integer *, doublecomplex *, integer *,
        doublecomplex *, doublereal *, integer *);

#define ztrqyc FORTRAN_WRAPPER(ztrqyc)
int ztrqyc(integer *, integer *, integer *, integer *, doublecomplex *,
        integer *, doublecomplex *, integer *, integer *, integer *,
        doublereal *, doublereal *, doublereal *, doublecomplex *,
        doublereal *, integer *);

#define ztrrnk FORTRAN_WRAPPER(ztrrnk)
int ztrrnk(integer *, doublecomplex *, integer *,
        doublereal *, integer *, doublecomplex *, integer *);
