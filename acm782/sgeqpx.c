/* sgeqpx.f -- translated by f2c (version 20041007).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "rrqr.h"

/* Subroutine */ int sgeqpx(integer *job, integer *m, integer *n, integer *k,
	 real *a, integer *lda, real *c__, integer *ldc, integer *jpvt, real *
	ircond, real *orcond, integer *rank, real *svlues, real *work,
	integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    static integer mn, wkmin;
    static real wsize;


/*     This code is part of a release of the package for computing */
/*     rank-revealing QR Factorizations written by: */
/*     ================================================================== */
/*     Christian H. Bischof        and   Gregorio Quintana-Orti */
/*     Math. and Comp. Sci. Div.         Departamento de Informatica */
/*     Argonne National Lab.             Universidad Jaime I */
/*     Argonne, IL 60439                 Campus P. Roja, 12071 Castellon */
/*     USA                               Spain */
/*     bischof@mcs.anl.gov               gquintan@inf.uji.es */
/*     ================================================================== */
/*     $Revision: 1.84 $ */
/*     $Date: 96/12/30 16:59:13 $ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGEQPX computes a QR factorization */
/*       A*P = Q*[ R11 R12 ] */
/*               [  0  R22 ] */
/*  of a real m by n matrix A. The permutation P is */
/*  chosen with the goal to reveal the rank of A by a */
/*  suitably dimensioned trailing submatrix R22 with norm(R22) */
/*  being small. */

/*  Based on methods related to Chandrasekaran&Ipsen's algorithms. */

/*  Arguments */
/*  ========= */

/*  JOB     (input) INTEGER */
/*          The job to do: */
/*          = 1: The orthogonal transformations needed in the */
/*               triangularization are only applied to matrix A. */
/*               Thus, matrix C is not updated. */
/*          = 2: The same orthogonal transformations needed in the */
/*               triangularization of matrix A are applied to */
/*               matrix C from the left. */
/*               That is, if Q'*A*P=R, then C := Q'*C. */
/*               In this case, matrix C is m-by-k. */
/*          = 3: The transpose of the orthogonal transformations needed */
/*               in the triangularization of matrix A are applied */
/*               to matrix C from the right. */
/*               That is, if Q'*A*P=R, then C := C*Q. */
/*               In this case, matrix C is k-by-m. */
/*          In these three cases, the permutations are always stored */
/*          in vector JPVT. */

/*  M       (input) INTEGER */
/*          The number of rows of matrices A. M >= 0. */
/*          If JOB=2, M is the number of rows of matrix C. */
/*          If JOB=3, M is the number of columns of matrix C. */

/*  N       (input) INTEGER */
/*          The number of columns of matrix A.  N >= 0. */

/*  K       (input) INTEGER */
/*          It defines the dimension of matrix C. K >= 0. */
/*          If JOB=2, K is the number of columns of matrix C. */
/*          If JOB=3, K is the number of rows of matrix C. */

/*  A       (input/output) REAL array, dimension (LDA,N) */
/*          On entry, the m by n matrix A. */
/*          On exit, the upper triangle of the array contains the */
/*          min(m,n) by n upper trapezoidal matrix R; the lower triangle */
/*          array is filled with zeros. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of array A. LDA >= max(1,M). */

/*  C       (input/output) REAL array, dimension */
/*                ( LDC, K ) if JOB=2. */
/*                ( LDC, M ) if JOB=3. */
/*          If argument JOB asks, all the orthogonal transformations */
/*          applied to matrix A are also applied to matrix C. */

/*  LDC     (input) INTEGER */
/*          The leading dimension of array C. */
/*          If JOB=2, then LDC >= MAX(1,M). */
/*          If JOB=3, then LDC >= MAX(1,K). */

/*  JPVT    (output) INTEGER array, dimension (N) */
/*          JPVT(I) = J <==> Column J of A has been permuted into */
/*                           position I in AP. */
/*          JPVT(1:RANK) contains the indices of the columns considered */
/*          linearly independent. */
/*          JPVT(RANK+1:N) contains the indices of the columns considered */
/*          linearly dependent from the previous ones. */

/*  IRCOND  (input) REAL */
/*          1/IRCOND specifies an upper bound on the condition number */
/*          of R11. If IRCOND == 0, IRCOND = machine precision is chosen */
/*          as default. IRCOND must be >= 0. */

/*  ORCOND  (output) REAL */
/*          1/ORCOND is an estimate for the condition number of R11. */

/*  RANK    (output) INTEGER */
/*          RANK is an estimate for the numerical rank of A with respect */
/*          to the threshold 1/IRCOND in the sense that */
/*               RANK = arg_max(cond_no(R(1:r,1:r))<1/IRCOND) */

/*  SVLUES  (output) REAL array, dimension (4) */
/*          On exit, SVLUES contains estimates of some of the */
/*          singular values of the triangular factor R. */
/*          SVLUES(1): largest singular value of R(1:RANK,1:RANK) */
/*          SVLUES(2): smallest singular value of R(1:RANK,1:RANK) */
/*          SVLUES(3): smallest singular value of R(1:RANK+1,1:RANK+1) */
/*          SVLUES(4): smallest singular value of R */
/*          If the triangular factorization is a rank-revealing one */
/*          (which will be the case if the leading columns were well- */
/*          conditioned), then SVLUES(1) will also be an estimate */
/*          for the largest singular value of A, SVLUES(2) and SVLUES(3) */
/*          will be estimates for the RANK-th and (RANK+1)-st singular */
/*          value of A, and SVLUES(4) will be an estimate for the */
/*          smallest singular value of A. */
/*          By examining these values, one can confirm that the rank is */
/*          well defined with respect to the threshold chosen. */

/*  WORK    (workspace) REAL array, dimension (LWORK) */
/*          On exit: work(1) is the size of the storage array needed */
/*                   for optimal performance */

/*  LWORK   (input) INTEGER */
/*          The dimension of array WORK. */
/*          If JOB=1: */
/*             The unblocked strategy requires that: */
/*                 LWORK >= 2*MN+3*N. */
/*             The block algorithm requires that: */
/*                 LWORK >= 2*MN+N*NB. */
/*          If JOB<>1: */
/*             The unblocked strategy requires that: */
/*                 LWORK >= 2*MN+2*N+MAX(K,N). */
/*             The block algorithm requires that: */
/*                 LWORK >= 2*MN+NB*NB+NB*MAX(K,N). */
/*          Where MN = min(M,N) and NB is the block size for this */
/*          environment. */
/*          In both cases, the minimum required workspace is the */
/*          one for the unblocked strategy. */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit. */
/*          < 0: If INFO = -i, the i-th argument had an illegal value */
/*          > 0: Problems in the computation of the rank. */
/*                   1: Exceeded the allowed maximum number of steps. */
/*                   2: Rank not well defined. */
/*          In adition, vector SVLUES tell if rank is not well defined. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --jpvt;
    --svlues;
    --work;

    /* Function Body */
    mn = min(*m,*n);
    if (*job == 1) {
	wkmin = (mn << 1) + *n * 3;
    } else {
	wkmin = (mn << 1) + (*n << 1) + max(*k,*n);
    }

/*     Test input arguments */
/*     ==================== */

    *info = 0;
    if (*job < 1 || *job > 3) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*k < 0) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -6;
    } else if (*job == 1 && *ldc < 1 || *job == 2 && *ldc < max(1,*m) || *job
	    == 3 && *ldc < max(1,*k)) {
	*info = -8;
    } else if (*ircond < 0.f) {
	*info = -10;
    } else if (*lwork < max(1,wkmin)) {
	*info = -15;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("SGEQPX", &i__1, 6);
	return 0;
    }

/*     Preprocessing */
/*     ============= */

    sgeqpb(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &jpvt[1],
	    ircond, orcond, rank, &svlues[1], &work[1], lwork, info);
    wsize = work[1];

/*     Postprocessing */
/*     ============== */

    if (*rank > 0) {
	strqpx(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &jpvt[1]
		, ircond, orcond, rank, &svlues[1], &work[1], lwork, info);
    }

    work[1] = wsize;
    return 0;

/*     End of SGEQPX */

} /* sgeqpx_ */

