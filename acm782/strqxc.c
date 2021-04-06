/* strqxc.f -- translated by f2c (version 20041007).
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

/* Common Block Declarations */

struct {
    integer nb;
} bsprqr_;

#define bsprqr_1 bsprqr_

/* Table of constant values */

static integer c__1 = 1;
static real c_b5 = .5f;
static integer c__2 = 2;

/* Subroutine */ int strqxc(integer *job, integer *m, integer *n, integer *k,
	 real *a, integer *lda, real *c__, integer *ldc, integer *jpvt,
	integer *rank, real *svlues, real *rcnr, real *rcnrp1, real *work,
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer j, mn, ns;
    static real sine, smin, smax;
    static real smxrp1;
    static integer nacptd;
    static real cosine, sminpr, smaxpr;
    static logical permut;
    static integer mxstps;


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
/*     $Date: 96/12/30 16:59:20 $ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STRQXC carries out an algorithm related to algorithm Hybrid-III */
/*  by Chandrasekaran and Ipsen for the stage RANK. The algorithm used */
/*  here offers the following advantages: */
/*  o It is faster since it is based on Chan-II instead of Stewart-II. */
/*  o This algorithm uses the F factor technique to reduce the number of */
/*    cycling problems due to roundoff errors. */
/*  o The final steps that do not improve the ordering are saved. */

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

/*  JPVT    (input/output) INTEGER array, dimension ( N ) */
/*          If JPVT(I) = K, then the Ith column of the permuted */
/*          A was the Kth column of the original A (just before */
/*          the preprocessing). If a permutation occurs, JPVT will */
/*          be updated correctly. */

/*  RANK    (input) INTEGER */
/*          Th estimate of the rank. 1 <= RANK <= MIN(M,N). */

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

/*  RCNR    (output) REAL */
/*          Th estimate for the inverse of the condition number of block */
/*          R(1:RANK,1:RANK). */

/*  RCNRP1  (output) REAL */
/*          Th estimate for the inverse of the condition number of block */
/*          R(1:RANK+1,1:RANK+1). */

/*  WORK    (workspace) REAL array, */
/*             dimension ( MN+MAX(N,2*MN) ), where MN=MIN(M,N). */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit. */
/*          < 0: If info = -i, the i-th argument had an illegal value. */
/*          = 1: Exceeded the allowed maximum number of steps. That is, */
/*               the matrix presents a slow convergence. */


/*  =================================================================== */

/*     .. Parameters .. */

/*     Indices into the `svlues' array. */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
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
    ns = 0;
    mxstps = *n + 25;
    *info = 0;

/*     Quick return if possible. */

    if (mn == 0) {
	return 0;
    }

/*     Inicialization of variable NACPTD, which controls main loop. */

    nacptd = 0;

/*     Compute the norms of block A(1:RANK,1:RANK) and store them */
/*     in vector WORK(1:RANK). It is computed only once at the */
/*     beginning and updated every iteration. It is used to estimate */
/*     the largest singular value in order to pass it to Chan-II. */

    i__1 = *rank;
    for (j = 1; j <= i__1; ++j) {
	work[j] = snrm2(&j, &a[j * a_dim1 + 1], &c__1);
/* L10: */
    }

/*     ***************** */
/*     * start of loop * */
/*     ***************** */

L20:

/*     *-*-*-*-*-*-*-*-*-*-*-*-* */
/*     * call to Golub-I(rank) * */
/*     *-*-*-*-*-*-*-*-*-*-*-*-* */

    if (nacptd < 4) {

/*        Apply Golub-I for the stage RANK. */

	sglbif(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &jpvt[1]
		, &c_b5, rank, &permut, &work[mn + 1], info);

/*        If necessary, update the contents of WORK(RANK). */

	if (permut) {
	    work[*rank] = snrm2(rank, &a[*rank * a_dim1 + 1], &c__1);
	}

/*        Update variables NACPTD and NS. */

	if (permut) {
	    nacptd = 1;
	} else {
	    ++nacptd;
	}
	++ns;
    }

/*     *-*-*-*-*-*-*-*-*-*-*-*-*-* */
/*     * call to Golub-I(rank+1) * */
/*     *-*-*-*-*-*-*-*-*-*-*-*-*-* */

    if (nacptd < 4) {

/*        Determine if the application of Golub-I(rank+1) is necessary. */

	if (*rank == mn) {

/*           Not necessary. Therefore, no permutation occurs. */

	    permut = FALSE_;
	} else {

/*           Apply Golub-I for the stage RANK+1. */

	    i__1 = *rank + 1;
	    sglbif(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &
		    jpvt[1], &c_b5, &i__1, &permut, &work[mn + 1], info);

/*           Update variable NS. */

	    ++ns;
	}

/*        Update variable NACPTD. */

	if (permut) {
	    nacptd = 1;
	} else {
	    ++nacptd;
	}
    }

/*     *-*-*-*-*-*-*-*-*-*-*-*-*-* */
/*     * call to Chan-II (rank+1)* */
/*     *-*-*-*-*-*-*-*-*-*-*-*-*-* */

    if (nacptd < 4) {

/*        Determine if the application of Chan-II(rank+1) is necessary. */

	if (*rank == mn) {

/*           Not necessary. Therefore, no permutation occurs. */

	    permut = FALSE_;
	} else {

/*           Extend vector WORK(1:RANK) to vector WORK(1:RANK+1). */
/*           So, pivoting vector WORK(1:N) inside Chan-II will be */
/*           easier. */

	    i__1 = *rank + 1;
	    work[*rank + 1] = snrm2(&i__1, &a[(*rank + 1) * a_dim1 + 1], &
		    c__1);

/*           Apply Chan-II for the stage RANK+1 */
/*           on block A(1:RANK+1,1:RANK+1). */

	    i__1 = *rank + 1;
	    scniif(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &
		    jpvt[1], &work[1], &c_b5, &i__1, &permut, &work[mn + 1],
		    info);

/*           Update variable NS. */

	    ++ns;
	}

/*        Update variable NACPTD. */

	if (permut) {
	    nacptd = 1;
	} else {
	    ++nacptd;
	}
    }

/*     *-*-*-*-*-*-*-*-*-*-*-*-* */
/*     * call to Chan-II(rank) * */
/*     *-*-*-*-*-*-*-*-*-*-*-*-* */

    if (nacptd < 4) {

/*        Apply Chan-II for the stage RANK on block A(1:RANK,1:RANK). */

	scniif(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &jpvt[1]
		, &work[1], &c_b5, rank, &permut, &work[mn + 1], info);

/*        Update variables NACPTD and NS. */

	if (permut) {
	    nacptd = 1;
	} else {
	    ++nacptd;
	}
	++ns;
    }

/*     Check if loop must finish. */

    if (ns >= mxstps) {
	*info = 1;
    } else if (nacptd < 4) {
	goto L20;
    }

/*     *************** */
/*     * end of loop * */
/*     *************** */

/*     ************************************************************** */
/*     * Computation of vector SVLUES and variables RCNR and RCNRP1 * */
/*     ************************************************************** */

/*     Computation of the largest singular value and the smallest */
/*     singular value of A(1:RANK,1:RANK). */

    smax = (r__1 = a[a_dim1 + 1], dabs(r__1));
    work[1] = 1.f;
    smin = smax;
    work[mn + 1] = 1.f;

    i__1 = *rank;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	slaic1(&c__1, &i__2, &work[1], &smax, &a[j * a_dim1 + 1], &a[j + j *
		a_dim1], &smaxpr, &sine, &cosine);
	i__2 = j - 1;
	sscal(&i__2, &sine, &work[1], &c__1);
	work[j] = cosine;
	smax = smaxpr;
	i__2 = j - 1;
	slaic1(&c__2, &i__2, &work[mn + 1], &smin, &a[j * a_dim1 + 1], &a[j
		+ j * a_dim1], &sminpr, &sine, &cosine);
	i__2 = j - 1;
	sscal(&i__2, &sine, &work[mn + 1], &c__1);
	work[mn + j] = cosine;
	smin = sminpr;
/* L30: */
    }
    svlues[1] = smax;
    svlues[2] = smin;

/*     Computation of the largest singular value and the smallest */
/*     singular value of A(1:RANK+1,1:RANK+1). */

    if (*rank < mn) {
	slaic1(&c__1, rank, &work[1], &smax, &a[(*rank + 1) * a_dim1 + 1], &
		a[*rank + 1 + (*rank + 1) * a_dim1], &smaxpr, &sine, &cosine);
	smax = smaxpr;
	slaic1(&c__2, rank, &work[mn + 1], &smin, &a[(*rank + 1) * a_dim1 +
		1], &a[*rank + 1 + (*rank + 1) * a_dim1], &sminpr, &sine, &
		cosine);
	sscal(rank, &sine, &work[mn + 1], &c__1);
	work[mn + *rank + 1] = cosine;
	smin = sminpr;
    }
    smxrp1 = smax;
    svlues[3] = smin;

/*     Computation of the smallest singular value of A(1:MN,1:MN). */

    i__1 = mn;
    for (j = *rank + 2; j <= i__1; ++j) {
	i__2 = j - 1;
	slaic1(&c__2, &i__2, &work[mn + 1], &smin, &a[j * a_dim1 + 1], &a[j
		+ j * a_dim1], &sminpr, &sine, &cosine);
	i__2 = j - 1;
	sscal(&i__2, &sine, &work[mn + 1], &c__1);
	work[mn + j] = cosine;
	smin = sminpr;
/* L40: */
    }
    svlues[4] = smin;

/*     Computation of RCNR and RCNRP1. */

    *rcnr = svlues[2] / svlues[1];
    *rcnrp1 = svlues[3] / smxrp1;
    return 0;

/*     End of STRQXC */

} /* strqxc_ */

/* Subroutine */ int sglbif(integer *job, integer *m, integer *n, integer *k,
	 real *a, integer *lda, real *c__, integer *ldc, integer *jpvt, real *
	f, integer *rank, logical *permut, real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Local variables */
    static integer j, jj, mn;
    static integer itemp;
    static real rdummy[1]	/* was [1][1] */;



/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGLBIF computes the column index of A(RANK:M,RANK:N) with largest */
/*  norm and determines if pivoting is necessary. If so, it pivots it into */
/*  column RANK, permuts vector JPVT and permuts and retriangularizes */
/*  matrix A. It does only one permutation. */

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

/*  JPVT    (input/output) INTEGER array, dimension ( N ) */
/*          If JPVT(I) = K, then the Ith column of the permuted */
/*          A was the Kth column of the original A (just before the */
/*          preprocessing). If a permutation occurs, it will be */
/*          updated correctly. */

/*  F       (input) REAL */
/*          F factor for the pivoting. It must be always 0 < f <= 1. */

/*  RANK    (input) INTEGER */
/*          Th estimate for the rank. 1 <= RANK <= MIN(M,N). */

/*  PERMUT  (output) LOGICAL */
/*          Tells if a permutation occurred. */

/*  WORK    (workspace) REAL array, */
/*             dimension ( MAX( N, 2*MIN(M,N) ) ) */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit. */
/*          < 0: If info = -i, the i-th argument had an illegal value. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
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
    --work;

    /* Function Body */
    mn = min(*m,*n);
    *info = 0;

/*     Quick return if possible. */

    if (mn == 0 || *rank == *n || *rank == 0) {
	*permut = FALSE_;
	return 0;
    }

/*     Compute the norms of the columns of block A(RANK:M,RANK:N) */
/*     and store them in vector WORK(RANK:N). */

    i__1 = *n;
    for (j = *rank; j <= i__1; ++j) {
	i__2 = min(*m,j) - *rank + 1;
	work[j] = snrm2(&i__2, &a[*rank + j * a_dim1], &c__1);
/* L10: */
    }

/*     Find column with largest two-norm of upper triangular block */
/*     A(RANK:M,RANK:N). Use the data stored in vector WORK(RANK:N). */

    i__1 = *n - *rank + 1;
    jj = *rank - 1 + isamax(&i__1, &work[*rank], &c__1);

/*     Determine if a permutation must occur. */

    *permut = jj > *rank && (r__1 = work[jj], dabs(r__1)) * *f > (r__2 = work[
	    *rank], dabs(r__2));

    if (*permut) {

/*        Exchage cyclically to the right the columns of matrix A */
/*        between RANK and JJ. That is, RANK->RANK+1, */
/*        RANK+1->RANK+2,...,JJ-1->JJ,JJ->K. Use vector WORK(1:MN) */
/*        to store temporal data. */

	i__1 = min(mn,jj);
	scopy(&i__1, &a[jj * a_dim1 + 1], &c__1, &work[1], &c__1);
	i__1 = *rank;
	for (j = jj - 1; j >= i__1; --j) {
/* Computing MIN */
	    i__3 = mn, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    scopy(&i__2, &a[j * a_dim1 + 1], &c__1, &a[(j + 1) * a_dim1 + 1],
		     &c__1);
/* L20: */
	}
	i__1 = min(mn,jj);
	scopy(&i__1, &work[1], &c__1, &a[*rank * a_dim1 + 1], &c__1);

/*        Exchange in the same way vector JPVT. */

	itemp = jpvt[jj];
	i__1 = *rank;
	for (j = jj - 1; j >= i__1; --j) {
	    jpvt[j + 1] = jpvt[j];
/* L30: */
	}
	jpvt[*rank] = itemp;

/*        Retriangularize matrix A after the permutation. */

	if (*job == 1) {
	    i__1 = min(*m,jj) - *rank + 1;
	    i__2 = *n - *rank + 1;
	    sgret(job, &i__1, &i__2, k, &a[*rank + *rank * a_dim1], lda,
		    rdummy, &c__1, &work[1], info);
	} else if (*job == 2) {
	    i__1 = min(*m,jj) - *rank + 1;
	    i__2 = *n - *rank + 1;
	    sgret(job, &i__1, &i__2, k, &a[*rank + *rank * a_dim1], lda, &
		    c__[*rank + c_dim1], ldc, &work[1], info);
	} else if (*job == 3) {
	    i__1 = min(*m,jj) - *rank + 1;
	    i__2 = *n - *rank + 1;
	    sgret(job, &i__1, &i__2, k, &a[*rank + *rank * a_dim1], lda, &
		    c__[*rank * c_dim1 + 1], ldc, &work[1], info);
	}
    }
    return 0;

/*     End of SGLBIF */

} /* sglbif_ */

/* Subroutine */ int scniif(integer *job, integer *m, integer *n, integer *k,
	 real *a, integer *lda, real *c__, integer *ldc, integer *jpvt, real *
	vnorm, real *f, integer *rank, logical *permut, real *work, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer j, jj, mn;
    static real sine, temp, smin, smax;
    static integer itemp;
    static real cosine;
    static real sminpr, rdummy[1];



/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SCNIIF computes the "worst" column in A(1:RANK,1:RANK) and */
/*  determines if pivoting is necessary. If so, it pivots it into column */
/*  RANK, permuts vector JPVT, adjusts vector VNORM and permuts and */
/*  retriangularizes matrix A. It does only one permutation. */

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

/*  JPVT    (input/output) INTEGER array, dimension (N) */
/*          If JPVT(I) = K, then the Ith column of the permuted */
/*          A was the Kth column of the original A (just before the */
/*          preprocessing). If a permutation occurs, this vector will */
/*          be updated correctly. */

/*  VNORM   (input/output) REAL array, dimension ( N ) */
/*          VNORM(1:RANK) contains the norms of the columns of upper */
/*          triangular block A(1:RANK,1:RANK). If a permutation occurs, */
/*          this vector will be updated correctly. */

/*  F       (input) REAL */
/*          F factor for the pivoting. It must be always 0 < f <= 1. */

/*  RANK    (input) INTEGER */
/*          Th estimate for the rank. 1 <= RANK <= MIN(M,N). */

/*  PERMUT  (output) LOGICAL */
/*          Tells if a permutation occurred. */

/*  WORK    (workspace) REAL array, dimension ( 2*MIN(M,N) ) */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit. */
/*          < 0: If info = -i, the i-th argument had an illegal value. */

/*  Further Details */
/*  =============== */

/*    If block R(1:RANK,1:RANK) is singular or near singular, there will */
/*  be no permutation because in that case the right (and left) singular */
/*  vectors are the canonical ones ((0,0,...0,1)^T). */
/*    That is, there will not be permutation if */
/*  RCOND <= SF * SLAMCH('Safe Minimum'), where SF (Safe Factor) is */
/*  a security factor to avoid arithmetic exceptions. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
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
    --vnorm;
    --work;

    /* Function Body */
    mn = min(*m,*n);
    *info = 0;

/*     Quick return if possible. */

    if (mn == 0 || *rank == 1) {
	*permut = FALSE_;
	return 0;
    }

/*     Estimation of the largest singular value of block */
/*     A(1:RANK,1:RANK) by using the contents of vector */
/*     VNORM. */

    itemp = isamax(rank, &vnorm[1], &c__1);
    smax = slasmx(rank) * vnorm[itemp];

/*     Estimation of the smallest singular value of block */
/*     A(1:RANK,1:RANK) and its corresponding left singular vector. */
/*     Save left singular vector in vector WORK(1:RANK). */

    smin = (r__1 = a[a_dim1 + 1], dabs(r__1));
    work[1] = 1.f;
    i__1 = *rank;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	slaic1(&c__2, &i__2, &work[1], &smin, &a[j * a_dim1 + 1], &a[j + j *
		a_dim1], &sminpr, &sine, &cosine);
	i__2 = j - 1;
	sscal(&i__2, &sine, &work[1], &c__1);
	work[j] = cosine;
	smin = sminpr;
/* L10: */
    }

/*     Determine if matrix A(1:RANK,1:RANK) is singular or nearly */
/*     singular. SF (Safe Factor) is used to say if it is singular or not. */

    if (smin <= smax * 100.f * slamch("Safe minimum")) {

/*        Singular or nearly singular matrix. Its right singular */
/*        vector is (0,0,...,0,1)^T. So, no pivoting is needed. */

	*permut = FALSE_;
    } else {

/*        Follow usual method: Estimate the right singular vector */
/*        corresponding to the smallest singular value of upper */
/*        triangular block A(1:RANK,1:RANK) and put into vector */
/*        WORK(1:RANK). */

	strsv("Upper", "No transpose", "No unit", rank, &a[a_offset], lda, &
		work[1], &c__1);

/*        Find the index with largest absolute value in vector */
/*        WORK(1:RANK). */

	jj = isamax(rank, &work[1], &c__1);

/*        Determine if a permutation must occur. */

	*permut = jj < *rank && (r__1 = work[jj], dabs(r__1)) * *f > (r__2 =
		work[*rank], dabs(r__2));

	if (*permut) {

/*           Exchange cyclically to the left the colums of matrix A */
/*           between JJ and RANK. That is, JJ->RANK,JJ+1->JJ,..., */
/*           RANK->RANK-1. Use vector WORK to store temporal data. */

	    scopy(rank, &a[jj * a_dim1 + 1], &c__1, &work[1], &c__1);
	    i__1 = *rank;
	    for (j = jj + 1; j <= i__1; ++j) {
		scopy(&j, &a[j * a_dim1 + 1], &c__1, &a[(j - 1) * a_dim1 + 1]
			, &c__1);
/* L20: */
	    }
	    scopy(rank, &work[1], &c__1, &a[*rank * a_dim1 + 1], &c__1);

/*           Exchange in the same way vector JPVT. */

	    itemp = jpvt[jj];
	    i__1 = *rank;
	    for (j = jj + 1; j <= i__1; ++j) {
		jpvt[j - 1] = jpvt[j];
/* L30: */
	    }
	    jpvt[*rank] = itemp;

/*           Adjust the contents of VNORM. */

	    temp = vnorm[jj];
	    i__1 = *rank;
	    for (j = jj + 1; j <= i__1; ++j) {
		vnorm[j - 1] = vnorm[j];
/* L40: */
	    }
	    vnorm[*rank] = temp;

/*           Retriangularize matrix A after the permutation. */

	    if (*job == 1) {
		i__1 = *rank - jj + 1;
		i__2 = *n - jj + 1;
		shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1], lda,
			rdummy, &c__1, &work[1], info);
	    } else if (*job == 2) {
		i__1 = *rank - jj + 1;
		i__2 = *n - jj + 1;
		shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1], lda, &c__[
			jj + c_dim1], ldc, &work[1], info);
	    } else if (*job == 3) {
		i__1 = *rank - jj + 1;
		i__2 = *n - jj + 1;
		shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1], lda, &c__[
			jj * c_dim1 + 1], ldc, &work[1], info);
	    }
	}
    }
    return 0;

/*     End of SCNIIF */

} /* scniif_ */

/* Subroutine */ int sgret(integer *job, integer *m, integer *n, integer *k,
	real *a, integer *lda, real *c__, integer *ldc, real *work, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j;
    static real r__;
    static integer jb;
    static real sine;
    static integer itemp;
    static real cosine;



/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGRET retriangularizes a special matrix. This has the following */
/*  features: its first column is non-zero and its main diagonal is zero */
/*  except the first element. Now it is showed a 4 by 8 small example: */
/*                           x x x x x x x x */
/*                           x 0 x x x x x x */
/*                           x 0 0 x x x x x */
/*                           x 0 0 0 x x x x */
/*  This subroutine assumes that in all cases N>=M. */
/*  The orthogonal transformations applied to matrix A can be also */
/*  applied to matrix C. */

/*  Parameters */
/*  ========== */

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

/*  WORK    (workspace) REAL array, dimension ( 2*M ) */
/*          If the block algorithm is used, the size of this workspace */
/*          must be ( 2*M ). */
/*          In this case this vector will contain the sines and cosines */
/*          of the angles of the Givens Rotations to be applied. */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit. */
/*          < 0: If info = -i, the i-th argument had an illegal value. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Common Blocks .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;

/*     Quick return if possible. */

    if (*m == 0 || *m == 1 || *n == 0) {
	return 0;
    }
    if (bsprqr_1.nb > 1) {

/*        Block Algorithm */
/*        =============== */

/*        Compute Givens Rotations needed to nullify the first column */
/*        of matrix A and apply on the fly to that column. Store */
/*        temporally the sine and cosine of the angles of the applied */
/*        Givens Rotations in vector WORK. */

	for (i__ = *m; i__ >= 2; --i__) {
	    slartg(&a[i__ - 1 + a_dim1], &a[i__ + a_dim1], &work[i__], &work[
		    *m + i__], &r__);
	    a[i__ - 1 + a_dim1] = r__;
	    a[i__ + a_dim1] = 0.f;
/* L10: */
	}

/*        Apply the previously computed Givens Rotations to the rest */
/*        of matrix A. */

	i__1 = *n;
	i__2 = bsprqr_1.nb;
	for (j = 2; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = bsprqr_1.nb, i__4 = *n - j + 1;
	    jb = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *m, i__4 = j + jb - 1;
	    i__5 = j;
	    for (i__ = min(i__3,i__4); i__ >= i__5; --i__) {
		i__3 = j + jb - i__;
		srot(&i__3, &a[i__ - 1 + i__ * a_dim1], lda, &a[i__ + i__ *
			a_dim1], lda, &work[i__], &work[*m + i__]);
/* L30: */
	    }
/* Computing MIN */
	    i__5 = *m, i__3 = j - 1;
	    for (i__ = min(i__5,i__3); i__ >= 2; --i__) {
		srot(&jb, &a[i__ - 1 + j * a_dim1], lda, &a[i__ + j * a_dim1]
			, lda, &work[i__], &work[*m + i__]);
/* L40: */
	    }
/* L20: */
	}

/*        Update the corresponding part of matrix C. */

	if (*job == 2 && *k > 0) {

/*           Apply the previously computed rotations from the left. */

	    i__2 = *k;
	    i__1 = bsprqr_1.nb;
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
/* Computing MIN */
		i__5 = bsprqr_1.nb, i__3 = *k - j + 1;
		jb = min(i__5,i__3);
		for (i__ = *m; i__ >= 2; --i__) {
		    srot(&jb, &c__[i__ - 1 + j * c_dim1], ldc, &c__[i__ + j *
			     c_dim1], ldc, &work[i__], &work[*m + i__]);
/* L60: */
		}
/* L50: */
	    }
	} else if (*job == 3 && *k > 0) {

/*           Apply the transpose of the previously computed rotations */
/*           from the right. */

	    for (i__ = *m; i__ >= 2; --i__) {
		srot(k, &c__[(i__ - 1) * c_dim1 + 1], &c__1, &c__[i__ *
			c_dim1 + 1], &c__1, &work[i__], &work[*m + i__]);
/* L70: */
	    }
	}
    } else {

/*        Non-Block Algorithm */
/*        =================== */

	for (i__ = *m; i__ >= 2; --i__) {
	    itemp = i__ - 1;

/*           Compute the rotation parameters and update column 1 of A. */

	    slartg(&a[itemp + a_dim1], &a[i__ + a_dim1], &cosine, &sine, &
		    r__);
	    a[itemp + a_dim1] = r__;
	    a[i__ + a_dim1] = 0.f;

/*           Update columns I:N of matrix A. */

	    i__1 = *n - i__ + 1;
	    srot(&i__1, &a[itemp + i__ * a_dim1], lda, &a[i__ + i__ * a_dim1]
		    , lda, &cosine, &sine);

/*           Update the corresponding part of matrix C. */

	    if (*job == 2 && *k > 0) {

/*              Apply the previously computed rotations from the left. */

		srot(k, &c__[itemp + c_dim1], ldc, &c__[i__ + c_dim1], ldc, &
			cosine, &sine);
	    } else if (*job == 3 && *k > 0) {

/*              Apply the transpose of the previously computed rotations */
/*              from the right. */

		srot(k, &c__[itemp * c_dim1 + 1], &c__1, &c__[i__ * c_dim1 +
			1], &c__1, &cosine, &sine);
	    }
/* L90: */
	}
    }
    return 0;

/*     End of SGRET */

} /* sgret_ */

/* Subroutine */ int shess(integer *job, integer *m, integer *n, integer *k,
	real *a, integer *lda, real *c__, integer *ldc, real *work, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j;
    static real r__;
    static integer jb;
    static real sine;
    static integer itemp;
    static real cosine;



/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SHESS reduces the upper hessemberg matrix A to upper triangular form. */
/*  applied to matrix C if argument JOB asks. */
/*  This subroutine assumes that in all cases N>=M. */

/*  Parameters */
/*  ========== */

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

/*  WORK    (workspace) REAL array, dimension ( 2*M ) */
/*          If the block algorithm is used, the size of this workspace */
/*          must be ( 2*M ). */
/*          In this case this vector will contain the sines and cosines */
/*          of the angles of the Givens Rotations to be applied. */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit. */
/*          < 0: If info = -i, the i-th argument had an illegal value. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Common Blocks .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;

/*     Quick return if possible. */

    if (*m == 0 || *m == 1 || *n == 0) {
	return 0;
    }
    if (bsprqr_1.nb > 1) {

/*        Block Algorithm */
/*        =============== */

/*        Compute Givens Rotations needed to reduce upper hessenberg */
/*        matrix A to triangular form, and apply on the fly those */
/*        rotations to matrix. Store temporally the sine and cosine */
/*        of the angles of the applied Givens Rotations in vector WORK. */

	i__1 = *n;
	i__2 = bsprqr_1.nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = bsprqr_1.nb, i__4 = *n - j + 1;
	    jb = min(i__3,i__4);
	    i__3 = min(*m,j);
	    for (i__ = 2; i__ <= i__3; ++i__) {
		srot(&jb, &a[i__ - 1 + j * a_dim1], lda, &a[i__ + j * a_dim1]
			, lda, &work[i__], &work[*m + i__]);
/* L20: */
	    }
/* Computing MIN */
	    i__4 = *m, i__5 = j + jb;
	    i__3 = min(i__4,i__5);
	    for (i__ = j + 1; i__ <= i__3; ++i__) {
		itemp = i__ - 1;
		slartg(&a[itemp + itemp * a_dim1], &a[i__ + itemp * a_dim1],
			&work[i__], &work[*m + i__], &r__);
		a[itemp + itemp * a_dim1] = r__;
		a[i__ + itemp * a_dim1] = 0.f;
		i__4 = j + jb - i__;
		srot(&i__4, &a[itemp + i__ * a_dim1], lda, &a[i__ + i__ *
			a_dim1], lda, &work[i__], &work[*m + i__]);
/* L30: */
	    }
/* L10: */
	}

/*        Update the corresponding part of matrix C. */

	if (*job == 2 && *k > 0) {

/*           Apply the previously computed rotations from the left. */

	    i__2 = *k;
	    i__1 = bsprqr_1.nb;
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
/* Computing MIN */
		i__3 = bsprqr_1.nb, i__4 = *k - j + 1;
		jb = min(i__3,i__4);
		i__3 = *m;
		for (i__ = 2; i__ <= i__3; ++i__) {
		    srot(&jb, &c__[i__ - 1 + j * c_dim1], ldc, &c__[i__ + j *
			     c_dim1], ldc, &work[i__], &work[*m + i__]);
/* L50: */
		}
/* L40: */
	    }
	} else if (*job == 3 && *k > 0) {

/*           Apply the transpose of the previously computed rotations */
/*           from the right. */

	    i__1 = *m;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		srot(k, &c__[(i__ - 1) * c_dim1 + 1], &c__1, &c__[i__ *
			c_dim1 + 1], &c__1, &work[i__], &work[*m + i__]);
/* L60: */
	    }
	}
    } else {

/*        Non-Block Algorithm */
/*        =================== */

	i__1 = *m;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    itemp = i__ - 1;

/*           Compute the rotation parameters. */

	    slartg(&a[itemp + itemp * a_dim1], &a[i__ + itemp * a_dim1], &
		    cosine, &sine, &r__);

/*           Update columns I-1:N of matrix A. */

	    a[itemp + itemp * a_dim1] = r__;
	    a[i__ + itemp * a_dim1] = 0.f;
	    i__2 = *n - i__ + 1;
	    srot(&i__2, &a[itemp + i__ * a_dim1], lda, &a[i__ + i__ * a_dim1]
		    , lda, &cosine, &sine);

/*           Update the corresponding part of matrix C. */

	    if (*job == 2 && *k > 0) {

/*              Apply the previously computed rotations from the left. */

		srot(k, &c__[itemp + c_dim1], ldc, &c__[i__ + c_dim1], ldc, &
			cosine, &sine);
	    } else if (*job == 3 && *k > 0) {

/*              Apply the transpose of the previously computed rotations */
/*              from the right. */

		srot(k, &c__[itemp * c_dim1 + 1], &c__1, &c__[i__ * c_dim1 +
			1], &c__1, &cosine, &sine);
	    }
/* L80: */
	}
    }
    return 0;

/*     End of SHESS */

} /* shess_ */

