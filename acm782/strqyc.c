/* strqyc.f -- translated by f2c (version 20041007).
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

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int strqyc(integer *job, integer *m, integer *n, integer *k,
	 real *a, integer *lda, real *c__, integer *ldc, integer *jpvt,
	integer *rank, real *svlues, real *rcnr, real *rcnrp1, real *work,
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real f;
    static integer i__, j, ii, jj, mn, ns, nca;
    static real diag, sine, temp, smin, smax;
    static integer nctba;
    static integer itemp;
    static real smnrp1, smxrp1;
    static real cosine;
    static real sminpr, smaxpr, rdummy[1]	/* was [1][1] */;
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
/*     $Date: 96/12/30 16:59:21 $ */

/*     .. Scalars Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STRQYC carries out Pan-Tang Algorithm 3 for the stage RANK. */
/*  This is a mofified version of the original algorithm. The improved */
/*  features are the following: */
/*  o Use of Bischof's ICE to reduce the computational cost. */
/*  o Reorganization of the main loop to save computations. */
/*  o No permutation is carried out if not strictly needed. */

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
/*          The number of rows of matrices A and C. M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of matrix A. N >= 0. */

/*  K       (input) INTEGER */
/*          The number of columns of matrix C. K >= 0. */

/*  A       (input/output) REAL array (LDA,N) */
/*          Upper triangular m by n matrix. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of array A. LDA >= MAX( 1, M ). */

/*  C       (input/output) REAL array (LDC,K) */
/*          Matrix of dimension m x k where to accumulate */
/*          orthogonal transformations from the left. */

/*  LDC     (input) INTEGER */
/*          The leading dimension of array C. LDC >= MAX( 1, M ). */

/*  JPVT    (input/output) INTEGER array (N) */
/*          Vector with the actual permutation of matrix A. */

/*  RANK    (input) INTEGER */
/*          The estimate for the rank. 1 <= RANK <= MIN(M,N). */

/*  SVLUES  (output) REAL array, dimension (4) */
/*          On exit, SVLUES contains estimates of some of the singular */
/*          values of the triangular factor R. */
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
/*          The estimate for the inverse of the condition number of */
/*          block R(1:RANK,1:RANK). */

/*  RCNRP1  (output) REAL */
/*          The estimate for the inverse of the condition number of */
/*          block R(1:RANK+1,1:RANK+1). */

/*  WORK    (workspace) REAL array, dimension (N+3*MIN(M,N)) */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit. */
/*          < 0: If info = -i, the i-th argument had an illegal value. */
/*          = 4: Exceeded the allowed maximum number of steps. That is, */
/*               the matrix presents a slow convergence. */


/*  Further Details */
/*  =============== */

/*    If the leading block of R is singular or near singular, there will */
/*  be no permutation because in that case the right (and left) singular */
/*  vectors are the canonical ones ((0,0,...0,1)^T). */
/*    That is, there will not be permutation if */
/*  RCOND <= SF * SLAMCH('Safe Minimum'), where SF (Safe Factor) is */
/*  a security factor to avoid arithmetic exceptions. */

/*  ===================================================================== */

/*     .. Parameters .. */

/*     Indices into the `svlues' array. */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
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
    mxstps = *n + 25;
    ns = 0;

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
    } else if (*rank < 1 || *rank > mn) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("STRQYC", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (mn == 0) {
	return 0;
    }

    if (*rank == mn) {

/*        ************************ */
/*        ************************ */
/*        * Apply Chan Algorithm * */
/*        ************************ */
/*        ************************ */

	f = .9f;

/*        Move the best column of A(1:M,M:N) to position M-th. */

	i__1 = *n - mn + 1;
	jj = mn - 1 + isamax(&i__1, &a[mn + mn * a_dim1], lda);
	if (jj > mn) {
	    sswap(m, &a[mn * a_dim1 + 1], &c__1, &a[jj * a_dim1 + 1], &c__1);
	    itemp = jpvt[mn];
	    jpvt[mn] = jpvt[jj];
	    jpvt[jj] = itemp;
	}

/*        Estimation of the largest singular value, the smallest singular */
/*        value, and its corresponding left singular vector. */

	smax = (r__1 = a[a_dim1 + 1], dabs(r__1));
	work[1] = 1.f;
	smin = smax;
	work[mn + 1] = 1.f;
	i__1 = *rank;
	for (j = 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    slaic1(&c__1, &i__2, &work[1], &smax, &a[j * a_dim1 + 1], &a[j +
		    j * a_dim1], &smaxpr, &sine, &cosine);
	    i__2 = j - 1;
	    sscal(&i__2, &sine, &work[1], &c__1);
	    work[j] = cosine;
	    smax = smaxpr;

	    i__2 = j - 1;
	    slaic1(&c__2, &i__2, &work[mn + 1], &smin, &a[j * a_dim1 + 1], &
		    a[j + j * a_dim1], &sminpr, &sine, &cosine);
	    i__2 = j - 1;
	    sscal(&i__2, &sine, &work[mn + 1], &c__1);
	    work[mn + j] = cosine;
	    smin = sminpr;
/* L10: */
	}

/*        Determine if matrix A is singular or nearly singular. */
/*        SF (Safe Factor) is used to say whether or not it is. */

	if (smin > smax * 100.f * slamch("Safe Minimum")) {

/*           Matrix is not singular or not nearly singular. */
/*           Follow usual method: Estimate the right singular vector */
/*           corresponding to the smallest singular value of upper */
/*           triangular block A(1:RANK,1:RANK). */

	    strsv("Upper", "No transpose", "No unit", rank, &a[a_offset],
		    lda, &work[mn + 1], &c__1);

/*           Find the index with largest absolute value in vector */
/*           WORK( MN+1:2*MN ). */

	    jj = isamax(rank, &work[mn + 1], &c__1);

/*           Permut if necessary. */

	    if (jj < *rank && (r__1 = work[mn + jj], dabs(r__1)) * f > (r__2 =
		     work[mn + *rank], dabs(r__2))) {

		ns = 1;

/*              Exchange cyclically to the left the columns of A between */
/*              JJ and RANK, that is: JJ->RANK, JJ+1->JJ, JJ+2->JJ+1,..., */
/*              RANK->RANK-1. */

		scopy(rank, &a[jj * a_dim1 + 1], &c__1, &work[1], &c__1);
		i__1 = *rank;
		for (j = jj + 1; j <= i__1; ++j) {
		    scopy(&j, &a[j * a_dim1 + 1], &c__1, &a[(j - 1) * a_dim1
			    + 1], &c__1);
/* L20: */
		}
		scopy(rank, &work[1], &c__1, &a[*rank * a_dim1 + 1], &c__1);

/*              Exchange in the same way vector JPVT. */

		itemp = jpvt[jj];
		i__1 = *rank;
		for (j = jj + 1; j <= i__1; ++j) {
		    jpvt[j - 1] = jpvt[j];
/* L30: */
		}
		jpvt[*rank] = itemp;

/*              Retriangularization of matrix A after the permutation. */

		if (*job == 1) {
		    i__1 = *rank - jj + 1;
		    i__2 = *n - jj + 1;
		    shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1], lda,
			    rdummy, &c__1, &work[1], info);
		} else if (*job == 2) {
		    i__1 = *rank - jj + 1;
		    i__2 = *n - jj + 1;
		    shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1], lda, &
			    c__[jj + c_dim1], ldc, &work[1], info);
		} else if (*job == 3) {
		    i__1 = *rank - jj + 1;
		    i__2 = *n - jj + 1;
		    shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1], lda, &
			    c__[jj * c_dim1 + 1], ldc, &work[1], info);
		}
	    }
	}

/*        Computation of the contents of vector SVLUES, RCNR and RCNRP1. */

	svlues[1] = smax;
	svlues[2] = smin;
	svlues[3] = smin;
	svlues[4] = smin;
	*rcnr = svlues[2] / svlues[1];
	*rcnrp1 = *rcnr;
    } else {

/*        *************************************** */
/*        *************************************** */
/*        * Apply Modified Pan&Tang Algorithm 3 * */
/*        *************************************** */
/*        *************************************** */

/*        Adjust the value of f. */

	f = .9f / sqrt((real) (*rank + 1));

/*        Compute the norms of columns of matrix A. Store them into */
/*        vector WORK(1:N). */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = min(*m,j);
	    work[j] = snrm2(&i__2, &a[j * a_dim1 + 1], &c__1);
/* L100: */
	}

/*        Estimate the smallest singular value of A(1:RANK,1:RANK) and */
/*        its corresponding left singular vector. */
/*        SMIN will contain the smallest singular value and */
/*        WORK(N+1:N+MN) will contain the left singular vector. */

	smin = (r__1 = a[a_dim1 + 1], dabs(r__1));
	work[*n + 1] = 1.f;
	i__1 = *rank;
	for (j = 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    slaic1(&c__2, &i__2, &work[*n + 1], &smin, &a[j * a_dim1 + 1], &
		    a[j + j * a_dim1], &sminpr, &sine, &cosine);
	    i__2 = j - 1;
	    sscal(&i__2, &sine, &work[*n + 1], &c__1);
	    work[*n + j] = cosine;
	    smin = sminpr;
/* L110: */
	}

/*        Initialize loop variables. */

	nca = 0;
	nctba = *n - *rank;
	ii = *rank + 1;

/*        *********************** */
/*        * Start of Loop WHILE * */
/*        *********************** */

L1000:
	if (nca < nctba && ns < mxstps) {

/*           Estimate the smallest singular value of A(1:RANK+1,1:RANK+1) */
/*           and its corresponding left singular vector as if column II */
/*           of matrix A were on column RANK+1. */

	    diag = a[min(mn,ii) + ii * a_dim1];
	    i__1 = *rank + 1;
	    for (i__ = min(mn,ii) - 1; i__ >= i__1; --i__) {
		slartg(&a[i__ + ii * a_dim1], &diag, &cosine, &sine, &temp);
		diag = temp;
/* L120: */
	    }

	    slaic1(&c__2, rank, &work[*n + 1], &smin, &a[ii * a_dim1 + 1], &
		    diag, &smnrp1, &sine, &cosine);
	    if (smnrp1 >= f * dabs(diag)) {

/*              Column II accepted on the right part of matrix A. */

		++nca;
		if (ii == *n) {
		    ii = *rank + 1;
		} else {
		    ++ii;
		}
	    } else {

/*              Column II not accepted on the right part of matrix A. */

/*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */
/*              * Permut column II to position RANK+1 * */
/*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */

/*              Exchange cyclically to the right the columns of A between */
/*              RANK+1 and II, that is, RANK+1->RANK+2, RANK+2->RANK+3, */
/*              ...,II-1->II,II->RANK+1. */

		i__1 = min(mn,ii);
		scopy(&i__1, &a[ii * a_dim1 + 1], &c__1, &work[*n + mn + 1],
			&c__1);
		i__1 = *rank + 1;
		for (j = ii - 1; j >= i__1; --j) {
/* Computing MIN */
		    i__3 = mn, i__4 = j + 1;
		    i__2 = min(i__3,i__4);
		    scopy(&i__2, &a[j * a_dim1 + 1], &c__1, &a[(j + 1) *
			    a_dim1 + 1], &c__1);
/* L130: */
		}
		i__1 = min(mn,ii);
		scopy(&i__1, &work[*n + mn + 1], &c__1, &a[(*rank + 1) *
			a_dim1 + 1], &c__1);

/*              Exchange in the same way vector JPVT. */

		itemp = jpvt[ii];
		i__1 = *rank + 1;
		for (j = ii - 1; j >= i__1; --j) {
		    jpvt[j + 1] = jpvt[j];
/* L140: */
		}
		jpvt[*rank + 1] = itemp;

/*              Exchange in the same way vector WORK(1:N). */

		temp = work[ii];
		i__1 = *rank + 1;
		for (j = ii - 1; j >= i__1; --j) {
		    work[j + 1] = work[j];
/* L150: */
		}
		work[*rank + 1] = temp;

/*              Retriangularize matrix A after permutation. */

		if (*job == 1) {
		    i__1 = min(*m,ii) - *rank;
		    i__2 = *n - *rank;
		    sgret(job, &i__1, &i__2, k, &a[*rank + 1 + (*rank + 1) *
			    a_dim1], lda, rdummy, &c__1, &work[*n + mn + 1],
			    info);
		} else if (*job == 2) {
		    i__1 = min(*m,ii) - *rank;
		    i__2 = *n - *rank;
		    sgret(job, &i__1, &i__2, k, &a[*rank + 1 + (*rank + 1) *
			    a_dim1], lda, &c__[*rank + 1 + c_dim1], ldc, &
			    work[*n + mn + 1], info);
		} else if (*job == 3) {
		    i__1 = min(*m,ii) - *rank;
		    i__2 = *n - *rank;
		    sgret(job, &i__1, &i__2, k, &a[*rank + 1 + (*rank + 1) *
			    a_dim1], lda, &c__[(*rank + 1) * c_dim1 + 1], ldc,
			     &work[*n + mn + 1], info);
		}

/*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */
/*              * Estimate the largest singular value * */
/*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */

		i__1 = *rank + 1;
		itemp = isamax(&i__1, &work[1], &c__1);
		i__1 = *rank + 1;
		smxrp1 = slasmx(&i__1) * work[itemp];

/*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */
/*              * Estimate the right singular vector  * */
/*              *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */

		if (smnrp1 > smxrp1 * 100.f * slamch("Safe minimum")) {

/*                 Matrix is not singular or not nearly singular. */

/*                 First, end the estimation of the left singular vector. */

		    scopy(rank, &work[*n + 1], &c__1, &work[*n + mn + 1], &
			    c__1);
		    sscal(rank, &sine, &work[*n + mn + 1], &c__1);
		    work[*n + mn + *rank + 1] = cosine;

/*                 Obtain the right singular vector from the left one. */

		    i__1 = *rank + 1;
		    strsv("Upper", "No transpose", "No unit", &i__1, &a[
			    a_offset], lda, &work[*n + mn + 1], &c__1);

		    i__1 = *rank + 1;
		    jj = isamax(&i__1, &work[*n + mn + 1], &c__1);

/*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */
/*                 * Permut column JJ to position RANK+1 * */
/*                 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */

		    if (jj < *rank + 1) {

/*                    Exchange cyclically to the left the columns of A */
/*                    between JJ and RANK+1, that is, JJ->RANK+1,JJ+1->JJ, */
/*                    JJ+2->JJ+1,...,RANK+1->RANK. */

			i__1 = *rank + 1;
			scopy(&i__1, &a[jj * a_dim1 + 1], &c__1, &work[*n +
				mn + 1], &c__1);
			i__1 = *rank + 1;
			for (j = jj + 1; j <= i__1; ++j) {
			    scopy(&j, &a[j * a_dim1 + 1], &c__1, &a[(j - 1) *
				     a_dim1 + 1], &c__1);
/* L160: */
			}
			i__1 = *rank + 1;
			scopy(&i__1, &work[*n + mn + 1], &c__1, &a[(*rank +
				1) * a_dim1 + 1], &c__1);

/*                    Exchange in the same way vector JPVT. */

			itemp = jpvt[jj];
			i__1 = *rank + 1;
			for (j = jj + 1; j <= i__1; ++j) {
			    jpvt[j - 1] = jpvt[j];
/* L170: */
			}
			jpvt[*rank + 1] = itemp;

/*                    Exchange in the same way vector WORK. */

			temp = work[jj];
			i__1 = *rank + 1;
			for (j = jj + 1; j <= i__1; ++j) {
			    work[j - 1] = work[j];
/* L180: */
			}
			work[*rank + 1] = temp;

/*                    Retriangularize matrix A after the permutation. */

			if (*job == 1) {
			    i__1 = *rank - jj + 2;
			    i__2 = *n - jj + 1;
			    shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1],
				     lda, rdummy, &c__1, &work[*n + mn + 1],
				    info);
			} else if (*job == 2) {
			    i__1 = *rank - jj + 2;
			    i__2 = *n - jj + 1;
			    shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1],
				     lda, &c__[jj + c_dim1], ldc, &work[*n +
				    mn + 1], info);
			} else if (*job == 3) {
			    i__1 = *rank - jj + 2;
			    i__2 = *n - jj + 1;
			    shess(job, &i__1, &i__2, k, &a[jj + jj * a_dim1],
				     lda, &c__[jj * c_dim1 + 1], ldc, &work[*
				    n + mn + 1], info);
			}

/*                    Estimate the smallest singular value of */
/*                    A(1:RANK,1:RANK) and its corresponding left */
/*                    singular vector. */
/*                    SMIN will contain the smallest singular value and */
/*                    WORK(N+1:N+MN) will contain the left singular */
/*                    vector. */

			smin = (r__1 = a[a_dim1 + 1], dabs(r__1));
			work[*n + 1] = 1.f;
			i__1 = *rank;
			for (j = 2; j <= i__1; ++j) {
			    i__2 = j - 1;
			    slaic1(&c__2, &i__2, &work[*n + 1], &smin, &a[j *
				     a_dim1 + 1], &a[j + j * a_dim1], &sminpr,
				     &sine, &cosine);
			    i__2 = j - 1;
			    sscal(&i__2, &sine, &work[*n + 1], &c__1);
			    work[*n + j] = cosine;
			    smin = sminpr;
/* L190: */
			}
		    }
		}

/*              Update loop variables. */

		nca = 0;
		++ns;
		if (ii == *n) {
		    ii = *rank + 1;
		} else {
		    ++ii;
		}
	    }
	    goto L1000;
	}

/*        ********************* */
/*        * End of Loop WHILE * */
/*        ********************* */

/*        ****************** */
/*        * Final Pivoting * */
/*        ****************** */

/*        Exchange column in R(RANK+1:M,RANK+1:N) with largest norm to */
/*        position RANK+1. */

	i__1 = *n - *rank;
	jj = *rank + isamax(&i__1, &work[*rank + 1], &c__1);
	if (jj > *rank + 1 && f * (r__1 = work[jj], dabs(r__1)) > (r__2 =
		work[*rank + 1], dabs(r__2))) {

/*           Exchange column JJ to position RANK+1. */

	    i__1 = min(mn,jj);
	    scopy(&i__1, &a[jj * a_dim1 + 1], &c__1, &work[*n + mn + 1], &
		    c__1);
	    i__1 = *rank + 1;
	    for (j = jj - 1; j >= i__1; --j) {
/* Computing MIN */
		i__3 = mn, i__4 = j + 1;
		i__2 = min(i__3,i__4);
		scopy(&i__2, &a[j * a_dim1 + 1], &c__1, &a[(j + 1) * a_dim1
			+ 1], &c__1);
/* L200: */
	    }
	    i__1 = min(mn,jj);
	    scopy(&i__1, &work[*n + mn + 1], &c__1, &a[(*rank + 1) * a_dim1
		    + 1], &c__1);

/*           Exchange in the same way vector JPVT. */

	    itemp = jpvt[jj];
	    i__1 = *rank + 1;
	    for (j = jj - 1; j >= i__1; --j) {
		jpvt[j + 1] = jpvt[j];
/* L210: */
	    }
	    jpvt[*rank + 1] = itemp;

/*           Retriangularize matrix A after permutation. */

	    if (*job == 1) {
		i__1 = min(*m,jj) - *rank;
		i__2 = *n - *rank;
		sgret(job, &i__1, &i__2, k, &a[*rank + 1 + (*rank + 1) *
			a_dim1], lda, rdummy, &c__1, &work[*n + mn + 1], info)
			;
	    } else if (*job == 2) {
		i__1 = min(*m,jj) - *rank;
		i__2 = *n - *rank;
		sgret(job, &i__1, &i__2, k, &a[*rank + 1 + (*rank + 1) *
			a_dim1], lda, &c__[*rank + 1 + c_dim1], ldc, &work[*n
			+ mn + 1], info);
	    } else if (*job == 3) {
		i__1 = min(*m,jj) - *rank;
		i__2 = *n - *rank;
		sgret(job, &i__1, &i__2, k, &a[*rank + 1 + (*rank + 1) *
			a_dim1], lda, &c__[(*rank + 1) * c_dim1 + 1], ldc, &
			work[*n + mn + 1], info);
	    }
	}

/*        ************************************************************** */
/*        * Computation of vector SVLUES and variables RCNR and RCNRP1 * */
/*        ************************************************************** */

/*        Computation of the largest singular value of A(1:RANK,1:RANK). */

	smax = (r__1 = a[a_dim1 + 1], dabs(r__1));
	work[1] = 1.f;

	i__1 = *rank;
	for (j = 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    slaic1(&c__1, &i__2, &work[1], &smax, &a[j * a_dim1 + 1], &a[j +
		    j * a_dim1], &smaxpr, &sine, &cosine);
	    i__2 = j - 1;
	    sscal(&i__2, &sine, &work[1], &c__1);
	    work[j] = cosine;
	    smax = smaxpr;
/* L220: */
	}
	svlues[1] = smax;
	svlues[2] = smin;

/*        Computation of the largest singular value and the smallest */
/*        singular value of A(1:RANK+1,1:RANK+1). */

	slaic1(&c__1, rank, &work[1], &smax, &a[(*rank + 1) * a_dim1 + 1], &
		a[*rank + 1 + (*rank + 1) * a_dim1], &smxrp1, &sine, &cosine);
	slaic1(&c__2, rank, &work[*n + 1], &smin, &a[(*rank + 1) * a_dim1 +
		1], &a[*rank + 1 + (*rank + 1) * a_dim1], &sminpr, &sine, &
		cosine);
	sscal(rank, &sine, &work[*n + 1], &c__1);
	work[*n + *rank + 1] = cosine;
	smin = sminpr;
	svlues[3] = smin;

/*        Computation of the smallest singular value of A(1:MN,1:MN). */

	i__1 = mn;
	for (j = *rank + 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    slaic1(&c__2, &i__2, &work[*n + 1], &smin, &a[j * a_dim1 + 1], &
		    a[j + j * a_dim1], &sminpr, &sine, &cosine);
	    i__2 = j - 1;
	    sscal(&i__2, &sine, &work[*n + 1], &c__1);
	    work[*n + j] = cosine;
	    smin = sminpr;
/* L230: */
	}
	svlues[4] = smin;

/*        Computation of RCNR and RCNRP1. */

	*rcnr = svlues[2] / svlues[1];
	*rcnrp1 = svlues[3] / smxrp1;

	if (ns >= mxstps) {
	    *info = 1;
	}
    }
    return 0;

/*     End of STRQYC */

} /* strqyc_ */

