/* dtrqpy.f -- translated by f2c (version 20041007).
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
static integer c__0 = 0;

/* Subroutine */ int dtrqpy(integer *job, integer *m, integer *n, integer *k,
	 doublereal *a, integer *lda, doublereal *c__, integer *ldc, integer *
	jpvt, doublereal *ircond, doublereal *orcond, integer *rank,
	doublereal *svlues, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer mn;
    static doublereal rcnr, rcond;
    static integer oinfo;
    static doublereal rcnrp1;
    static logical goleft, rnkdtd;


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
/*     $Date: 96/12/30 16:59:19 $ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRQPY detects the right rank for upper triangular matrix A. */
/*  The algorithm used here is an version of Pan and Tang's RRQR */
/*  algorithm number 3. */
/*  This algorithm is applied to matrix A until the right rank is */
/*  obtained. If the input ordering of matrix A is not accepted, the */
/*  matrix will be permuted and retriangularized until the rank is */
/*  revealed. */

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

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the m by n matrix A. */
/*          On exit, the upper triangle of the array contains the */
/*          min(m,n) by n upper trapezoidal matrix R; the lower triangle */
/*          array is filled with zeros. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of array A. LDA >= max(1,M). */

/*  C       (input/output) DOUBLE PRECISION array, dimension */
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
/*          JPVT(1:RANK) contains the indices of the columns considered */
/*          linearly independent. */
/*          JPVT(RANK+1:N) contains the indices of the columns considered */
/*          linearly dependent from the previous ones. */

/*  IRCOND  (input) DOUBLE PRECISION */
/*          1/IRCOND specifies an upper bound on the condition number */
/*          of R11. If IRCOND == 0, IRCOND = machine precision is chosen */
/*          as default. IRCOND must be >= 0. */

/*  ORCOND  (output) DOUBLE PRECISION */
/*          1/ORCOND is an estimate for the condition number of R11. */

/*  RANK    (output) INTEGER */
/*          An estimate of the rank offered by this algorithm. */
/*          0 <= RANK <= MIN(M,N). */

/*  SVLUES  (output) DOUBLE PRECISION array, dimension (4) */
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

/*  WORK    (workspace) DOUBLE PRECISION array, dimension ( LWORK ) */

/*  LWORK   (input) INTEGER */
/*          The dimension of array WORK. LWORK >= N+3*MN, where */
/*          MN = min(M,N). */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit. */
/*          < 0: If INFO = -i, the i-th argument had an illegal value */
/*          > 0: Problems in the computation of the rank. */
/*                   1: Exceeded the allowed maximum number of steps. */
/*                   2: Rank not well defined. */
/*          In adition, vector SVLUES tell if rank is not well defined. */
/*          When INFO.NE.0, the contents of ORCOND may be not the right */
/*          one. */


/*  =================================================================== */

/*     .. Parameters .. */

/*     Indices into the `svlues' array. */

/*     .. */
/*     .. */
/*     .. Common Block .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
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
    --jpvt;
    --svlues;
    --work;

    /* Function Body */
    mn = min(*m,*n);
    bsprqr_1.nb = ilaenv(&c__1, "DGEQRF", " ", m, n, &c__0, &c__0, (ftnlen)6,
	     (ftnlen)1);

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
    } else if (*ircond < 0.) {
	*info = -10;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n + mn * 3;
	if (*lwork < max(i__1,i__2)) {
	    *info = -15;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("DTRQPY", &i__1, 6);
	return 0;
    }

/*     Quick return if possible. */

    if (mn == 0) {
	*rank = 0;
	*orcond = 0.;
	svlues[1] = 0.;
	svlues[2] = 0.;
	svlues[3] = 0.;
	svlues[4] = 0.;
	return 0;
    }

/*     Check whether Threshold for condition number was supplied. */
/*     If not, choose machine precision as default for RCOND. */

    if (*ircond > 0.) {
	rcond = *ircond;
    } else {
	rcond = dlamch("Epsilon");
    }

/*     Compute the initial estimate for the rank. */

    dtrrnk(&mn, &a[a_offset], lda, &rcond, rank, &work[1], info);

/*     ************************ */
/*     * First call to xTRQYC * */
/*     ************************ */

/*     Get tighter bounds for the value RANK. */

    dtrqyc(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &jpvt[1],
	    rank, &svlues[1], &rcnr, &rcnrp1, &work[1], info);
    oinfo = 0;
    if (*info != 0) {
	oinfo = *info;
    }

/*     Check if the numerical rank is larger, equal or smaller than */
/*     the contents of RANK. */

    if (rcnr >= rcond && *rank == mn || rcnr >= rcond && rcnrp1 < rcond) {
	rnkdtd = TRUE_;
    } else if (rcnr >= rcond && rcnrp1 >= rcond) {
	rnkdtd = FALSE_;
	goleft = FALSE_;
	++(*rank);
    } else if (rcnr < rcond && rcnrp1 < rcond) {
	if (*rank == 1) {
	    rnkdtd = TRUE_;
	    if ((d__1 = a[a_dim1 + 1], abs(d__1)) == 0.) {
		*rank = 0;
		*orcond = 0.;
		svlues[1] = 0.;
		svlues[2] = 0.;
		svlues[3] = 0.;
		svlues[4] = 0.;
	    } else {
		*rank = 1;
	    }
	} else {
	    rnkdtd = FALSE_;
	    goleft = TRUE_;
	    --(*rank);
	}
    } else {
	rnkdtd = TRUE_;
	*info = 2;
    }

/*     ***************** */
/*     * Start of Loop * */
/*     ***************** */

/*     Loop for the detection of the actual rank. The variable RANK is */
/*     updated until the rank is found. To avoid infinite loops, the */
/*     variable RANK either increases or decreases. */

L10:
    if (! rnkdtd) {

/*        Call to xTRQYC to get tighter bounds for the value RANK. */

	dtrqyc(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &jpvt[1]
		, rank, &svlues[1], &rcnr, &rcnrp1, &work[1], info);
	if (*info != 0) {
	    oinfo = *info;
	}

/*        Check if the numerical rank is larger, equal or smaller than */
/*        the contents of RANK. */

	if (rcnr >= rcond && *rank == mn || rcnr >= rcond && rcnrp1 < rcond) {
	    rnkdtd = TRUE_;
	} else if (rcnr >= rcond && rcnrp1 >= rcond) {
	    if (! goleft) {
		++(*rank);
	    } else {
		rnkdtd = TRUE_;
		*info = 2;
	    }
	} else if (rcnr < rcond && rcnrp1 < rcond) {
	    if (*rank == 1) {
		rnkdtd = TRUE_;
		if ((d__1 = a[a_dim1 + 1], abs(d__1)) == 0.) {
		    *rank = 0;
		    *orcond = 0.;
		    svlues[1] = 0.;
		    svlues[2] = 0.;
		    svlues[3] = 0.;
		    svlues[4] = 0.;
		} else {
		    *rank = 1;
		}
	    } else {
		goleft = TRUE_;
		--(*rank);
	    }
	} else {
	    rnkdtd = TRUE_;
	    *info = 2;
	}

/*        Jump to the beginning of the loop. */

	goto L10;
    }

/*     *************** */
/*     * end of loop * */
/*     *************** */

/*     Give back the obtained value of RCOND and check the value of INFO. */

    *orcond = rcnr;
    if (oinfo != 0) {
	*info = oinfo;
    }

    return 0;

/*     End of DTRQPY */

} /* dtrqpy_ */

