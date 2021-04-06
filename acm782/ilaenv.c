/* ilaenv.f -- translated by f2c (version 20041007).
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
    integer nenvir[10];
} cenvir_;

#define cenvir_1 cenvir_

integer ilaenv(integer *ispec, char *subnam, char *opts, integer *n1,
	integer *n2, integer *n3, integer *n4, ftnlen subnam_len, ftnlen
	opts_len)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char name1[1], name23[2], name46[3], namwrk[6];


/*  -- LAPACK auxiliary routine (preliminary version) -- */
/*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, */
/*     Courant Institute, NAG Ltd., and Rice University */
/*     August 17, 1990 */

/*     ** TEST VERSION ** */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ILAENV returns machine and problem-dependent parameters. */

/*  Arguments */
/*  ========= */

/*  ISPEC   (input) INTEGER */
/*          Specifies what quantity is to be returned (as the function's */
/*          value): */
/*          = 1: The optimum blocksize.  ILAENV(1,...)=1 is a flag that */
/*               no blocking should be done. */
/*          = 2: The minimum blocksize. */
/*          = 3: The "crossover point".  When two versions of a solution */
/*               method are implemented, one of which is faster for */
/*               problems and the other faster for smaller problems, */
/*               this is the largest problem for which the small-problem */
/*               method is preferred. */
/*          = 4: the number of shifts to use.  At present, this is only */
/*               appropriate for the nonsymmetric eigenvalue and */
/*               generalized eigenvalue routines, and -1 will be returned */
/*               if SUBNAM(2:6) is not 'HSEQR' or 'HGEQZ'. */
/*          = 5: The minimum second dimension for blocked updates.  This */
/*               is primarily intended for methods which update */
/*               rectangular blocks using Householder transformations */
/*               with short Householder vectors (e.g., xLAEBC and */
/*               xLAGBC).  If a  k x k  Householder transformation is */
/*               used to update a  k x m  block, then blocking (i.e., */
/*               use of xGEMM) will only be done if k is at least */
/*               ILAENV(2,...) and m is at least ILAENV(5,...) */
/*          = 6: (Used only by the SVD drivers.)  When reducing an m x n */
/*               matrix to bidiagonal form, if the larger dimension is */
/*               less than ILAENV(6,...,m,n,,), then the usual procedure */
/*               is preferred.  If the larger dimension is larger than */
/*               ILAENV(6,...,m,n,,), then it is preferred to first use */
/*               a QR (or LQ) factorization to make it triangular, and */
/*               then reduce it to bidiagonal form. */

/*          = 7: Number of processors. */

/*  SUBNAM  (input) CHARACTER*(*) */
/*          The name of the calling routine, or the routine which is */
/*          expected to use the value returned. */

/*  OPTS    (input) CHARACTER*(*) */
/*          The values of the CHARACTER*1 options passed to the routine */
/*          whose name is in SUBNAM, all run together.  For example, */
/*          a subroutine called with "CALL SGEZZZ('T','Y',N4,'C',A,LDA)" */
/*          would pass 'SGEZZZ' as SUBNAM and 'TYC' as OPTS. */

/*  N1,N2,N3,N4 (input) INTEGER */
/*          The problem dimensions, in the order that they appear in */
/*          the calling sequence.  If there is only one dimension */
/*          (customarily called N), then that value should be passed */
/*          as N1 and N2, N3, and N4 will be ignored, etc. */

/* (ILAENV) (output) INTEGER */
/*          The function value returned will be the value specified */
/*          by ISPEC.  If no reasonable value is available, a negative */
/*          value will be returned, otherwise the returned value will */
/*          be non-negative.  Note that the value returned may be */
/*          unreasonably large for the problem, e.g., a blocksize of */
/*          32 for an 8 x 8 problem, and thus should be restricted */
/*          to whatever range is appropriate. */

/*  Further Details */
/*  ======= ======= */

/*  The calling sequence is intended to match up to the arguments of */
/*  the routine needing the value in a simple and mindless way. */
/*  For example, if SGEZZZ were defined as */

/*       SUBROUTINE SGEZZZ( OPT1, OPT2, N4, OPT3, N1, A, LDA ) */
/*       CHARACTER*1  OPT1, OPT2, OPT3 */
/*       INTEGER      N1, LDA, N4 */
/*       REAL         A(LDA,*) */

/*  then in SGEZZ, the blocksize would be found by a call like: */

/*       NBLOCK = ILAENV( 1, 'SGEZZZ', OPT1//OPT2//OPT3, N4, N1, 0, 0 ) */

/*  It would be further checked and restricted by code like the */
/*  following: */

/*       NBLOCK = MAX( 1, MIN( N4, NBLOCK ) ) */
/*       IF( NBLOCK.EQ.1 ) THEN */
/*  c */
/*  c    unblocked method */
/*  c */
/*    . . . */

/* ======================================================================= */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Arrays in Common .. */
/*     .. */
/*     .. Common blocks .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ISPEC > 7 or ISPEC < 1: Error (return -1) */

    if (*ispec < 1 || *ispec > 7) {
	ret_val = -1;
	return ret_val;
    }

/*     ISPEC=7: Number of processors */

    if (*ispec == 7) {
	ret_val = max(1,cenvir_1.nenvir[6]);
	return ret_val;
    }

/*     ISPEC=6: Jim's crossover. */

    if (*ispec == 6) {
	ret_val = (integer) ((real) max(*n1,*n2) * 1.6f);
	return ret_val;
    }

/*     ISPEC=1 through ISPEC=6: split up name into components */

    s_copy(namwrk, subnam, (ftnlen)6, subnam_len);
    *(unsigned char *)name1 = *(unsigned char *)namwrk;
    s_copy(name23, namwrk + 1, (ftnlen)2, (ftnlen)2);
    s_copy(name46, namwrk + 3, (ftnlen)3, (ftnlen)3);

/*     Test version: just use number from common block */

/* Computing MAX */
    i__1 = 0, i__2 = cenvir_1.nenvir[*ispec - 1];
    ret_val = max(i__1,i__2);
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */

