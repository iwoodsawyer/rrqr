/* clasmx.f -- translated by f2c (version 20041007).
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

static doublereal c_b2 = .33333333333333331;

doublereal clasmx(integer *i__)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);


    d__1 = (doublereal) ((real) (*i__));
    ret_val = pow_dd(&d__1, &c_b2);
    return ret_val;
} /* slasmx_ */

