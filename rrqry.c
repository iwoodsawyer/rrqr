/*
 * Rank Revealing QR Factorization
 *
 * Based on Pan&Tang's algorithm number 3.
 *
 * R = rrqry(A)
 * [Q,R,p] = rrqry(A)
 * [Q,R,p,r] = rrqry(A)
 * [Q,R,p,r] = rrqry(A,0)
 * [Q,R,p,r] = rrqry(A,job)
 * [Q,R,p,r] = rrqry(A,job,0)
 *
 * compile command:
 * mex -O rrqry.c libmwlapack.lib rrqr.lib
 * or
 * mex -O rrqry.c libmwblas.lib libmwlapack.lib rrqr.lib (>= R2007B)
 *
 * calls the DGEQPY/SGEQPY/CGEQPY/ZGEQPY functions from RRQR library
 */

#include "mex.h"
#include "matrix.h"
#include "rrqr.h"

void rrqry_double(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    ptrdiff_t *Pp, ipspec, lwork, rank, job, info = 1;
    double p, orcond, ircond = 0;
    double *Qpr, *Rpr, *Ipr, *Ppr, *Ap, *Cp, *slvalues, *pwork, *rwork;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    double *Qpi, *Rpi, *Ipi;
    #endif
    char name[] = "DGEQRF", opts[] = " ";
    size_t m, n, n2, k, ldc, ldr, max_mn, min_mn, nb;
    size_t element_size = sizeof(double), econ = 0, cplx = 0, dc = 1;
    mwIndex i, j, limit;
    mxClassID classid = mxDOUBLE_CLASS;
    mxComplexity cplxflag = mxREAL;
    
    /* check complex */
    if (mxIsComplex(prhs[0])) {
        cplxflag = mxCOMPLEX;
        cplx = 1;
        dc = 2;
    }
    
    /* determine job and economy size */
    if (nrhs == 3) {
        job = (ptrdiff_t)mxGetScalar(prhs[1]);
        if (mxGetScalar(prhs[2]) == 0) {
            econ = 1;
        }
    }
    else if (nrhs == 2) {
        job = (ptrdiff_t)mxGetScalar(prhs[1]);
        if (job == 0) {
            econ = 1;
            job = 3;
        }
    }
    else if (nrhs == 1) {
        if (nlhs == 1) {
            job = 1;
        }
        else {
            job = 3;
        }
    }
    
    /* get matrix data */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (m == 0 || n == 0) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        }
        else {
            if (nlhs >= 3) {
                plhs[2] = mxCreateNumericMatrix(n,1,classid,cplxflag);
                if (n != 0) {
                    Ppr = mxGetData(plhs[2]);
                    p = 1;
                    for (i=0; i<n; i++) {
                        Ppr[dc*i] = p;
                        p = p + 1;
                    }
                }
            }
            if (nlhs == 4) {
                plhs[3] = mxCreateNumericMatrix(0,0,classid,cplxflag);
            }
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
            if (m != 0) {
                Qpr = mxGetData(plhs[0]);
                for (i=0; i<m; i++) {
                    Qpr[i*dc*m+dc*i] = 1;
                }
            }
        }
        return;
    }
    
    /* minimal size */
    #if defined(_MSC_VER)
    min_mn = min(m,n);
    #else
    min_mn = m < n ? m : n;
    #endif
    
    /* maximal size */
    #if defined(_MSC_VER)
    max_mn = max(m,n);
    #else
    max_mn = m > n ? m : n;
    #endif
    
    /* allocate slvalues */
    slvalues = mxMalloc(4*element_size);
    
    /* allocate A matrix */
    Ap = mxMalloc(dc*m*n*element_size);
    
    /* allocate C matrix */
    switch (job) {
        case 1:
            k = 0;
            ldc = 0;
            break;
        case 2:
            k = m;
            ldc = m;
            Cp = mxCalloc(dc*ldc*k,element_size);
            for (i=0; i<m; i++) {
                Cp[i*dc*m+dc*i] = 1;
            }
            break;
        case 3:
            k = m;
            ldc = k;
            Cp = mxCalloc(dc*ldc*m,element_size);
            for (i=0; i<m; i++) {
                Cp[i*dc*m+dc*i] = 1;
            }
            break;
    }
    
    /* allocate P matrix */
    Pp = mxMalloc(n*sizeof(ptrdiff_t));
    
    /* copy input to A matrix */
    Ipr = mxGetData(prhs[0]);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        Ipi = mxGetImagData(prhs[0]);
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                Ap[j*2*m+2*i] = Ipr[j*m+i];
                Ap[j*2*m+2*i+1] = Ipi[j*m+i];
            }
        }
    }
    else {
    #endif
        memcpy(Ap,Ipr,dc*m*n*element_size);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif
    
    /* determine blocksize */
    if (cplx) {
        name[0] = 'Z';
    }
    ipspec = 1;
    nb = ilaenv(&ipspec, name, opts, &m, &n, &m, &n, 6, 1);
    
    /* allocate workspace */
    if (job == 1) {
        if (nb < 3) {
            lwork = 2*min_mn + 3*n;
        }
        else {
            lwork = 2*min_mn + n*nb;
        }
    }
    else {
        if (nb < 3) {
            lwork = 2*min_mn + 2*n + max_mn;
        }
        else {
            lwork = 2*min_mn + nb*nb + nb*max_mn;
        }
    }
    pwork = mxMalloc(dc*lwork*element_size);
    if (cplx) {
        rwork = mxMalloc(dc*n*element_size);
    }
    
    /* calls the DGEQPY function */
    if (cplx) {
        zgeqpy(&job, &m, &n, &k, Ap, &m, Cp, &ldc, Pp, &ircond, &orcond, &rank, slvalues, pwork, &lwork, rwork, &info);
    }
    else {
        dgeqpy(&job, &m, &n, &k, Ap, &m, Cp, &ldc, Pp, &ircond, &orcond, &rank, slvalues, pwork, &lwork, &info);
    }
    mxFree(pwork);
    mxFree(slvalues);
    if (cplx) {
        mxFree(rwork);
    }
    
    /* check info for errors */
    if (info != 0) {
        mxFree(Ap);
        if (info == 1) {
            mexErrMsgTxt("Exceeded the allowed maximum number of steps.");
        }
        else if (info == 1) {
            mexErrMsgTxt("Rank not well defined.");
        }
        else {
            if (cplx) {
                mexErrMsgTxt("ZGEQPY not succesful");
            }
            else {
                mexErrMsgTxt("DGEQPY not succesful");
            }
        }
    }
    
    /* extract lower triangular part of Ap */
    if ((econ == 1) && m > n) {
        ldr = n;
    }
    else {
        ldr = m;
    }
    if (nlhs == 1) {
        plhs[0] = mxCreateNumericMatrix(ldr,n,classid,cplxflag);
        Rpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Rpi = mxGetImagData(plhs[0]);
        }
        #endif
    }
    else {
        plhs[1] = mxCreateNumericMatrix(ldr,n,classid,cplxflag);
        Rpr = mxGetData(plhs[1]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Rpi = mxGetImagData(plhs[1]);
        }
        #endif
    }
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        for (j=0; j<n; j++) {
            limit = j < min_mn-1 ? j+1 : min_mn;
            for (i=0; i<limit; i++) {
                Rpr[j*ldr+i] = Ap[j*2*m+2*i];
                Rpi[j*ldr+i] = Ap[j*2*m+2*i+1];
            }
        }
    }
    else {
    #endif
        for (j=0; j<n; j++) {
            limit = j < min_mn-1 ? j+1 : min_mn;
            for (i=0; i<dc*limit; i++) {
                Rpr[j*dc*ldr+i] = Ap[j*dc*m+i];
            }
        }
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif
    mxFree(Ap);
    
    /* extract the orthogonal matrix */
    if (nlhs >= 2) {
        n2 = (econ == 1) ? min_mn : m;
        
        /* allocate Q matrix */
        switch (job) {
            case 1:
                plhs[0] = mxCreateNumericMatrix(m,n2,classid,cplxflag);
                break;
            case 2:
                plhs[0] = mxCreateNumericMatrix(n2,m,classid,cplxflag);
                Qpr = mxGetData(plhs[0]);
                #if !(MX_HAS_INTERLEAVED_COMPLEX)
                if (cplx) {
                    Qpi = mxGetImagData(plhs[0]);
                    
                    /* copy Cp to Qp */
                    for (j=0; j<m; j++) {
                        for (i=0; i<n2; i++) {
                            Qpr[j*n2+i] = Cp[j*2*m+2*i];
                            Qpi[j*n2+i] = Cp[j*2*m+2*i+1];
                        }
                    }
                }
                else {
                #endif
                    /* copy Cp to Qp */
                    for (j=0; j<m; j++) {
                        for (i=0; i<dc*n2; i++) {
                            Qpr[j*dc*n2+i] = Cp[j*dc*m+i];
                        }
                    }
                #if !(MX_HAS_INTERLEAVED_COMPLEX)
                }
                #endif
                break;
            case 3:
                plhs[0] = mxCreateNumericMatrix(m,n2,classid,cplxflag);
                Qpr = mxGetData(plhs[0]);
                #if !(MX_HAS_INTERLEAVED_COMPLEX)
                if (cplx) {
                    Qpi = mxGetImagData(plhs[0]);

                    /* copy Cp to Qp */
                    for (i=0; i<m; i++) {
                        for (j=0; j<n2; j++) {
                            Qpr[j*m+i] = Cp[j*2*m+2*i];
                            Qpi[j*m+i] = Cp[j*2*m+2*i+1];
                        }
                    }
                }
                else {
                 #endif
                    /* copy Cp to Qp */
                    memcpy(Qpr,Cp,dc*m*n2*element_size);
                #if !(MX_HAS_INTERLEAVED_COMPLEX)
                }
                #endif
                break;
        }
    }
    mxFree(Cp);
    
    /* extract the permutation vector */
    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleMatrix(n,1,mxREAL);
        Ppr = mxGetData(plhs[2]);      
        for (i=0; i<n; i++) {
            Ppr[i] = (double)Pp[i];
        }
    }
    mxFree(Pp);
    
    /* extract the rank */
    if (nlhs >= 4) {
        plhs[3] = mxCreateDoubleScalar((double)rank);
    }
}


void rrqry_single(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    ptrdiff_t *Pp, ipspec, lwork, rank, job, info = 1;
    float p, orcond, ircond = 0;
    float *Qpr, *Rpr, *Ipr, *Ap, *Cp, *slvalues, *pwork, *rwork;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    float *Qpi, *Rpi, *Ipi;
    #endif
    double *Ppr;
    char name[] = "SGEQRF", opts[] = " ";
    size_t m, n, n2, k, ldc, ldr, max_mn, min_mn, nb;
    size_t element_size = sizeof(float), econ = 0, cplx = 0, dc = 1;
    mwIndex i, j, limit;
    mxClassID classid = mxSINGLE_CLASS;
    mxComplexity cplxflag = mxREAL;
    
    /* check complex */
    if (mxIsComplex(prhs[0])) {
        cplxflag = mxCOMPLEX;
        cplx = 1;
        dc = 2;
    }
    
    /* determine job and economy size */
    if (nrhs == 3) {
        job = (ptrdiff_t)mxGetScalar(prhs[1]);
        if (mxGetScalar(prhs[2]) == 0) {
            econ = 1;
        }
    }
    else if (nrhs == 2) {
        job = (ptrdiff_t)mxGetScalar(prhs[1]);
        if (job == 0) {
            econ = 1;
            job = 3;
        }
    }
    else if (nrhs == 1) {
        if (nlhs == 1) {
            job = 1;
        }
        else {
            job = 3;
        }
    }
    
    /* get matrix data */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (m == 0 || n == 0) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        }
        else {
            if (nlhs >= 3) {
                plhs[2] = mxCreateNumericMatrix(n,1,classid,cplxflag);
                if (n != 0) {
                    Ppr = mxGetData(plhs[2]);
                    p = 1;
                    for (i=0; i<n; i++) {
                        Ppr[dc*i] = p;
                        p = p + 1;
                    }
                }
            }
            if (nlhs == 4) {
                plhs[3] = mxCreateNumericMatrix(0,0,classid,cplxflag);
            }
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
            if (m != 0) {
                Qpr = mxGetData(plhs[0]);
                for (i=0; i<m; i++) {
                    Qpr[i*dc*m+dc*i] = 1;
                }
            }
        }
        return;
    }
    
    /* minimal size */
    #if defined(_MSC_VER)
    min_mn = min(m,n);
    #else
    min_mn = m < n ? m : n;
    #endif
    
    /* maximal size */
    #if defined(_MSC_VER)
    max_mn = max(m,n);
    #else
    max_mn = m > n ? m : n;
    #endif
    
    /* allocate slvalues */
    slvalues = mxMalloc(4*element_size);
    
    /* allocate A matrix */
    Ap = mxMalloc(dc*m*n*element_size);
    
    /* allocate C matrix */
    switch (job) {
        case 1:
            k = 0;
            ldc = 0;
            break;
        case 2:
            k = m;
            ldc = m;
            Cp = mxCalloc(dc*ldc*k,element_size);
            for (i=0; i<m; i++) {
                Cp[i*dc*m+dc*i] = 1;
            }
            break;
        case 3:
            k = m;
            ldc = k;
            Cp = mxCalloc(dc*ldc*m,element_size);
            for (i=0; i<m; i++) {
                Cp[i*dc*m+dc*i] = 1;
            }
            break;
    }
    
    /* allocate P matrix */
    Pp = mxMalloc(n*sizeof(ptrdiff_t));
    
    /* copy input to A matrix */
    Ipr = mxGetData(prhs[0]);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        Ipi = mxGetImagData(prhs[0]);
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                Ap[j*2*m+2*i] = Ipr[j*m+i];
                Ap[j*2*m+2*i+1] = Ipi[j*m+i];
            }
        }
    }
    else {
    #endif
        memcpy(Ap,Ipr,dc*m*n*element_size);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif
    
    /* determine blocksize */
    if (cplx) {
        name[0] = 'C';
    }
    ipspec = 1;
    nb = ilaenv(&ipspec, name, opts, &m, &n, &m, &n, 6, 1);
    
    /* allocate workspace */
    if (job == 1) {
        if (nb < 3) {
            lwork = 2*min_mn + 3*n;
        }
        else {
            lwork = 2*min_mn + n*nb;
        }
    }
    else {
        if (nb < 3) {
            lwork = 2*min_mn + 2*n + max_mn;
        }
        else {
            lwork = 2*min_mn + nb*nb + nb*max_mn;
        }
    }
    pwork = mxMalloc(dc*lwork*element_size);
    if (cplx) {
        rwork = mxMalloc(dc*n*element_size);
    }
    
    /* calls the DGEQPY function */
    if (cplx) {
        cgeqpy(&job, &m, &n, &k, Ap, &m, Cp, &ldc, Pp, &ircond, &orcond, &rank, slvalues, pwork, &lwork, rwork, &info);
    }
    else {
        sgeqpy(&job, &m, &n, &k, Ap, &m, Cp, &ldc, Pp, &ircond, &orcond, &rank, slvalues, pwork, &lwork, &info);
    }
    mxFree(pwork);
    mxFree(slvalues);
    if (cplx) {
        mxFree(rwork);
    }
    
    /* check info for errors */
    if (info != 0) {
        mxFree(Ap);
        if (info == 1) {
            mexErrMsgTxt("Exceeded the allowed maximum number of steps.");
        }
        else if (info == 1) {
            mexErrMsgTxt("Rank not well defined.");
        }
        else {
            if (cplx) {
                mexErrMsgTxt("CGEQPY not succesful");
            }
            else {
                mexErrMsgTxt("SGEQPY not succesful");
            }
        }
    }
    
    /* extract lower triangular part of Ap */
    if ((econ == 1) && m > n) {
        ldr = n;
    }
    else {
        ldr = m;
    }
    if (nlhs == 1) {
        plhs[0] = mxCreateNumericMatrix(ldr,n,classid,cplxflag);
        Rpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Rpi = mxGetImagData(plhs[0]);
        }
        #endif
    }
    else {
        plhs[1] = mxCreateNumericMatrix(ldr,n,classid,cplxflag);
        Rpr = mxGetData(plhs[1]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Rpi = mxGetImagData(plhs[1]);
        }
        #endif
    }
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        for (j=0; j<n; j++) {
            limit = j < min_mn-1 ? j+1 : min_mn;
            for (i=0; i<limit; i++) {
                Rpr[j*ldr+i] = Ap[j*2*m+2*i];
                Rpi[j*ldr+i] = Ap[j*2*m+2*i+1];
            }
        }
    }
    else {
    #endif
        for (j=0; j<n; j++) {
            limit = j < min_mn-1 ? j+1 : min_mn;
            for (i=0; i<dc*limit; i++) {
                Rpr[j*dc*ldr+i] = Ap[j*dc*m+i];
            }
        }
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif
    mxFree(Ap);
    
    /* extract the orthogonal matrix */
    if (nlhs >= 2) {
        n2 = (econ == 1) ? min_mn : m;
        
        /* allocate Q matrix */
        switch (job) {
            case 1:
                plhs[0] = mxCreateNumericMatrix(m,n2,classid,cplxflag);
                break;
            case 2:
                plhs[0] = mxCreateNumericMatrix(n2,m,classid,cplxflag);
                Qpr = mxGetData(plhs[0]);
                #if !(MX_HAS_INTERLEAVED_COMPLEX)
                if (cplx) {
                    Qpi = mxGetImagData(plhs[0]);
                    
                    /* copy Cp to Qp */
                    for (j=0; j<m; j++) {
                        for (i=0; i<n2; i++) {
                            Qpr[j*n2+i] = Cp[j*2*m+2*i];
                            Qpi[j*n2+i] = Cp[j*2*m+2*i+1];
                        }
                    }
                }
                else {
                #endif
                    /* copy Cp to Qp */
                    for (j=0; j<m; j++) {
                        for (i=0; i<dc*n2; i++) {
                            Qpr[j*dc*n2+i] = Cp[j*dc*m+i];
                        }
                    }
                #if !(MX_HAS_INTERLEAVED_COMPLEX)
                }
                #endif
                break;
            case 3:
                plhs[0] = mxCreateNumericMatrix(m,n2,classid,cplxflag);
                Qpr = mxGetData(plhs[0]);
                #if !(MX_HAS_INTERLEAVED_COMPLEX)
                if (cplx) {
                    Qpi = mxGetImagData(plhs[0]);

                    /* copy Cp to Qp */
                    for (i=0; i<m; i++) {
                        for (j=0; j<n2; j++) {
                            Qpr[j*m+i] = Cp[j*2*m+2*i];
                            Qpi[j*m+i] = Cp[j*2*m+2*i+1];
                        }
                    }
                }
                else {
                 #endif
                    /* copy Cp to Qp */
                    memcpy(Qpr,Cp,dc*m*n2*element_size);
                #if !(MX_HAS_INTERLEAVED_COMPLEX)
                }
                #endif
                break;
        }
    }
    mxFree(Cp);
    
    /* extract the permutation vector */
    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleMatrix(n,1,mxREAL);
        Ppr = mxGetData(plhs[2]);      
        for (i=0; i<n; i++) {
            Ppr[i] = (double)Pp[i];
        }
    }
    mxFree(Pp);
    
    /* extract the rank */
    if (nlhs >= 4) {
        plhs[3] = mxCreateDoubleScalar((double)rank);
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if (nrhs != 1 && nrhs != 2 && nrhs != 3) {
        mexErrMsgTxt("RRQRY requires one, two or three arguments.");
    }
    if (nlhs > 4) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (!mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgTxt( "Input must be a full matrix." );
    }
    if (mxIsDouble(prhs[0])) {
        rrqry_double(nlhs, plhs, nrhs, prhs);
    }
    else if (mxIsSingle(prhs[0])) {
        rrqry_single(nlhs, plhs, nrhs, prhs);
    }
    else {
        mexErrMsgTxt( "Class is not supported." );
    }
}
