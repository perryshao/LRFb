// MATLAB MEX gateway for a triangle circumcenter relative to its first vertex.
/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter3d()   Find the circumcenter of a triangle in 3D.         */
/*                                                                           */
/*  The result is returned both in terms of xyz coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
/*                                                                           */
/*****************************************************************************/
#include <mex.h>
#include <matrix.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INEXACT /* Nothing */
/* #define INEXACT volatile */
#define EXACT /* Nothing */

#define REAL double /* float or double */
#define REALPRINT doubleprint
#define REALRAND doublerand
#define NARROWRAND narrowdoublerand
#define UNIFORMRAND uniformdoublerand

#define Absolute(a) ((a) >= 0.0 ? (a) : -(a))
/* #define Absolute(a)  fabs(a) */

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
/*   Split(), and Two_Product() are all implemented as described in the      */
/*   reference.  Each of these macros requires certain variables to be       */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `INEXACT'.                                             */

#define Fast_Two_Sum_Tail(a, b, x, y)                                                              \
    bvirt = x - a;                                                                                 \
    y = b - bvirt

#define Fast_Two_Sum(a, b, x, y)                                                                   \
    x = (REAL)(a + b);                                                                             \
    Fast_Two_Sum_Tail(a, b, x, y)

#define Fast_Two_Diff_Tail(a, b, x, y)                                                             \
    bvirt = a - x;                                                                                 \
    y = bvirt - b

#define Fast_Two_Diff(a, b, x, y)                                                                  \
    x = (REAL)(a - b);                                                                             \
    Fast_Two_Diff_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y)                                                                   \
    bvirt = (REAL)(x - a);                                                                         \
    avirt = x - bvirt;                                                                             \
    bround = b - bvirt;                                                                            \
    around = a - avirt;                                                                            \
    y = around + bround

#define Two_Sum(a, b, x, y)                                                                        \
    x = (REAL)(a + b);                                                                             \
    Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y)                                                                  \
    bvirt = (REAL)(a - x);                                                                         \
    avirt = x + bvirt;                                                                             \
    bround = bvirt - b;                                                                            \
    around = a - avirt;                                                                            \
    y = around + bround

#define Two_Diff(a, b, x, y)                                                                       \
    x = (REAL)(a - b);                                                                             \
    Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo)                                                                         \
    c = (REAL)(splitter * a);                                                                      \
    abig = (REAL)(c - a);                                                                          \
    ahi = c - abig;                                                                                \
    alo = a - ahi

#define Two_Product_Tail(a, b, x, y)                                                               \
    Split(a, ahi, alo);                                                                            \
    Split(b, bhi, blo);                                                                            \
    err1 = x - (ahi * bhi);                                                                        \
    err2 = err1 - (alo * bhi);                                                                     \
    err3 = err2 - (ahi * blo);                                                                     \
    y = (alo * blo) - err3

#define Two_Product(a, b, x, y)                                                                    \
    x = (REAL)(a * b);                                                                             \
    Two_Product_Tail(a, b, x, y)

/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_Presplit(a, b, bhi, blo, x, y)                                                 \
    x = (REAL)(a * b);                                                                             \
    Split(a, ahi, alo);                                                                            \
    err1 = x - (ahi * bhi);                                                                        \
    err2 = err1 - (alo * bhi);                                                                     \
    err3 = err2 - (ahi * blo);                                                                     \
    y = (alo * blo) - err3

/* Two_Product_2Presplit() is Two_Product() where both of the inputs have    */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_2Presplit(a, ahi, alo, b, bhi, blo, x, y)                                      \
    x = (REAL)(a * b);                                                                             \
    err1 = x - (ahi * bhi);                                                                        \
    err2 = err1 - (alo * bhi);                                                                     \
    err3 = err2 - (ahi * blo);                                                                     \
    y = (alo * blo) - err3

/* Square() can be done more quickly than Two_Product().                     */

#define Square_Tail(a, x, y)                                                                       \
    Split(a, ahi, alo);                                                                            \
    err1 = x - (ahi * ahi);                                                                        \
    err3 = err1 - ((ahi + ahi) * alo);                                                             \
    y = (alo * alo) - err3

#define Square(a, x, y)                                                                            \
    x = (REAL)(a * a);                                                                             \
    Square_Tail(a, x, y)

/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

#define Two_One_Sum(a1, a0, b, x2, x1, x0)                                                         \
    Two_Sum(a0, b, _i, x0);                                                                        \
    Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0)                                                        \
    Two_Diff(a0, b, _i, x0);                                                                       \
    Two_Sum(a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0)                                                \
    Two_One_Sum(a1, a0, b0, _j, _0, x0);                                                           \
    Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0)                                               \
    Two_One_Diff(a1, a0, b0, _j, _0, x0);                                                          \
    Two_One_Diff(_j, _0, b1, x3, x2, x1)

#define Four_One_Sum(a3, a2, a1, a0, b, x4, x3, x2, x1, x0)                                        \
    Two_One_Sum(a1, a0, b, _j, x1, x0);                                                            \
    Two_One_Sum(a3, a2, _j, x4, x3, x2)

#define Four_Two_Sum(a3, a2, a1, a0, b1, b0, x5, x4, x3, x2, x1, x0)                               \
    Four_One_Sum(a3, a2, a1, a0, b0, _k, _2, _1, _0, x0);                                          \
    Four_One_Sum(_k, _2, _1, _0, b1, x5, x4, x3, x2, x1)

#define Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0)              \
    Four_Two_Sum(a3, a2, a1, a0, b1, b0, _l, _2, _1, _0, x1, x0);                                  \
    Four_Two_Sum(_l, _2, _1, _0, b4, b3, x7, x6, x5, x4, x3, x2)

#define Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b, x8, x7, x6, x5, x4, x3, x2, x1, x0)       \
    Four_One_Sum(a3, a2, a1, a0, b, _j, x3, x2, x1, x0);                                           \
    Four_One_Sum(a7, a6, a5, a4, _j, x8, x7, x6, x5, x4)

#define Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, x9, x8, x7, x6, x5, x4, x3, x2, x1,  \
                      x0)                                                                          \
    Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, _k, _6, _5, _4, _3, _2, _1, _0, x0);         \
    Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1, x9, x8, x7, x6, x5, x4, x3, x2, x1)

#define Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0, x11, x10, x9, x8, x7, x6,   \
                       x5, x4, x3, x2, x1, x0)                                                     \
    Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, _l, _6, _5, _4, _3, _2, _1, _0, x1, x0); \
    Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2)

/* Macros for multiplying expansions of various fixed lengths.               */

#define Two_One_Product(a1, a0, b, x3, x2, x1, x0)                                                 \
    Split(b, bhi, blo);                                                                            \
    Two_Product_Presplit(a0, b, bhi, blo, _i, x0);                                                 \
    Two_Product_Presplit(a1, b, bhi, blo, _j, _0);                                                 \
    Two_Sum(_i, _0, _k, x1);                                                                       \
    Fast_Two_Sum(_j, _k, x3, x2)

#define Four_One_Product(a3, a2, a1, a0, b, x7, x6, x5, x4, x3, x2, x1, x0)                        \
    Split(b, bhi, blo);                                                                            \
    Two_Product_Presplit(a0, b, bhi, blo, _i, x0);                                                 \
    Two_Product_Presplit(a1, b, bhi, blo, _j, _0);                                                 \
    Two_Sum(_i, _0, _k, x1);                                                                       \
    Fast_Two_Sum(_j, _k, _i, x2);                                                                  \
    Two_Product_Presplit(a2, b, bhi, blo, _j, _0);                                                 \
    Two_Sum(_i, _0, _k, x3);                                                                       \
    Fast_Two_Sum(_j, _k, _i, x4);                                                                  \
    Two_Product_Presplit(a3, b, bhi, blo, _j, _0);                                                 \
    Two_Sum(_i, _0, _k, x5);                                                                       \
    Fast_Two_Sum(_j, _k, x7, x6)

#define Two_Two_Product(a1, a0, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0)                            \
    Split(a0, a0hi, a0lo);                                                                         \
    Split(b0, bhi, blo);                                                                           \
    Two_Product_2Presplit(a0, a0hi, a0lo, b0, bhi, blo, _i, x0);                                   \
    Split(a1, a1hi, a1lo);                                                                         \
    Two_Product_2Presplit(a1, a1hi, a1lo, b0, bhi, blo, _j, _0);                                   \
    Two_Sum(_i, _0, _k, _1);                                                                       \
    Fast_Two_Sum(_j, _k, _l, _2);                                                                  \
    Split(b1, bhi, blo);                                                                           \
    Two_Product_2Presplit(a0, a0hi, a0lo, b1, bhi, blo, _i, _0);                                   \
    Two_Sum(_1, _0, _k, x1);                                                                       \
    Two_Sum(_2, _k, _j, _1);                                                                       \
    Two_Sum(_l, _j, _m, _2);                                                                       \
    Two_Product_2Presplit(a1, a1hi, a1lo, b1, bhi, blo, _j, _0);                                   \
    Two_Sum(_i, _0, _n, _0);                                                                       \
    Two_Sum(_1, _0, _i, x2);                                                                       \
    Two_Sum(_2, _i, _k, _1);                                                                       \
    Two_Sum(_m, _k, _l, _2);                                                                       \
    Two_Sum(_j, _n, _k, _0);                                                                       \
    Two_Sum(_1, _0, _j, x3);                                                                       \
    Two_Sum(_2, _j, _i, _1);                                                                       \
    Two_Sum(_l, _i, _m, _2);                                                                       \
    Two_Sum(_1, _k, _i, x4);                                                                       \
    Two_Sum(_2, _i, _k, x5);                                                                       \
    Two_Sum(_m, _k, x7, x6)

/* An expansion of length two can be squared more quickly than finding the   */
/*   product of two different expansions of length two, and the result is    */
/*   guaranteed to have no more than six (rather than eight) components.     */

#define Two_Square(a1, a0, x5, x4, x3, x2, x1, x0)                                                 \
    Square(a0, _j, x0);                                                                            \
    _0 = a0 + a0;                                                                                  \
    Two_Product(a1, _0, _k, _1);                                                                   \
    Two_One_Sum(_k, _1, _j, _l, _2, x1);                                                           \
    Square(a1, _j, _1);                                                                            \
    Two_Two_Sum(_j, _1, _l, _2, x5, x4, x3, x2);
REAL splitter; /* = 2^ceiling(p / 2) + 1.  Used to split floats in half. */
REAL epsilon;  /* = 2^(-p).  Used to estimate roundoff errors. */
/* A set of coefficients used to calculate maximum roundoff errors.          */
REAL resulterrbound;
REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
REAL o3derrboundA, o3derrboundB, o3derrboundC;
REAL iccerrboundA, iccerrboundB, iccerrboundC;
REAL isperrboundA, isperrboundB, isperrboundC;

/*****************************************************************************/
/*                                                                           */
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*                                                                           */
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-point error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-point numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  I imagine that a highly optimizing compiler might be too smart for its   */
/*  own good, and somehow cause this routine to fail, if it pretends that    */
/*  floating-point arithmetic is too much like real arithmetic.              */
/*                                                                           */
/*  Don't change this routine unless you fully understand it.                */
/*                                                                           */
/*****************************************************************************/

void exactinit()
{
    REAL half;
    REAL check, lastcheck;
    int every_other;

    every_other = 1;
    half = 0.5;
    epsilon = 1.0;
    splitter = 1.0;
    check = 1.0;
    /* Repeatedly divide `epsilon' by two until it is too small to add to    */
    /*   one without causing roundoff.  (Also check if the sum is equal to   */
    /*   the previous sum, for machines that round up instead of using exact */
    /*   rounding.  Not that this library will work on such machines anyway. */
    do
    {
        lastcheck = check;
        epsilon *= half;
        if (every_other)
        {
            splitter *= 2.0;
        }
        every_other = !every_other;
        check = 1.0 + epsilon;
    } while ((check != 1.0) && (check != lastcheck));
    splitter += 1.0;

    /* Error bounds for orientation and incircle tests. */
    resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
    ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
    ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
    ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
    o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
    o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
    o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
    iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
    iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
    iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
    isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
    isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
    isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;
}

REAL estimate(int elen, REAL *e)
{
    REAL Q;
    int eindex;

    Q = e[0];
    for (eindex = 1; eindex < elen; eindex++)
    {
        Q += e[eindex];
    }
    return Q;
}

REAL orient2dadapt(REAL *pa, REAL *pb, REAL *pc, REAL detsum)
{
    INEXACT REAL acx, acy, bcx, bcy;
    REAL acxtail, acytail, bcxtail, bcytail;
    INEXACT REAL detleft, detright;
    REAL detlefttail, detrighttail;
    REAL det, errbound;
    REAL B[4], C1[8], C2[12], D[16];
    INEXACT REAL B3;
    int C1length, C2length, Dlength;
    REAL u[4];
    INEXACT REAL u3;
    INEXACT REAL s1, t1;
    REAL s0, t0;

    INEXACT REAL bvirt;
    REAL avirt, bround, around;
    INEXACT REAL c;
    INEXACT REAL abig;
    REAL ahi, alo, bhi, blo;
    REAL err1, err2, err3;
    INEXACT REAL _i, _j;
    REAL _0;

    acx = (REAL)(pa[0] - pc[0]);
    bcx = (REAL)(pb[0] - pc[0]);
    acy = (REAL)(pa[1] - pc[1]);
    bcy = (REAL)(pb[1] - pc[1]);

    Two_Product(acx, bcy, detleft, detlefttail);
    Two_Product(acy, bcx, detright, detrighttail);

    Two_Two_Diff(detleft, detlefttail, detright, detrighttail, B3, B[2], B[1], B[0]);
    B[3] = B3;

    det = estimate(4, B);
    errbound = ccwerrboundB * detsum;
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
    Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
    Two_Diff_Tail(pa[1], pc[1], acy, acytail);
    Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

    if ((acxtail == 0.0) && (acytail == 0.0) && (bcxtail == 0.0) && (bcytail == 0.0))
    {
        return det;
    }
}

REAL orient2d(REAL *pa, REAL *pb, REAL *pc)
{

    REAL detleft, detright, det;
    REAL detsum, errbound;

    detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    det = detleft - detright;
    exactinit();
    if (detleft > 0.0)
    {
        if (detright <= 0.0)
        {
            return det;
        }
        else
        {
            detsum = detleft + detright;
        }
    }
    else if (detleft < 0.0)
    {
        if (detright >= 0.0)
        {
            return det;
        }
        else
        {
            detsum = -detleft - detright;
        }
    }
    else
    {
        return det;
    }

    errbound = ccwerrboundA * detsum;
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    return orient2dadapt(pa, pb, pc, detsum);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3 || nlhs < 1 || nlhs > 3)
        mexErrMsgIdAndTxt("LRFb:circumcenterArity", "Expected three points and one to three outputs.");
    for (int i = 0; i < 3; ++i)
    {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) ||
            mxGetNumberOfElements(prhs[i]) != 3)
            mexErrMsgIdAndTxt("LRFb:circumcenterInput", "Each point must be a full real double 3-vector.");
    }
    double *a, *b, *c;
    double pa[2], pb[2], pc[2];
    double *circumcenter;
    double *xi;
    double *eta;
    int m, n;

    double xba, yba, zba, xca, yca, zca;
    double balength, calength;
    double xcrossbc, ycrossbc, zcrossbc;
    double denominator;
    double xcirca, ycirca, zcirca;

    /* parse input arguments*/
    a = mxGetPr(prhs[0]); //a
    b = mxGetPr(prhs[1]); //b
    c = mxGetPr(prhs[2]); //c
    m = mxGetM(prhs[0]);  //rows for a
    n = mxGetN(prhs[0]);  //colums for a

    /* create output arguments*/
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); //circumcenter
    circumcenter = mxGetPr(plhs[0]);
    // Optional coordinates still participate in the shared formula. Keep
    // unrequested values local instead of writing beyond the output array.
    double local_xi = 0.0, local_eta = 0.0;
    xi = &local_xi;
    eta = &local_eta;
    if (nlhs >= 2)
    {
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        xi = mxGetPr(plhs[1]);
    }
    if (nlhs >= 3)
    {
        plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
        eta = mxGetPr(plhs[2]);
    }

    /* Use coordinates relative to point `a' of the triangle. */
    xba = b[0] - a[0];
    yba = b[1] - a[1];
    zba = b[2] - a[2];
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    zca = c[2] - a[2];
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba + zba * zba;
    calength = xca * xca + yca * yca + zca * zca;

    /* Cross product of these edges. */
#ifdef EXACT
    /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
    /*   to ensure a correctly signed (and reasonably accurate) result, */
    /*   avoiding any possibility of division by zero.                  */
    pa[0] = a[1];
    pa[1] = a[2];
    pb[0] = b[1];
    pb[1] = b[2];
    pc[0] = c[1];
    pc[1] = c[2];
    xcrossbc = orient2d(pb, pc, pa);
    pa[0] = a[2];
    pa[1] = a[0];
    pb[0] = b[2];
    pb[1] = b[0];
    pc[0] = c[2];
    pc[1] = c[0];
    ycrossbc = orient2d(pb, pc, pa);
    pa[0] = a[0];
    pa[1] = a[1];
    pb[0] = b[0];
    pb[1] = b[1];
    pc[0] = c[0];
    pc[1] = c[1];
    zcrossbc = orient2d(pb, pc, pa);
#else
    /* Take your chances with floating-point roundoff. */
    xcrossbc = yba * zca - yca * zba;
    ycrossbc = zba * xca - zca * xba;
    zcrossbc = xba * yca - xca * yba;
#endif

    /* Calculate the denominator of the formulae. */
    denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc + zcrossbc * zcrossbc);

    /* Calculate offset (from `a') of circumcenter. */
    xcirca = ((balength * yca - calength * yba) * zcrossbc -
              (balength * zca - calength * zba) * ycrossbc) *
             denominator;
    ycirca = ((balength * zca - calength * zba) * xcrossbc -
              (balength * xca - calength * xba) * zcrossbc) *
             denominator;
    zcirca = ((balength * xca - calength * xba) * ycrossbc -
              (balength * yca - calength * yba) * xcrossbc) *
             denominator;
    circumcenter[0] = xcirca;
    circumcenter[1] = ycirca;
    circumcenter[2] = zcirca;

    if (xi != (double *)NULL)
    {
        /* To interpolate a linear function at the circumcenter, define a     */
        /*   coordinate system with a xi-axis directed from `a' to `b' and    */
        /*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
        /*   are computed by Cramer's Rule for solving systems of linear      */
        /*   equations.                                                       */

        /* There are three ways to do this calculation - using xcrossbc, */
        /*   ycrossbc, or zcrossbc.  Choose whichever has the largest    */
        /*   magnitude, to improve stability and avoid division by zero. */
        if (((xcrossbc >= ycrossbc) ^ (-xcrossbc > ycrossbc)) &&
            ((xcrossbc >= zcrossbc) ^ (-xcrossbc > zcrossbc)))
        {
            *xi = (ycirca * zca - zcirca * yca) / xcrossbc;
            *eta = (zcirca * yba - ycirca * zba) / xcrossbc;
        }
        else if ((ycrossbc >= zcrossbc) ^ (-ycrossbc > zcrossbc))
        {
            *xi = (zcirca * xca - xcirca * zca) / ycrossbc;
            *eta = (xcirca * zba - zcirca * xba) / ycrossbc;
        }
        else
        {
            *xi = (xcirca * yca - ycirca * xca) / zcrossbc;
            *eta = (ycirca * xba - xcirca * yba) / zcrossbc;
        }
    }
}
