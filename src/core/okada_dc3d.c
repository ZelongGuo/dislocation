/* ../dc3d.c -- translated by f2c (version 20230428).
 *
 * **Finite Rectangular Source**
 *
 * The C version of Fortran codes of Okada dislocation (Okada, 1992).
 * Zelong Guo, 29.01.2024
 * @ Potsdam, Germany
 *
 */

/* #include "f2c.h" */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Common Block Declarations */
typedef double doublereal;
typedef long int integer;

static union {
    struct {
        doublereal dummy[5], sd, cd;
    } _1;
    struct {
        doublereal alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d,
            c2d;
    } _2;
} c0_;

#define c0_1 (c0_._1)
#define c0_2 (c0_._2)

struct {
    doublereal xi2, et2, q2, r__, r2, r3, r5, y, d__, tt, alx, ale, x11, y11,
        x32, y32, ey, ez, fy, fz, gy, gz, hy, hz;
} c2_;

#define c2_1 c2_

/* Subroutine */
int dc3d_(doublereal *alpha, doublereal *x, doublereal *y, doublereal *z__,
          doublereal *depth, doublereal *dip, doublereal *al1, doublereal *al2,
          doublereal *aw1, doublereal *aw2, doublereal *disl1,
          doublereal *disl2, doublereal *disl3, doublereal *ux, doublereal *uy,
          doublereal *uz, doublereal *uxx, doublereal *uyx, doublereal *uzx,
          doublereal *uxy, doublereal *uyy, doublereal *uzy, doublereal *uxz,
          doublereal *uyz, doublereal *uzz, int *iret) {
    /* Initialized data */
    static doublereal f0 = 0.;
    static doublereal eps = 1e-6;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;
    static doublereal p, q, u[12], r12, r21, r22, et[2], du[12];

    /* Subroutine */
    extern int ua_(doublereal *, doublereal *, doublereal *, doublereal *,
                   doublereal *, doublereal *, doublereal *),
        ub_(doublereal *, doublereal *, doublereal *, doublereal *,
            doublereal *, doublereal *, doublereal *),
        uc_(doublereal *, doublereal *, doublereal *, doublereal *,
            doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal xi[2], zz, dd1, dd2, dd3, dua[12], dub[12], duc[12];
    static integer ket[2], kxi[2];
    static doublereal ddip;
    /* Subroutine */
    extern int dccon0_(doublereal *, doublereal *),
        dccon2_(doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, integer *, integer *);
    static doublereal aalpha;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
    /*<        >*/

    /* ******************************************************************** */
    /* *****                                                          ***** */
    /* *****    DISPLACEMENT AND STRAIN AT DEPTH                      ***** */
    /* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
    /* *****              CODED BY  Y.OKADA ... SEP.1991              ***** */
    /* *****              REVISED ... NOV.1991, APR.1992, MAY.1993,   ***** */
    /* *****                          JUL.1993, MAY.2002              ***** */
    /* ******************************************************************** */

    /* ***** INPUT */
    /* *****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU) */
    /* *****   X,Y,Z : COORDINATE OF OBSERVING POINT */
    /* *****   DEPTH : DEPTH OF REFERENCE POINT */
    /* *****   DIP   : DIP-ANGLE (DEGREE) */
    /* *****   AL1,AL2   : FAULT LENGTH RANGE */
    /* *****   AW1,AW2   : FAULT WIDTH RANGE */
    /* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS */

    /* ***** OUTPUT */
    /* *****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL) */
    /* *****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) / */
    /* *****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )
     */
    /* *****   UXZ,UYZ,UZZ : Z-DERIVATIVE */
    /* *****   IRET        : RETURN CODE */
    /* *****               :   =0....NORMAL */
    /* *****               :   =1....SINGULAR */
    /* *****               :   =2....POSITIVE Z WAS GIVEN */

    /*<       COMMON /C0/DUMMY(5),SD,CD >*/
    /*<       DIMENSION  XI(2),ET(2),KXI(2),KET(2) >*/
    /*<       DIMENSION  U(12),DU(12),DUA(12),DUB(12),DUC(12) >*/
    /*<       DATA  F0,EPS/ 0.D0, 1.D-6 / >*/
    /* ----- */
    /*<       IRET=0 >*/
    *iret = 0;
    /*<       IF(Z.GT.0.) THEN >*/
    if (*z__ > 0.f) {
        /*<         IRET=2 >*/
        *iret = 2;
        /*<         GO TO 99 >*/
        goto L99;
        // printf("\n%s\n", STARS);
        // fprintf(stderr, "Error, positive z is specified!");
        // printf("\n%s\n", STARS);
        // exit(EXIT_FAILURE);
        // puts("**Note**: positive Z given!!!");
        // printf("\n%s\n", STARS);
        /*<       ENDIF >*/
    }
    /* ----- */
    /*<       DO 111 I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<         U  (I)=F0 >*/
        u[i__ - 1] = f0;
        /*<         DUA(I)=F0 >*/
        dua[i__ - 1] = f0;
        /*<         DUB(I)=F0 >*/
        dub[i__ - 1] = f0;
        /*<         DUC(I)=F0 >*/
        duc[i__ - 1] = f0;
        /*<   111 CONTINUE >*/
        /* L111: */
    }
    /*<       AALPHA=ALPHA >*/
    aalpha = *alpha;
    /*<       DDIP=DIP >*/
    ddip = *dip;
    /*<       CALL DCCON0(AALPHA,DDIP) >*/
    dccon0_(&aalpha, &ddip);
    /* ----- */
    /*<       ZZ=Z >*/
    zz = *z__;
    /*<       DD1=DISL1 >*/
    dd1 = *disl1;
    /*<       DD2=DISL2 >*/
    dd2 = *disl2;
    /*<       DD3=DISL3 >*/
    dd3 = *disl3;
    /*<       XI(1)=X-AL1 >*/
    xi[0] = *x - *al1;
    /*<       XI(2)=X-AL2 >*/
    xi[1] = *x - *al2;
    /*<       IF(DABS(XI(1)).LT.EPS) XI(1)=F0 >*/
    if (fabs(xi[0]) < eps) {
        xi[0] = f0;
    }
    /*<       IF(DABS(XI(2)).LT.EPS) XI(2)=F0 >*/
    if (fabs(xi[1]) < eps) {
        xi[1] = f0;
    }
    /* ====================================== */
    /* =====  REAL-SOURCE CONTRIBUTION  ===== */
    /* ====================================== */
    /*<       D=DEPTH+Z >*/
    d__ = *depth + *z__;
    /*<       P=Y*CD+D*SD >*/
    p = *y * c0_1.cd + d__ * c0_1.sd;
    /*<       Q=Y*SD-D*CD >*/
    q = *y * c0_1.sd - d__ * c0_1.cd;
    /*<       ET(1)=P-AW1 >*/
    et[0] = p - *aw1;
    /*<       ET(2)=P-AW2 >*/
    et[1] = p - *aw2;
    /*<       IF(DABS(Q).LT.EPS)  Q=F0 >*/
    if (fabs(q) < eps) {
        q = f0;
    }
    /*<       IF(DABS(ET(1)).LT.EPS) ET(1)=F0 >*/
    if (fabs(et[0]) < eps) {
        et[0] = f0;
    }
    /*<       IF(DABS(ET(2)).LT.EPS) ET(2)=F0 >*/
    if (fabs(et[1]) < eps) {
        et[1] = f0;
    }
    /* -------------------------------- */
    /* ----- REJECT SINGULAR CASE ----- */
    /* -------------------------------- */
    /* ----- ON FAULT EDGE */
    /*<        >*/
    if (q == f0 && ((xi[0] * xi[1] <= f0 && et[0] * et[1] == f0) ||
                    (et[0] * et[1] <= f0 && xi[0] * xi[1] == f0))) {
        /*<         IRET=1 >*/
        *iret = 1;
        /*<         GO TO 99 >*/
        goto L99;
        /*<       ENDIF >*/
    }
    /* ----- ON NEGATIVE EXTENSION OF FAULT EDGE */
    /*<       KXI(1)=0 >*/
    kxi[0] = 0;
    /*<       KXI(2)=0 >*/
    kxi[1] = 0;
    /*<       KET(1)=0 >*/
    ket[0] = 0;
    /*<       KET(2)=0 >*/
    ket[1] = 0;
    /*<       R12=DSQRT(XI(1)*XI(1)+ET(2)*ET(2)+Q*Q) >*/
    r12 = sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q);
    /*<       R21=DSQRT(XI(2)*XI(2)+ET(1)*ET(1)+Q*Q) >*/
    r21 = sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q);
    /*<       R22=DSQRT(XI(2)*XI(2)+ET(2)*ET(2)+Q*Q) >*/
    r22 = sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q);
    /*<       IF(XI(1).LT.F0 .AND. R21+XI(2).LT.EPS) KXI(1)=1 >*/
    if (xi[0] < f0 && r21 + xi[1] < eps) {
        kxi[0] = 1;
    }
    /*<       IF(XI(1).LT.F0 .AND. R22+XI(2).LT.EPS) KXI(2)=1 >*/
    if (xi[0] < f0 && r22 + xi[1] < eps) {
        kxi[1] = 1;
    }
    /*<       IF(ET(1).LT.F0 .AND. R12+ET(2).LT.EPS) KET(1)=1 >*/
    if (et[0] < f0 && r12 + et[1] < eps) {
        ket[0] = 1;
    }
    /*<       IF(ET(1).LT.F0 .AND. R22+ET(2).LT.EPS) KET(2)=1 >*/
    if (et[0] < f0 && r22 + et[1] < eps) {
        ket[1] = 1;
    }
    /* ===== */
    /*<       DO 223 K=1,2 >*/
    for (k = 1; k <= 2; ++k) {
        /*<       DO 222 J=1,2 >*/
        for (j = 1; j <= 2; ++j) {
            /*<         CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J)) >*/
            dccon2_(&xi[j - 1], &et[k - 1], &q, &c0_1.sd, &c0_1.cd, &kxi[k - 1],
                    &ket[j - 1]);
            /*<         CALL UA(XI(J),ET(K),Q,DD1,DD2,DD3,DUA) >*/
            ua_(&xi[j - 1], &et[k - 1], &q, &dd1, &dd2, &dd3, dua);
            /* ----- */
            /*<         DO 220 I=1,10,3 >*/
            for (i__ = 1; i__ <= 10; i__ += 3) {
                /*<           DU(I)  =-DUA(I) >*/
                du[i__ - 1] = -dua[i__ - 1];
                /*<           DU(I+1)=-DUA(I+1)*CD+DUA(I+2)*SD >*/
                du[i__] = -dua[i__] * c0_1.cd + dua[i__ + 1] * c0_1.sd;
                /*<           DU(I+2)=-DUA(I+1)*SD-DUA(I+2)*CD >*/
                du[i__ + 1] = -dua[i__] * c0_1.sd - dua[i__ + 1] * c0_1.cd;
                /*<           IF(I.LT.10) GO TO 220 >*/
                if (i__ < 10) {
                    goto L220;
                }
                /*<           DU(I)  =-DU(I) >*/
                du[i__ - 1] = -du[i__ - 1];
                /*<           DU(I+1)=-DU(I+1) >*/
                du[i__] = -du[i__];
                /*<           DU(I+2)=-DU(I+2) >*/
                du[i__ + 1] = -du[i__ + 1];
            /*<   220   CONTINUE >*/
            L220:;
            }
            /*<         DO 221 I=1,12 >*/
            for (i__ = 1; i__ <= 12; ++i__) {
                /*<           IF(J+K.NE.3) U(I)=U(I)+DU(I) >*/
                if (j + k != 3) {
                    u[i__ - 1] += du[i__ - 1];
                }
                /*<           IF(J+K.EQ.3) U(I)=U(I)-DU(I) >*/
                if (j + k == 3) {
                    u[i__ - 1] -= du[i__ - 1];
                }
                /*<   221   CONTINUE >*/
                /* L221: */
            }
            /* ----- */
            /*<   222 CONTINUE >*/
            /* L222: */
        }
        /*<   223 CONTINUE >*/
        /* L223: */
    }
    /* ======================================= */
    /* =====  IMAGE-SOURCE CONTRIBUTION  ===== */
    /* ======================================= */
    /*<       D=DEPTH-Z >*/
    d__ = *depth - *z__;
    /*<       P=Y*CD+D*SD >*/
    p = *y * c0_1.cd + d__ * c0_1.sd;
    /*<       Q=Y*SD-D*CD >*/
    q = *y * c0_1.sd - d__ * c0_1.cd;
    /*<       ET(1)=P-AW1 >*/
    et[0] = p - *aw1;
    /*<       ET(2)=P-AW2 >*/
    et[1] = p - *aw2;
    /*<       IF(DABS(Q).LT.EPS)  Q=F0 >*/
    if (fabs(q) < eps) {
        q = f0;
    }
    /*<       IF(DABS(ET(1)).LT.EPS) ET(1)=F0 >*/
    if (fabs(et[0]) < eps) {
        et[0] = f0;
    }
    /*<       IF(DABS(ET(2)).LT.EPS) ET(2)=F0 >*/
    if (fabs(et[1]) < eps) {
        et[1] = f0;
    }
    /* -------------------------------- */
    /* ----- REJECT SINGULAR CASE ----- */
    /* -------------------------------- */
    /* ----- ON FAULT EDGE */
    /*<        >*/
    if (q == f0 && ((xi[0] * xi[1] <= f0 && et[0] * et[1] == f0) ||
                    (et[0] * et[1] <= f0 && xi[0] * xi[1] == f0))) {
        /*<         IRET=1 >*/
        *iret = 1;
        /*<         GO TO 99 >*/
        goto L99;
        /*<       ENDIF >*/
    }
    /* ----- ON NEGATIVE EXTENSION OF FAULT EDGE */
    /*<       KXI(1)=0 >*/
    kxi[0] = 0;
    /*<       KXI(2)=0 >*/
    kxi[1] = 0;
    /*<       KET(1)=0 >*/
    ket[0] = 0;
    /*<       KET(2)=0 >*/
    ket[1] = 0;
    /*<       R12=DSQRT(XI(1)*XI(1)+ET(2)*ET(2)+Q*Q) >*/
    r12 = sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q);
    /*<       R21=DSQRT(XI(2)*XI(2)+ET(1)*ET(1)+Q*Q) >*/
    r21 = sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q);
    /*<       R22=DSQRT(XI(2)*XI(2)+ET(2)*ET(2)+Q*Q) >*/
    r22 = sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q);
    /*<       IF(XI(1).LT.F0 .AND. R21+XI(2).LT.EPS) KXI(1)=1 >*/
    if (xi[0] < f0 && r21 + xi[1] < eps) {
        kxi[0] = 1;
    }
    /*<       IF(XI(1).LT.F0 .AND. R22+XI(2).LT.EPS) KXI(2)=1 >*/
    if (xi[0] < f0 && r22 + xi[1] < eps) {
        kxi[1] = 1;
    }
    /*<       IF(ET(1).LT.F0 .AND. R12+ET(2).LT.EPS) KET(1)=1 >*/
    if (et[0] < f0 && r12 + et[1] < eps) {
        ket[0] = 1;
    }
    /*<       IF(ET(1).LT.F0 .AND. R22+ET(2).LT.EPS) KET(2)=1 >*/
    if (et[0] < f0 && r22 + et[1] < eps) {
        ket[1] = 1;
    }
    /* ===== */
    /*<       DO 334 K=1,2 >*/
    for (k = 1; k <= 2; ++k) {
        /*<       DO 333 J=1,2 >*/
        for (j = 1; j <= 2; ++j) {
            /*<         CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J)) >*/
            dccon2_(&xi[j - 1], &et[k - 1], &q, &c0_1.sd, &c0_1.cd, &kxi[k - 1],
                    &ket[j - 1]);
            /*<         CALL UA(XI(J),ET(K),Q,DD1,DD2,DD3,DUA) >*/
            ua_(&xi[j - 1], &et[k - 1], &q, &dd1, &dd2, &dd3, dua);
            /*<         CALL UB(XI(J),ET(K),Q,DD1,DD2,DD3,DUB) >*/
            ub_(&xi[j - 1], &et[k - 1], &q, &dd1, &dd2, &dd3, dub);
            /*<         CALL UC(XI(J),ET(K),Q,ZZ,DD1,DD2,DD3,DUC) >*/
            uc_(&xi[j - 1], &et[k - 1], &q, &zz, &dd1, &dd2, &dd3, duc);
            /* ----- */
            /*<         DO 330 I=1,10,3 >*/
            for (i__ = 1; i__ <= 10; i__ += 3) {
                /*<           DU(I)=DUA(I)+DUB(I)+Z*DUC(I) >*/
                du[i__ - 1] = dua[i__ - 1] + dub[i__ - 1] + *z__ * duc[i__ - 1];
                /*<        >*/
                du[i__] = (dua[i__] + dub[i__] + *z__ * duc[i__]) * c0_1.cd -
                          (dua[i__ + 1] + dub[i__ + 1] + *z__ * duc[i__ + 1]) *
                              c0_1.sd;
                /*<        >*/
                du[i__ + 1] =
                    (dua[i__] + dub[i__] - *z__ * duc[i__]) * c0_1.sd +
                    (dua[i__ + 1] + dub[i__ + 1] - *z__ * duc[i__ + 1]) *
                        c0_1.cd;
                /*<           IF(I.LT.10) GO TO 330 >*/
                if (i__ < 10) {
                    goto L330;
                }
                /*<           DU(10)=DU(10)+DUC(1) >*/
                du[9] += duc[0];
                /*<           DU(11)=DU(11)+DUC(2)*CD-DUC(3)*SD >*/
                du[10] = du[10] + duc[1] * c0_1.cd - duc[2] * c0_1.sd;
                /*<           DU(12)=DU(12)-DUC(2)*SD-DUC(3)*CD >*/
                du[11] = du[11] - duc[1] * c0_1.sd - duc[2] * c0_1.cd;
            /*<   330   CONTINUE >*/
            L330:;
            }
            /*<         DO 331 I=1,12 >*/
            for (i__ = 1; i__ <= 12; ++i__) {
                /*<           IF(J+K.NE.3) U(I)=U(I)+DU(I) >*/
                if (j + k != 3) {
                    u[i__ - 1] += du[i__ - 1];
                }
                /*<           IF(J+K.EQ.3) U(I)=U(I)-DU(I) >*/
                if (j + k == 3) {
                    u[i__ - 1] -= du[i__ - 1];
                }
                /*<   331   CONTINUE >*/
                /* L331: */
            }
            /* ----- */
            /*<   333 CONTINUE >*/
            /* L333: */
        }
        /*<   334 CONTINUE >*/
        /* L334: */
    }
    /* ===== */
    /*<       UX=U(1) >*/
    *ux = u[0];
    /*<       UY=U(2) >*/
    *uy = u[1];
    /*<       UZ=U(3) >*/
    *uz = u[2];
    /*<       UXX=U(4) >*/
    *uxx = u[3];
    /*<       UYX=U(5) >*/
    *uyx = u[4];
    /*<       UZX=U(6) >*/
    *uzx = u[5];
    /*<       UXY=U(7) >*/
    *uxy = u[6];
    /*<       UYY=U(8) >*/
    *uyy = u[7];
    /*<       UZY=U(9) >*/
    *uzy = u[8];
    /*<       UXZ=U(10) >*/
    *uxz = u[9];
    /*<       UYZ=U(11) >*/
    *uyz = u[10];
    /*<       UZZ=U(12) >*/
    *uzz = u[11];
    /*<       RETURN >*/
    return 0;
/* =========================================== */
/* =====  IN CASE OF SINGULAR (ON EDGE)  ===== */
/* =========================================== */
/*<    99 UX=F0 >*/
L99:
    *ux = f0;
    /*<       UY=F0 >*/
    *uy = f0;
    /*<       UZ=F0 >*/
    *uz = f0;
    /*<       UXX=F0 >*/
    *uxx = f0;
    /*<       UYX=F0 >*/
    *uyx = f0;
    /*<       UZX=F0 >*/
    *uzx = f0;
    /*<       UXY=F0 >*/
    *uxy = f0;
    /*<       UYY=F0 >*/
    *uyy = f0;
    /*<       UZY=F0 >*/
    *uzy = f0;
    /*<       UXZ=F0 >*/
    *uxz = f0;
    /*<       UYZ=F0 >*/
    *uyz = f0;
    /*<       UZZ=F0 >*/
    *uzz = f0;
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* dc3d_ */

/*<       SUBROUTINE  UA(XI,ET,Q,DISL1,DISL2,DISL3,U) >*/
/* Subroutine */ int ua_(doublereal *xi, doublereal *et, doublereal *q,
                         doublereal *disl1, doublereal *disl2,
                         doublereal *disl3, doublereal *u) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f2 = 2.;
    static doublereal pi2 = 6.283185307179586;

    static integer i__;
    static doublereal du[12], qx, qy, xy;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
    /*<       DIMENSION U(12),DU(12) >*/

    /* ******************************************************************** */
    /* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             ***** */
    /* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
    /* ******************************************************************** */

    /* ***** INPUT */
    /* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM */
    /* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS */
    /* ***** OUTPUT */
    /* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /*<       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
     * >*/
    /*<        >*/
    /*<       DATA F0,F2,PI2/0.D0,2.D0,6.283185307179586D0/ >*/
    /* Parameter adjustments */
    --u;

    /* Function Body */
    /* ----- */
    /*<       DO 111  I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<   111 U(I)=F0 >*/
        /* L111: */
        u[i__] = f0;
    }
    /*<       XY=XI*Y11 >*/
    xy = *xi * c2_1.y11;
    /*<       QX=Q *X11 >*/
    qx = *q * c2_1.x11;
    /*<       QY=Q *Y11 >*/
    qy = *q * c2_1.y11;
    /* ====================================== */
    /* =====  STRIKE-SLIP CONTRIBUTION  ===== */
    /* ====================================== */
    /*<       IF(DISL1.NE.F0) THEN >*/
    if (*disl1 != f0) {
        /*<         DU( 1)=    TT/F2 +ALP2*XI*QY >*/
        du[0] = c2_1.tt / f2 + c0_2.alp2 * *xi * qy;
        /*<         DU( 2)=           ALP2*Q/R >*/
        du[1] = c0_2.alp2 * *q / c2_1.r__;
        /*<         DU( 3)= ALP1*ALE -ALP2*Q*QY >*/
        du[2] = c0_2.alp1 * c2_1.ale - c0_2.alp2 * *q * qy;
        /*<         DU( 4)=-ALP1*QY  -ALP2*XI2*Q*Y32 >*/
        du[3] = -c0_2.alp1 * qy - c0_2.alp2 * c2_1.xi2 * *q * c2_1.y32;
        /*<         DU( 5)=          -ALP2*XI*Q/R3 >*/
        du[4] = -c0_2.alp2 * *xi * *q / c2_1.r3;
        /*<         DU( 6)= ALP1*XY  +ALP2*XI*Q2*Y32 >*/
        du[5] = c0_2.alp1 * xy + c0_2.alp2 * *xi * c2_1.q2 * c2_1.y32;
        /*<         DU( 7)= ALP1*XY*SD        +ALP2*XI*FY+D/F2*X11 >*/
        du[6] = c0_2.alp1 * xy * c0_2.sd + c0_2.alp2 * *xi * c2_1.fy +
                c2_1.d__ / f2 * c2_1.x11;
        /*<         DU( 8)=                    ALP2*EY >*/
        du[7] = c0_2.alp2 * c2_1.ey;
        /*<         DU( 9)= ALP1*(CD/R+QY*SD) -ALP2*Q*FY >*/
        du[8] = c0_2.alp1 * (c0_2.cd / c2_1.r__ + qy * c0_2.sd) -
                c0_2.alp2 * *q * c2_1.fy;
        /*<         DU(10)= ALP1*XY*CD        +ALP2*XI*FZ+Y/F2*X11 >*/
        du[9] = c0_2.alp1 * xy * c0_2.cd + c0_2.alp2 * *xi * c2_1.fz +
                c2_1.y / f2 * c2_1.x11;
        /*<         DU(11)=                    ALP2*EZ >*/
        du[10] = c0_2.alp2 * c2_1.ez;
        /*<         DU(12)=-ALP1*(SD/R-QY*CD) -ALP2*Q*FZ >*/
        du[11] = -c0_2.alp1 * (c0_2.sd / c2_1.r__ - qy * c0_2.cd) -
                 c0_2.alp2 * *q * c2_1.fz;
        /*<         DO 222 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   222   U(I)=U(I)+DISL1/PI2*DU(I) >*/
            /* L222: */
            u[i__] += *disl1 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ====================================== */
    /* =====    DIP-SLIP CONTRIBUTION   ===== */
    /* ====================================== */
    /*<       IF(DISL2.NE.F0) THEN >*/
    if (*disl2 != f0) {
        /*<         DU( 1)=           ALP2*Q/R >*/
        du[0] = c0_2.alp2 * *q / c2_1.r__;
        /*<         DU( 2)=    TT/F2 +ALP2*ET*QX >*/
        du[1] = c2_1.tt / f2 + c0_2.alp2 * *et * qx;
        /*<         DU( 3)= ALP1*ALX -ALP2*Q*QX >*/
        du[2] = c0_2.alp1 * c2_1.alx - c0_2.alp2 * *q * qx;
        /*<         DU( 4)=        -ALP2*XI*Q/R3 >*/
        du[3] = -c0_2.alp2 * *xi * *q / c2_1.r3;
        /*<         DU( 5)= -QY/F2 -ALP2*ET*Q/R3 >*/
        du[4] = -qy / f2 - c0_2.alp2 * *et * *q / c2_1.r3;
        /*<         DU( 6)= ALP1/R +ALP2*Q2/R3 >*/
        du[5] = c0_2.alp1 / c2_1.r__ + c0_2.alp2 * c2_1.q2 / c2_1.r3;
        /*<         DU( 7)=                      ALP2*EY >*/
        du[6] = c0_2.alp2 * c2_1.ey;
        /*<         DU( 8)= ALP1*D*X11+XY/F2*SD +ALP2*ET*GY >*/
        du[7] = c0_2.alp1 * c2_1.d__ * c2_1.x11 + xy / f2 * c0_2.sd +
                c0_2.alp2 * *et * c2_1.gy;
        /*<         DU( 9)= ALP1*Y*X11          -ALP2*Q*GY >*/
        du[8] = c0_2.alp1 * c2_1.y * c2_1.x11 - c0_2.alp2 * *q * c2_1.gy;
        /*<         DU(10)=                      ALP2*EZ >*/
        du[9] = c0_2.alp2 * c2_1.ez;
        /*<         DU(11)= ALP1*Y*X11+XY/F2*CD +ALP2*ET*GZ >*/
        du[10] = c0_2.alp1 * c2_1.y * c2_1.x11 + xy / f2 * c0_2.cd +
                 c0_2.alp2 * *et * c2_1.gz;
        /*<         DU(12)=-ALP1*D*X11          -ALP2*Q*GZ >*/
        du[11] = -c0_2.alp1 * c2_1.d__ * c2_1.x11 - c0_2.alp2 * *q * c2_1.gz;
        /*<         DO 333 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   333   U(I)=U(I)+DISL2/PI2*DU(I) >*/
            /* L333: */
            u[i__] += *disl2 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ======================================== */
    /* =====  TENSILE-FAULT CONTRIBUTION  ===== */
    /* ======================================== */
    /*<       IF(DISL3.NE.F0) THEN >*/
    if (*disl3 != f0) {
        /*<         DU( 1)=-ALP1*ALE -ALP2*Q*QY >*/
        du[0] = -c0_2.alp1 * c2_1.ale - c0_2.alp2 * *q * qy;
        /*<         DU( 2)=-ALP1*ALX -ALP2*Q*QX >*/
        du[1] = -c0_2.alp1 * c2_1.alx - c0_2.alp2 * *q * qx;
        /*<         DU( 3)=    TT/F2 -ALP2*(ET*QX+XI*QY) >*/
        du[2] = c2_1.tt / f2 - c0_2.alp2 * (*et * qx + *xi * qy);
        /*<         DU( 4)=-ALP1*XY  +ALP2*XI*Q2*Y32 >*/
        du[3] = -c0_2.alp1 * xy + c0_2.alp2 * *xi * c2_1.q2 * c2_1.y32;
        /*<         DU( 5)=-ALP1/R   +ALP2*Q2/R3 >*/
        du[4] = -c0_2.alp1 / c2_1.r__ + c0_2.alp2 * c2_1.q2 / c2_1.r3;
        /*<         DU( 6)=-ALP1*QY  -ALP2*Q*Q2*Y32 >*/
        du[5] = -c0_2.alp1 * qy - c0_2.alp2 * *q * c2_1.q2 * c2_1.y32;
        /*<         DU( 7)=-ALP1*(CD/R+QY*SD)  -ALP2*Q*FY >*/
        du[6] = -c0_2.alp1 * (c0_2.cd / c2_1.r__ + qy * c0_2.sd) -
                c0_2.alp2 * *q * c2_1.fy;
        /*<         DU( 8)=-ALP1*Y*X11         -ALP2*Q*GY >*/
        du[7] = -c0_2.alp1 * c2_1.y * c2_1.x11 - c0_2.alp2 * *q * c2_1.gy;
        /*<         DU( 9)= ALP1*(D*X11+XY*SD) +ALP2*Q*HY >*/
        du[8] = c0_2.alp1 * (c2_1.d__ * c2_1.x11 + xy * c0_2.sd) +
                c0_2.alp2 * *q * c2_1.hy;
        /*<         DU(10)= ALP1*(SD/R-QY*CD)  -ALP2*Q*FZ >*/
        du[9] = c0_2.alp1 * (c0_2.sd / c2_1.r__ - qy * c0_2.cd) -
                c0_2.alp2 * *q * c2_1.fz;
        /*<         DU(11)= ALP1*D*X11         -ALP2*Q*GZ >*/
        du[10] = c0_2.alp1 * c2_1.d__ * c2_1.x11 - c0_2.alp2 * *q * c2_1.gz;
        /*<         DU(12)= ALP1*(Y*X11+XY*CD) +ALP2*Q*HZ >*/
        du[11] = c0_2.alp1 * (c2_1.y * c2_1.x11 + xy * c0_2.cd) +
                 c0_2.alp2 * *q * c2_1.hz;
        /*<         DO 444 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   444   U(I)=U(I)+DISL3/PI2*DU(I) >*/
            /* L444: */
            u[i__] += *disl3 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* ua_ */

/*<       SUBROUTINE  UB(XI,ET,Q,DISL1,DISL2,DISL3,U) >*/
/* Subroutine */ int ub_(doublereal *xi, doublereal *et, doublereal *q,
                         doublereal *disl1, doublereal *disl2,
                         doublereal *disl3, doublereal *u) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal pi2 = 6.283185307179586;

    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal), log(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal x, d11, rd, du[12], qx, qy, xy, ai1, ai2, aj2, ai4, ai3,
        aj5, ak1, ak3, aj3, aj6, ak2, ak4, aj1, rd2, aj4;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
    /*<       DIMENSION U(12),DU(12) >*/

    /* ******************************************************************** */
    /* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             ***** */
    /* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
    /* ******************************************************************** */

    /* ***** INPUT */
    /* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM */
    /* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS */
    /* ***** OUTPUT */
    /* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /*<       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
     * >*/
    /*<        >*/
    /*<       DATA  F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/ >*/
    /* Parameter adjustments */
    --u;

    /* Function Body */
    /* ----- */
    /*<       RD=R+D >*/
    rd = c2_1.r__ + c2_1.d__;
    /*<       D11=F1/(R*RD) >*/
    d11 = f1 / (c2_1.r__ * rd);
    /*<       AJ2=XI*Y/RD*D11 >*/
    aj2 = *xi * c2_1.y / rd * d11;
    /*<       AJ5=-(D+Y*Y/RD)*D11 >*/
    aj5 = -(c2_1.d__ + c2_1.y * c2_1.y / rd) * d11;
    /*<       IF(CD.NE.F0) THEN >*/
    if (c0_2.cd != f0) {
        /*<         IF(XI.EQ.F0) THEN >*/
        if (*xi == f0) {
            /*<           AI4=F0 >*/
            ai4 = f0;
            /*<         ELSE >*/
        } else {
            /*<           X=DSQRT(XI2+Q2) >*/
            x = sqrt(c2_1.xi2 + c2_1.q2);
            /*<        >*/
            ai4 = f1 / c0_2.cdcd *
                  (*xi / rd * c0_2.sdcd +
                   f2 * atan((*et * (x + *q * c0_2.cd) +
                              x * (c2_1.r__ + x) * c0_2.sd) /
                             (*xi * (c2_1.r__ + x) * c0_2.cd)));
            /*<         ENDIF >*/
        }
        /*<         AI3=(Y*CD/RD-ALE+SD*DLOG(RD))/CDCD >*/
        ai3 =
            (c2_1.y * c0_2.cd / rd - c2_1.ale + c0_2.sd * log(rd)) / c0_2.cdcd;
        /*<         AK1=XI*(D11-Y11*SD)/CD >*/
        ak1 = *xi * (d11 - c2_1.y11 * c0_2.sd) / c0_2.cd;
        /*<         AK3=(Q*Y11-Y*D11)/CD >*/
        ak3 = (*q * c2_1.y11 - c2_1.y * d11) / c0_2.cd;
        /*<         AJ3=(AK1-AJ2*SD)/CD >*/
        aj3 = (ak1 - aj2 * c0_2.sd) / c0_2.cd;
        /*<         AJ6=(AK3-AJ5*SD)/CD >*/
        aj6 = (ak3 - aj5 * c0_2.sd) / c0_2.cd;
        /*<       ELSE >*/
    } else {
        /*<         RD2=RD*RD >*/
        rd2 = rd * rd;
        /*<         AI3=(ET/RD+Y*Q/RD2-ALE)/F2 >*/
        ai3 = (*et / rd + c2_1.y * *q / rd2 - c2_1.ale) / f2;
        /*<         AI4=XI*Y/RD2/F2 >*/
        ai4 = *xi * c2_1.y / rd2 / f2;
        /*<         AK1=XI*Q/RD*D11 >*/
        ak1 = *xi * *q / rd * d11;
        /*<         AK3=SD/RD*(XI2*D11-F1) >*/
        ak3 = c0_2.sd / rd * (c2_1.xi2 * d11 - f1);
        /*<         AJ3=-XI/RD2*(Q2*D11-F1/F2) >*/
        aj3 = -(*xi) / rd2 * (c2_1.q2 * d11 - f1 / f2);
        /*<         AJ6=-Y/RD2*(XI2*D11-F1/F2) >*/
        aj6 = -c2_1.y / rd2 * (c2_1.xi2 * d11 - f1 / f2);
        /*<       ENDIF >*/
    }
    /* ----- */
    /*<       XY=XI*Y11 >*/
    xy = *xi * c2_1.y11;
    /*<       AI1=-XI/RD*CD-AI4*SD >*/
    ai1 = -(*xi) / rd * c0_2.cd - ai4 * c0_2.sd;
    /*<       AI2= DLOG(RD)+AI3*SD >*/
    ai2 = log(rd) + ai3 * c0_2.sd;
    /*<       AK2= F1/R+AK3*SD >*/
    ak2 = f1 / c2_1.r__ + ak3 * c0_2.sd;
    /*<       AK4= XY*CD-AK1*SD >*/
    ak4 = xy * c0_2.cd - ak1 * c0_2.sd;
    /*<       AJ1= AJ5*CD-AJ6*SD >*/
    aj1 = aj5 * c0_2.cd - aj6 * c0_2.sd;
    /*<       AJ4=-XY-AJ2*CD+AJ3*SD >*/
    aj4 = -xy - aj2 * c0_2.cd + aj3 * c0_2.sd;
    /* ===== */
    /*<       DO 111  I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<   111 U(I)=F0 >*/
        /* L111: */
        u[i__] = f0;
    }
    /*<       QX=Q*X11 >*/
    qx = *q * c2_1.x11;
    /*<       QY=Q*Y11 >*/
    qy = *q * c2_1.y11;
    /* ====================================== */
    /* =====  STRIKE-SLIP CONTRIBUTION  ===== */
    /* ====================================== */
    /*<       IF(DISL1.NE.F0) THEN >*/
    if (*disl1 != f0) {
        /*<         DU( 1)=-XI*QY-TT -ALP3*AI1*SD >*/
        du[0] = -(*xi) * qy - c2_1.tt - c0_2.alp3 * ai1 * c0_2.sd;
        /*<         DU( 2)=-Q/R      +ALP3*Y/RD*SD >*/
        du[1] = -(*q) / c2_1.r__ + c0_2.alp3 * c2_1.y / rd * c0_2.sd;
        /*<         DU( 3)= Q*QY     -ALP3*AI2*SD >*/
        du[2] = *q * qy - c0_2.alp3 * ai2 * c0_2.sd;
        /*<         DU( 4)= XI2*Q*Y32 -ALP3*AJ1*SD >*/
        du[3] = c2_1.xi2 * *q * c2_1.y32 - c0_2.alp3 * aj1 * c0_2.sd;
        /*<         DU( 5)= XI*Q/R3   -ALP3*AJ2*SD >*/
        du[4] = *xi * *q / c2_1.r3 - c0_2.alp3 * aj2 * c0_2.sd;
        /*<         DU( 6)=-XI*Q2*Y32 -ALP3*AJ3*SD >*/
        du[5] = -(*xi) * c2_1.q2 * c2_1.y32 - c0_2.alp3 * aj3 * c0_2.sd;
        /*<         DU( 7)=-XI*FY-D*X11 +ALP3*(XY+AJ4)*SD >*/
        du[6] = -(*xi) * c2_1.fy - c2_1.d__ * c2_1.x11 +
                c0_2.alp3 * (xy + aj4) * c0_2.sd;
        /*<         DU( 8)=-EY          +ALP3*(F1/R+AJ5)*SD >*/
        du[7] = -c2_1.ey + c0_2.alp3 * (f1 / c2_1.r__ + aj5) * c0_2.sd;
        /*<         DU( 9)= Q*FY        -ALP3*(QY-AJ6)*SD >*/
        du[8] = *q * c2_1.fy - c0_2.alp3 * (qy - aj6) * c0_2.sd;
        /*<         DU(10)=-XI*FZ-Y*X11 +ALP3*AK1*SD >*/
        du[9] =
            -(*xi) * c2_1.fz - c2_1.y * c2_1.x11 + c0_2.alp3 * ak1 * c0_2.sd;
        /*<         DU(11)=-EZ          +ALP3*Y*D11*SD >*/
        du[10] = -c2_1.ez + c0_2.alp3 * c2_1.y * d11 * c0_2.sd;
        /*<         DU(12)= Q*FZ        +ALP3*AK2*SD >*/
        du[11] = *q * c2_1.fz + c0_2.alp3 * ak2 * c0_2.sd;
        /*<         DO 222 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   222   U(I)=U(I)+DISL1/PI2*DU(I) >*/
            /* L222: */
            u[i__] += *disl1 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ====================================== */
    /* =====    DIP-SLIP CONTRIBUTION   ===== */
    /* ====================================== */
    /*<       IF(DISL2.NE.F0) THEN >*/
    if (*disl2 != f0) {
        /*<         DU( 1)=-Q/R      +ALP3*AI3*SDCD >*/
        du[0] = -(*q) / c2_1.r__ + c0_2.alp3 * ai3 * c0_2.sdcd;
        /*<         DU( 2)=-ET*QX-TT -ALP3*XI/RD*SDCD >*/
        du[1] = -(*et) * qx - c2_1.tt - c0_2.alp3 * *xi / rd * c0_2.sdcd;
        /*<         DU( 3)= Q*QX     +ALP3*AI4*SDCD >*/
        du[2] = *q * qx + c0_2.alp3 * ai4 * c0_2.sdcd;
        /*<         DU( 4)= XI*Q/R3     +ALP3*AJ4*SDCD >*/
        du[3] = *xi * *q / c2_1.r3 + c0_2.alp3 * aj4 * c0_2.sdcd;
        /*<         DU( 5)= ET*Q/R3+QY  +ALP3*AJ5*SDCD >*/
        du[4] = *et * *q / c2_1.r3 + qy + c0_2.alp3 * aj5 * c0_2.sdcd;
        /*<         DU( 6)=-Q2/R3       +ALP3*AJ6*SDCD >*/
        du[5] = -c2_1.q2 / c2_1.r3 + c0_2.alp3 * aj6 * c0_2.sdcd;
        /*<         DU( 7)=-EY          +ALP3*AJ1*SDCD >*/
        du[6] = -c2_1.ey + c0_2.alp3 * aj1 * c0_2.sdcd;
        /*<         DU( 8)=-ET*GY-XY*SD +ALP3*AJ2*SDCD >*/
        du[7] = -(*et) * c2_1.gy - xy * c0_2.sd + c0_2.alp3 * aj2 * c0_2.sdcd;
        /*<         DU( 9)= Q*GY        +ALP3*AJ3*SDCD >*/
        du[8] = *q * c2_1.gy + c0_2.alp3 * aj3 * c0_2.sdcd;
        /*<         DU(10)=-EZ          -ALP3*AK3*SDCD >*/
        du[9] = -c2_1.ez - c0_2.alp3 * ak3 * c0_2.sdcd;
        /*<         DU(11)=-ET*GZ-XY*CD -ALP3*XI*D11*SDCD >*/
        du[10] =
            -(*et) * c2_1.gz - xy * c0_2.cd - c0_2.alp3 * *xi * d11 * c0_2.sdcd;
        /*<         DU(12)= Q*GZ        -ALP3*AK4*SDCD >*/
        du[11] = *q * c2_1.gz - c0_2.alp3 * ak4 * c0_2.sdcd;
        /*<         DO 333 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   333   U(I)=U(I)+DISL2/PI2*DU(I) >*/
            /* L333: */
            u[i__] += *disl2 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ======================================== */
    /* =====  TENSILE-FAULT CONTRIBUTION  ===== */
    /* ======================================== */
    /*<       IF(DISL3.NE.F0) THEN >*/
    if (*disl3 != f0) {
        /*<         DU( 1)= Q*QY           -ALP3*AI3*SDSD >*/
        du[0] = *q * qy - c0_2.alp3 * ai3 * c0_2.sdsd;
        /*<         DU( 2)= Q*QX           +ALP3*XI/RD*SDSD >*/
        du[1] = *q * qx + c0_2.alp3 * *xi / rd * c0_2.sdsd;
        /*<         DU( 3)= ET*QX+XI*QY-TT -ALP3*AI4*SDSD >*/
        du[2] = *et * qx + *xi * qy - c2_1.tt - c0_2.alp3 * ai4 * c0_2.sdsd;
        /*<         DU( 4)=-XI*Q2*Y32 -ALP3*AJ4*SDSD >*/
        du[3] = -(*xi) * c2_1.q2 * c2_1.y32 - c0_2.alp3 * aj4 * c0_2.sdsd;
        /*<         DU( 5)=-Q2/R3     -ALP3*AJ5*SDSD >*/
        du[4] = -c2_1.q2 / c2_1.r3 - c0_2.alp3 * aj5 * c0_2.sdsd;
        /*<         DU( 6)= Q*Q2*Y32  -ALP3*AJ6*SDSD >*/
        du[5] = *q * c2_1.q2 * c2_1.y32 - c0_2.alp3 * aj6 * c0_2.sdsd;
        /*<         DU( 7)= Q*FY -ALP3*AJ1*SDSD >*/
        du[6] = *q * c2_1.fy - c0_2.alp3 * aj1 * c0_2.sdsd;
        /*<         DU( 8)= Q*GY -ALP3*AJ2*SDSD >*/
        du[7] = *q * c2_1.gy - c0_2.alp3 * aj2 * c0_2.sdsd;
        /*<         DU( 9)=-Q*HY -ALP3*AJ3*SDSD >*/
        du[8] = -(*q) * c2_1.hy - c0_2.alp3 * aj3 * c0_2.sdsd;
        /*<         DU(10)= Q*FZ +ALP3*AK3*SDSD >*/
        du[9] = *q * c2_1.fz + c0_2.alp3 * ak3 * c0_2.sdsd;
        /*<         DU(11)= Q*GZ +ALP3*XI*D11*SDSD >*/
        du[10] = *q * c2_1.gz + c0_2.alp3 * *xi * d11 * c0_2.sdsd;
        /*<         DU(12)=-Q*HZ +ALP3*AK4*SDSD >*/
        du[11] = -(*q) * c2_1.hz + c0_2.alp3 * ak4 * c0_2.sdsd;
        /*<         DO 444 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   444   U(I)=U(I)+DISL3/PI2*DU(I) >*/
            /* L444: */
            u[i__] += *disl3 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* ub_ */

/*<       SUBROUTINE  UC(XI,ET,Q,Z,DISL1,DISL2,DISL3,U) >*/
/* Subroutine */ int uc_(doublereal *xi, doublereal *et, doublereal *q,
                         doublereal *z__, doublereal *disl1, doublereal *disl2,
                         doublereal *disl3, doublereal *u) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal f3 = 3.;
    static doublereal pi2 = 6.283185307179586;

    static doublereal c__, h__;
    static integer i__;
    static doublereal y0, z0, du[12], x53, y53, z32, z53, qq, qx, qy, qr, xy,
        yy0, cdr, cqx, ppy, ppz, qqy, qqz;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
    /*<       DIMENSION U(12),DU(12) >*/

    /* ******************************************************************** */
    /* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             ***** */
    /* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
    /* ******************************************************************** */

    /* ***** INPUT */
    /* *****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM */
    /* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS */
    /* ***** OUTPUT */
    /* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /*<       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
     * >*/
    /*<        >*/
    /*<       DATA F0,F1,F2,F3,PI2/0.D0,1.D0,2.D0,3.D0,6.283185307179586D0/ >*/
    /* Parameter adjustments */
    --u;

    /* Function Body */
    /* ----- */
    /*<       C=D+Z >*/
    c__ = c2_1.d__ + *z__;
    /*<       X53=(8.D0*R2+9.D0*R*XI+F3*XI2)*X11*X11*X11/R2 >*/
    x53 = (c2_1.r2 * 8. + c2_1.r__ * 9. * *xi + f3 * c2_1.xi2) * c2_1.x11 *
          c2_1.x11 * c2_1.x11 / c2_1.r2;
    /*<       Y53=(8.D0*R2+9.D0*R*ET+F3*ET2)*Y11*Y11*Y11/R2 >*/
    y53 = (c2_1.r2 * 8. + c2_1.r__ * 9. * *et + f3 * c2_1.et2) * c2_1.y11 *
          c2_1.y11 * c2_1.y11 / c2_1.r2;
    /*<       H=Q*CD-Z >*/
    h__ = *q * c0_2.cd - *z__;
    /*<       Z32=SD/R3-H*Y32 >*/
    z32 = c0_2.sd / c2_1.r3 - h__ * c2_1.y32;
    /*<       Z53=F3*SD/R5-H*Y53 >*/
    z53 = f3 * c0_2.sd / c2_1.r5 - h__ * y53;
    /*<       Y0=Y11-XI2*Y32 >*/
    y0 = c2_1.y11 - c2_1.xi2 * c2_1.y32;
    /*<       Z0=Z32-XI2*Z53 >*/
    z0 = z32 - c2_1.xi2 * z53;
    /*<       PPY=CD/R3+Q*Y32*SD >*/
    ppy = c0_2.cd / c2_1.r3 + *q * c2_1.y32 * c0_2.sd;
    /*<       PPZ=SD/R3-Q*Y32*CD >*/
    ppz = c0_2.sd / c2_1.r3 - *q * c2_1.y32 * c0_2.cd;
    /*<       QQ=Z*Y32+Z32+Z0 >*/
    qq = *z__ * c2_1.y32 + z32 + z0;
    /*<       QQY=F3*C*D/R5-QQ*SD >*/
    qqy = f3 * c__ * c2_1.d__ / c2_1.r5 - qq * c0_2.sd;
    /*<       QQZ=F3*C*Y/R5-QQ*CD+Q*Y32 >*/
    qqz = f3 * c__ * c2_1.y / c2_1.r5 - qq * c0_2.cd + *q * c2_1.y32;
    /*<       XY=XI*Y11 >*/
    xy = *xi * c2_1.y11;
    /*<       QX=Q*X11 >*/
    qx = *q * c2_1.x11;
    /*<       QY=Q*Y11 >*/
    qy = *q * c2_1.y11;
    /*<       QR=F3*Q/R5 >*/
    qr = f3 * *q / c2_1.r5;
    /*<       CQX=C*Q*X53 >*/
    cqx = c__ * *q * x53;
    /*<       CDR=(C+D)/R3 >*/
    cdr = (c__ + c2_1.d__) / c2_1.r3;
    /*<       YY0=Y/R3-Y0*CD >*/
    yy0 = c2_1.y / c2_1.r3 - y0 * c0_2.cd;
    /* ===== */
    /*<       DO 111  I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<   111 U(I)=F0 >*/
        /* L111: */
        u[i__] = f0;
    }
    /* ====================================== */
    /* =====  STRIKE-SLIP CONTRIBUTION  ===== */
    /* ====================================== */
    /*<       IF(DISL1.NE.F0) THEN >*/
    if (*disl1 != f0) {
        /*<         DU( 1)= ALP4*XY*CD           -ALP5*XI*Q*Z32 >*/
        du[0] = c0_2.alp4 * xy * c0_2.cd - c0_2.alp5 * *xi * *q * z32;
        /*<         DU( 2)= ALP4*(CD/R+F2*QY*SD) -ALP5*C*Q/R3 >*/
        du[1] = c0_2.alp4 * (c0_2.cd / c2_1.r__ + f2 * qy * c0_2.sd) -
                c0_2.alp5 * c__ * *q / c2_1.r3;
        /*<         DU( 3)= ALP4*QY*CD           -ALP5*(C*ET/R3-Z*Y11+XI2*Z32)
         * >*/
        du[2] = c0_2.alp4 * qy * c0_2.cd -
                c0_2.alp5 *
                    (c__ * *et / c2_1.r3 - *z__ * c2_1.y11 + c2_1.xi2 * z32);
        /*<         DU( 4)= ALP4*Y0*CD                  -ALP5*Q*Z0 >*/
        du[3] = c0_2.alp4 * y0 * c0_2.cd - c0_2.alp5 * *q * z0;
        /*<         DU( 5)=-ALP4*XI*(CD/R3+F2*Q*Y32*SD) +ALP5*C*XI*QR >*/
        du[4] = -c0_2.alp4 * *xi *
                    (c0_2.cd / c2_1.r3 + f2 * *q * c2_1.y32 * c0_2.sd) +
                c0_2.alp5 * c__ * *xi * qr;
        /*<         DU( 6)=-ALP4*XI*Q*Y32*CD            +ALP5*XI*(F3*C*ET/R5-QQ)
         * >*/
        du[5] = -c0_2.alp4 * *xi * *q * c2_1.y32 * c0_2.cd +
                c0_2.alp5 * *xi * (f3 * c__ * *et / c2_1.r5 - qq);
        /*<         DU( 7)=-ALP4*XI*PPY*CD    -ALP5*XI*QQY >*/
        du[6] = -c0_2.alp4 * *xi * ppy * c0_2.cd - c0_2.alp5 * *xi * qqy;
        /*<        >*/
        du[7] = c0_2.alp4 * f2 * (c2_1.d__ / c2_1.r3 - y0 * c0_2.sd) * c0_2.sd -
                c2_1.y / c2_1.r3 * c0_2.cd -
                c0_2.alp5 * (cdr * c0_2.sd - *et / c2_1.r3 - c__ * c2_1.y * qr);
        /*<         DU( 9)=-ALP4*Q/R3+YY0*SD
         * +ALP5*(CDR*CD+C*D*QR-(Y0*CD+Q*Z0)*SD)
         * >*/
        du[8] = -c0_2.alp4 * *q / c2_1.r3 + yy0 * c0_2.sd +
                c0_2.alp5 * (cdr * c0_2.cd + c__ * c2_1.d__ * qr -
                             (y0 * c0_2.cd + *q * z0) * c0_2.sd);
        /*<         DU(10)= ALP4*XI*PPZ*CD    -ALP5*XI*QQZ >*/
        du[9] = c0_2.alp4 * *xi * ppz * c0_2.cd - c0_2.alp5 * *xi * qqz;
        /*<         DU(11)= ALP4*F2*(Y/R3-Y0*CD)*SD+D/R3*CD
         * -ALP5*(CDR*CD+C*D*QR)
         * >*/
        du[10] = c0_2.alp4 * f2 * (c2_1.y / c2_1.r3 - y0 * c0_2.cd) * c0_2.sd +
                 c2_1.d__ / c2_1.r3 * c0_2.cd -
                 c0_2.alp5 * (cdr * c0_2.cd + c__ * c2_1.d__ * qr);
        /*<         DU(12)=         YY0*CD -ALP5*(CDR*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD)
         * >*/
        du[11] =
            yy0 * c0_2.cd - c0_2.alp5 * (cdr * c0_2.sd - c__ * c2_1.y * qr -
                                         y0 * c0_2.sdsd + *q * z0 * c0_2.cd);
        /*<         DO 222 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   222   U(I)=U(I)+DISL1/PI2*DU(I) >*/
            /* L222: */
            u[i__] += *disl1 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ====================================== */
    /* =====    DIP-SLIP CONTRIBUTION   ===== */
    /* ====================================== */
    /*<       IF(DISL2.NE.F0) THEN >*/
    if (*disl2 != f0) {
        /*<         DU( 1)= ALP4*CD/R -QY*SD -ALP5*C*Q/R3 >*/
        du[0] = c0_2.alp4 * c0_2.cd / c2_1.r__ - qy * c0_2.sd -
                c0_2.alp5 * c__ * *q / c2_1.r3;
        /*<         DU( 2)= ALP4*Y*X11       -ALP5*C*ET*Q*X32 >*/
        du[1] = c0_2.alp4 * c2_1.y * c2_1.x11 -
                c0_2.alp5 * c__ * *et * *q * c2_1.x32;
        /*<         DU( 3)=     -D*X11-XY*SD -ALP5*C*(X11-Q2*X32) >*/
        du[2] = -c2_1.d__ * c2_1.x11 - xy * c0_2.sd -
                c0_2.alp5 * c__ * (c2_1.x11 - c2_1.q2 * c2_1.x32);
        /*<         DU( 4)=-ALP4*XI/R3*CD +ALP5*C*XI*QR +XI*Q*Y32*SD >*/
        du[3] = -c0_2.alp4 * *xi / c2_1.r3 * c0_2.cd +
                c0_2.alp5 * c__ * *xi * qr + *xi * *q * c2_1.y32 * c0_2.sd;
        /*<         DU( 5)=-ALP4*Y/R3     +ALP5*C*ET*QR >*/
        du[4] = -c0_2.alp4 * c2_1.y / c2_1.r3 + c0_2.alp5 * c__ * *et * qr;
        /*<         DU( 6)=    D/R3-Y0*SD +ALP5*C/R3*(F1-F3*Q2/R2) >*/
        du[5] = c2_1.d__ / c2_1.r3 - y0 * c0_2.sd +
                c0_2.alp5 * c__ / c2_1.r3 * (f1 - f3 * c2_1.q2 / c2_1.r2);
        /*<         DU( 7)=-ALP4*ET/R3+Y0*SDSD -ALP5*(CDR*SD-C*Y*QR) >*/
        du[6] = -c0_2.alp4 * *et / c2_1.r3 + y0 * c0_2.sdsd -
                c0_2.alp5 * (cdr * c0_2.sd - c__ * c2_1.y * qr);
        /*<         DU( 8)= ALP4*(X11-Y*Y*X32)
         * -ALP5*C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)
         * >*/
        du[7] = c0_2.alp4 * (c2_1.x11 - c2_1.y * c2_1.y * c2_1.x32) -
                c0_2.alp5 * c__ *
                    ((c2_1.d__ + f2 * *q * c0_2.cd) * c2_1.x32 -
                     c2_1.y * *et * *q * x53);
        /*<         DU( 9)=  XI*PPY*SD+Y*D*X32
         * +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)
         * >*/
        du[8] = *xi * ppy * c0_2.sd + c2_1.y * c2_1.d__ * c2_1.x32 +
                c0_2.alp5 * c__ *
                    ((c2_1.y + f2 * *q * c0_2.sd) * c2_1.x32 -
                     c2_1.y * c2_1.q2 * x53);
        /*<         DU(10)=      -Q/R3+Y0*SDCD -ALP5*(CDR*CD+C*D*QR) >*/
        du[9] = -(*q) / c2_1.r3 + y0 * c0_2.sdcd -
                c0_2.alp5 * (cdr * c0_2.cd + c__ * c2_1.d__ * qr);
        /*<         DU(11)= ALP4*Y*D*X32 -ALP5*C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)
         * >*/
        du[10] = c0_2.alp4 * c2_1.y * c2_1.d__ * c2_1.x32 -
                 c0_2.alp5 * c__ *
                     ((c2_1.y - f2 * *q * c0_2.sd) * c2_1.x32 +
                      c2_1.d__ * *et * *q * x53);
        /*< DU(12)=-XI*PPZ*SD+X11-D*D*X32-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53)
         * >*/
        du[11] = -(*xi) * ppz * c0_2.sd + c2_1.x11 -
                 c2_1.d__ * c2_1.d__ * c2_1.x32 -
                 c0_2.alp5 * c__ *
                     ((c2_1.d__ - f2 * *q * c0_2.cd) * c2_1.x32 -
                      c2_1.d__ * c2_1.q2 * x53);
        /*<         DO 333 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   333   U(I)=U(I)+DISL2/PI2*DU(I) >*/
            /* L333: */
            u[i__] += *disl2 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ======================================== */
    /* =====  TENSILE-FAULT CONTRIBUTION  ===== */
    /* ======================================== */
    /*<       IF(DISL3.NE.F0) THEN >*/
    if (*disl3 != f0) {
        /*<         DU( 1)=-ALP4*(SD/R+QY*CD)   -ALP5*(Z*Y11-Q2*Z32) >*/
        du[0] = -c0_2.alp4 * (c0_2.sd / c2_1.r__ + qy * c0_2.cd) -
                c0_2.alp5 * (*z__ * c2_1.y11 - c2_1.q2 * z32);
        /*<         DU( 2)= ALP4*F2*XY*SD+D*X11 -ALP5*C*(X11-Q2*X32) >*/
        du[1] = c0_2.alp4 * f2 * xy * c0_2.sd + c2_1.d__ * c2_1.x11 -
                c0_2.alp5 * c__ * (c2_1.x11 - c2_1.q2 * c2_1.x32);
        /*<         DU( 3)= ALP4*(Y*X11+XY*CD)  +ALP5*Q*(C*ET*X32+XI*Z32) >*/
        du[2] = c0_2.alp4 * (c2_1.y * c2_1.x11 + xy * c0_2.cd) +
                c0_2.alp5 * *q * (c__ * *et * c2_1.x32 + *xi * z32);
        /*<         DU( 4)=
         * ALP4*XI/R3*SD+XI*Q*Y32*CD+ALP5*XI*(F3*C*ET/R5-F2*Z32-Z0)
         * >*/
        du[3] = c0_2.alp4 * *xi / c2_1.r3 * c0_2.sd +
                *xi * *q * c2_1.y32 * c0_2.cd +
                c0_2.alp5 * *xi * (f3 * c__ * *et / c2_1.r5 - f2 * z32 - z0);
        /*<         DU( 5)= ALP4*F2*Y0*SD-D/R3 +ALP5*C/R3*(F1-F3*Q2/R2) >*/
        du[4] = c0_2.alp4 * f2 * y0 * c0_2.sd - c2_1.d__ / c2_1.r3 +
                c0_2.alp5 * c__ / c2_1.r3 * (f1 - f3 * c2_1.q2 / c2_1.r2);
        /*<         DU( 6)=-ALP4*YY0           -ALP5*(C*ET*QR-Q*Z0) >*/
        du[5] = -c0_2.alp4 * yy0 - c0_2.alp5 * (c__ * *et * qr - *q * z0);
        /*<         DU( 7)= ALP4*(Q/R3+Y0*SDCD)   +ALP5*(Z/R3*CD+C*D*QR-Q*Z0*SD)
         * >*/
        du[6] = c0_2.alp4 * (*q / c2_1.r3 + y0 * c0_2.sdcd) +
                c0_2.alp5 * (*z__ / c2_1.r3 * c0_2.cd + c__ * c2_1.d__ * qr -
                             *q * z0 * c0_2.sd);
        /*<        >*/
        du[7] = -c0_2.alp4 * f2 * *xi * ppy * c0_2.sd -
                c2_1.y * c2_1.d__ * c2_1.x32 +
                c0_2.alp5 * c__ *
                    ((c2_1.y + f2 * *q * c0_2.sd) * c2_1.x32 -
                     c2_1.y * c2_1.q2 * x53);
        /*<        >*/
        du[8] = -c0_2.alp4 * (*xi * ppy * c0_2.cd - c2_1.x11 +
                              c2_1.y * c2_1.y * c2_1.x32) +
                c0_2.alp5 * (c__ * ((c2_1.d__ + f2 * *q * c0_2.cd) * c2_1.x32 -
                                    c2_1.y * *et * *q * x53) +
                             *xi * qqy);
        /*<         DU(10)=  -ET/R3+Y0*CDCD
         * -ALP5*(Z/R3*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD)
         * >*/
        du[9] = -(*et) / c2_1.r3 + y0 * c0_2.cdcd -
                c0_2.alp5 * (*z__ / c2_1.r3 * c0_2.sd - c__ * c2_1.y * qr -
                             y0 * c0_2.sdsd + *q * z0 * c0_2.cd);
        /*<        >*/
        du[10] = c0_2.alp4 * f2 * *xi * ppz * c0_2.sd - c2_1.x11 +
                 c2_1.d__ * c2_1.d__ * c2_1.x32 -
                 c0_2.alp5 * c__ *
                     ((c2_1.d__ - f2 * *q * c0_2.cd) * c2_1.x32 -
                      c2_1.d__ * c2_1.q2 * x53);
        /*<        >*/
        du[11] =
            c0_2.alp4 * (*xi * ppz * c0_2.cd + c2_1.y * c2_1.d__ * c2_1.x32) +
            c0_2.alp5 * (c__ * ((c2_1.y - f2 * *q * c0_2.sd) * c2_1.x32 +
                                c2_1.d__ * *et * *q * x53) +
                         *xi * qqz);
        /*<         DO 444 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   444   U(I)=U(I)+DISL3/PI2*DU(I) >*/
            /* L444: */
            u[i__] += *disl3 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* uc_ */

/*<       SUBROUTINE  DCCON0(ALPHA,DIP) >*/
/* Subroutine */ int dccon0_(doublereal *alpha, doublereal *dip) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal pi2 = 6.283185307179586;
    static doublereal eps = 1e-6;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal p18;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/

    /* ******************************************************************* */
    /* *****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    ***** */
    /* ******************************************************************* */

    /* ***** INPUT */
    /* *****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU) */
    /* *****   DIP   : DIP-ANGLE (DEGREE) */
    /* ### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO */

    /*<       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
     * >*/
    /*<       DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/ >*/
    /*<       DATA EPS/1.D-6/ >*/
    /* ----- */
    /*<       ALP1=(F1-ALPHA)/F2 >*/
    c0_2.alp1 = (f1 - *alpha) / f2;
    /*<       ALP2= ALPHA/F2 >*/
    c0_2.alp2 = *alpha / f2;
    /*<       ALP3=(F1-ALPHA)/ALPHA >*/
    c0_2.alp3 = (f1 - *alpha) / *alpha;
    /*<       ALP4= F1-ALPHA >*/
    c0_2.alp4 = f1 - *alpha;
    /*<       ALP5= ALPHA >*/
    c0_2.alp5 = *alpha;
    /* ----- */
    /*<       P18=PI2/360.D0 >*/
    p18 = pi2 / 360.;
    /*<       SD=DSIN(DIP*P18) >*/
    c0_2.sd = sin(*dip * p18);
    /*<       CD=DCOS(DIP*P18) >*/
    c0_2.cd = cos(*dip * p18);
    /*<       IF(DABS(CD).LT.EPS) THEN >*/
    if (fabs(c0_2.cd) < eps) {
        /*<         CD=F0 >*/
        c0_2.cd = f0;
        /*<         IF(SD.GT.F0) SD= F1 >*/
        if (c0_2.sd > f0) {
            c0_2.sd = f1;
        }
        /*<         IF(SD.LT.F0) SD=-F1 >*/
        if (c0_2.sd < f0) {
            c0_2.sd = -f1;
        }
        /*<       ENDIF >*/
    }
    /*<       SDSD=SD*SD >*/
    c0_2.sdsd = c0_2.sd * c0_2.sd;
    /*<       CDCD=CD*CD >*/
    c0_2.cdcd = c0_2.cd * c0_2.cd;
    /*<       SDCD=SD*CD >*/
    c0_2.sdcd = c0_2.sd * c0_2.cd;
    /*<       S2D=F2*SDCD >*/
    c0_2.s2d = f2 * c0_2.sdcd;
    /*<       C2D=CDCD-SDSD >*/
    c0_2.c2d = c0_2.cdcd - c0_2.sdsd;
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* dccon0_ */

/*<       SUBROUTINE  DCCON2(XI,ET,Q,SD,CD,KXI,KET) >*/
/* Subroutine */ int dccon2_(doublereal *xi, doublereal *et, doublereal *q,
                             doublereal *sd, doublereal *cd, integer *kxi,
                             integer *ket) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal eps = 1e-6;

    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal), log(doublereal);

    /* Local variables */
    static doublereal ret, rxi;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/

    /* ********************************************************************** */
    /* *****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** */
    /* ********************************************************************** */

    /* ***** INPUT */
    /* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM */
    /* *****   SD,CD   : SIN, COS OF DIP-ANGLE */
    /* *****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY */

    /* ### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER0
     */

    /*<        >*/
    /*<       DATA  F0,F1,F2,EPS/0.D0,1.D0,2.D0,1.D-6/ >*/
    /* ----- */
    /*<       IF(DABS(XI).LT.EPS) XI=F0 >*/
    if (fabs(*xi) < eps) {
        *xi = f0;
    }
    /*<       IF(DABS(ET).LT.EPS) ET=F0 >*/
    if (fabs(*et) < eps) {
        *et = f0;
    }
    /*<       IF(DABS( Q).LT.EPS)  Q=F0 >*/
    if (fabs(*q) < eps) {
        *q = f0;
    }
    /*<       XI2=XI*XI >*/
    c2_1.xi2 = *xi * *xi;
    /*<       ET2=ET*ET >*/
    c2_1.et2 = *et * *et;
    /*<       Q2=Q*Q >*/
    c2_1.q2 = *q * *q;
    /*<       R2=XI2+ET2+Q2 >*/
    c2_1.r2 = c2_1.xi2 + c2_1.et2 + c2_1.q2;
    /*<       R =DSQRT(R2) >*/
    c2_1.r__ = sqrt(c2_1.r2);
    /*<       IF(R.EQ.F0) RETURN >*/
    if (c2_1.r__ == f0) {
        return 0;
    }
    /*<       R3=R *R2 >*/
    c2_1.r3 = c2_1.r__ * c2_1.r2;
    /*<       R5=R3*R2 >*/
    c2_1.r5 = c2_1.r3 * c2_1.r2;
    /*<       Y =ET*CD+Q*SD >*/
    c2_1.y = *et * *cd + *q * *sd;
    /*<       D =ET*SD-Q*CD >*/
    c2_1.d__ = *et * *sd - *q * *cd;
    /* ----- */
    /*<       IF(Q.EQ.F0) THEN >*/
    if (*q == f0) {
        /*<         TT=F0 >*/
        c2_1.tt = f0;
        /*<       ELSE >*/
    } else {
        /*<         TT=DATAN(XI*ET/(Q*R)) >*/
        c2_1.tt = atan(*xi * *et / (*q * c2_1.r__));
        /*<       ENDIF >*/
    }
    /* ----- */
    /*<       IF(KXI.EQ.1) THEN >*/
    if (*kxi == 1) {
        /*<         ALX=-DLOG(R-XI) >*/
        c2_1.alx = -log(c2_1.r__ - *xi);
        /*<         X11=F0 >*/
        c2_1.x11 = f0;
        /*<         X32=F0 >*/
        c2_1.x32 = f0;
        /*<       ELSE >*/
    } else {
        /*<         RXI=R+XI >*/
        rxi = c2_1.r__ + *xi;
        /*<         ALX=DLOG(RXI) >*/
        c2_1.alx = log(rxi);
        /*<         X11=F1/(R*RXI) >*/
        c2_1.x11 = f1 / (c2_1.r__ * rxi);
        /*<         X32=(R+RXI)*X11*X11/R >*/
        c2_1.x32 = (c2_1.r__ + rxi) * c2_1.x11 * c2_1.x11 / c2_1.r__;
        /*<       ENDIF >*/
    }
    /* ----- */
    /*<       IF(KET.EQ.1) THEN >*/
    if (*ket == 1) {
        /*<         ALE=-DLOG(R-ET) >*/
        c2_1.ale = -log(c2_1.r__ - *et);
        /*<         Y11=F0 >*/
        c2_1.y11 = f0;
        /*<         Y32=F0 >*/
        c2_1.y32 = f0;
        /*<       ELSE >*/
    } else {
        /*<         RET=R+ET >*/
        ret = c2_1.r__ + *et;
        /*<         ALE=DLOG(RET) >*/
        c2_1.ale = log(ret);
        /*<         Y11=F1/(R*RET) >*/
        c2_1.y11 = f1 / (c2_1.r__ * ret);
        /*<         Y32=(R+RET)*Y11*Y11/R >*/
        c2_1.y32 = (c2_1.r__ + ret) * c2_1.y11 * c2_1.y11 / c2_1.r__;
        /*<       ENDIF >*/
    }
    /* ----- */
    /*<       EY=SD/R-Y*Q/R3 >*/
    c2_1.ey = *sd / c2_1.r__ - c2_1.y * *q / c2_1.r3;
    /*<       EZ=CD/R+D*Q/R3 >*/
    c2_1.ez = *cd / c2_1.r__ + c2_1.d__ * *q / c2_1.r3;
    /*<       FY=D/R3+XI2*Y32*SD >*/
    c2_1.fy = c2_1.d__ / c2_1.r3 + c2_1.xi2 * c2_1.y32 * *sd;
    /*<       FZ=Y/R3+XI2*Y32*CD >*/
    c2_1.fz = c2_1.y / c2_1.r3 + c2_1.xi2 * c2_1.y32 * *cd;
    /*<       GY=F2*X11*SD-Y*Q*X32 >*/
    c2_1.gy = f2 * c2_1.x11 * *sd - c2_1.y * *q * c2_1.x32;
    /*<       GZ=F2*X11*CD+D*Q*X32 >*/
    c2_1.gz = f2 * c2_1.x11 * *cd + c2_1.d__ * *q * c2_1.x32;
    /*<       HY=D*Q*X32+XI*Q*Y32*SD >*/
    c2_1.hy = c2_1.d__ * *q * c2_1.x32 + *xi * *q * c2_1.y32 * *sd;
    /*<       HZ=Y*Q*X32+XI*Q*Y32*CD >*/
    c2_1.hz = c2_1.y * *q * c2_1.x32 + *xi * *q * c2_1.y32 * *cd;
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* dccon2_ */

/* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
 */
/* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
 */
