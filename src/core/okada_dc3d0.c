/* ../dc3d0.c -- translated by f2c (version 20230428).
 *
 * **Point Source**
 *
 * The C version of Fortran codes of Okada dislocation (Okada, 1992).
 * Zelong Guo, 29.01.2024
 * @ Potsdam, Germany
 */

// #include "f2c.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Common Block Declarations */

typedef double doublereal;
typedef long int integer;

union {
    struct {
        doublereal dummy[8], r__;
    } _1;
    struct {
        doublereal p, q, s, t, xy, x2, y2, d2, r__, r2, r3, r5, qr, qrx, a3, a5, b3, c3, uy, vy, wy,
            uz, vz, wz;
    } _2;
    struct {
        doublereal p, q, s, t, xy, x2, y2, d2, r__, r2, r3, r5, qr, qrx, a3, a5, b3, c3;
    } _3;
} c1_;

#define c1_1 (c1_._1)
#define c1_2 (c1_._2)
#define c1_3 (c1_._3)

static union {
    struct {
        doublereal alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, c2d;
    } _1;
    struct {
        doublereal dummy[5], sd, cd;
    } _2;
} c0_;

#define c0_1 (c0_._1)
#define c0_2 (c0_._2)

/*<        >*/
/* Subroutine */ int dc3d0_(doublereal *alpha, doublereal *x, doublereal *y, doublereal *z__,
                            doublereal *depth, doublereal *dip, doublereal *pot1, doublereal *pot2,
                            doublereal *pot3, doublereal *pot4, doublereal *ux, doublereal *uy,
                            doublereal *uz, doublereal *uxx, doublereal *uyx, doublereal *uzx,
                            doublereal *uxy, doublereal *uyy, doublereal *uzy, doublereal *uxz,
                            doublereal *uyz, doublereal *uzz, int *iret) {
    /* Initialized data */

    static doublereal f0 = 0.;

    static integer i__;
    static doublereal u[12], dd, du, xx, yy, zz;
    extern /* Subroutine */ int ua0_(doublereal *, doublereal *, doublereal *, doublereal *,
                                     doublereal *, doublereal *, doublereal *, doublereal *),
        ub0_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
             doublereal *, doublereal *, doublereal *),
        uc0_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
             doublereal *, doublereal *, doublereal *);
    static doublereal pp1, pp2, pp3, pp4, dua[12], dub[12], duc[12], ddip;
    extern /* Subroutine */ int dccon00_(doublereal *, doublereal *),
        dccon1_(doublereal *, doublereal *, doublereal *);
    static doublereal aalpha;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
    /*<        >*/

    /* ******************************************************************** */
    /* *****                                                          ***** */
    /* *****    DISPLACEMENT AND STRAIN AT DEPTH                      ***** */
    /* *****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   ***** */
    /* *****                         CODED BY  Y.OKADA ... SEP.1991   ***** */
    /* *****                         REVISED     NOV.1991, MAY.2002   ***** */
    /* *****                                                          ***** */
    /* ******************************************************************** */

    /* ***** INPUT */
    /* *****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU) */
    /* *****   X,Y,Z : COORDINATE OF OBSERVING POINT */
    /* *****   DEPTH : SOURCE DEPTH */
    /* *****   DIP   : DIP-ANGLE (DEGREE) */
    /* *****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY */
    /* *****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2 */
    /* *****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3 */
    /* *****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4 */

    /* ***** OUTPUT */
    /* *****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) / */
    /* *****               :                     (UNIT OF X,Y,Z,DEPTH)**2  ) */
    /* *****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) / */
    /* *****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  ) */
    /* *****   UXZ,UYZ,UZZ : Z-DERIVATIVE */
    /* *****   IRET        : RETURN CODE */
    /* *****               :   =0....NORMAL */
    /* *****               :   =1....SINGULAR */
    /* *****               :   =2....POSITIVE Z WAS GIVEN */

    /*<       COMMON /C1/DUMMY(8),R >*/
    /*<       DIMENSION  U(12),DUA(12),DUB(12),DUC(12) >*/
    /*<       DATA  F0/0.D0/ >*/
    /* ----- */
    /*<       IRET=0 >*/
    *iret = 0;
    /*<       IF(Z.GT.0.) THEN >*/
    if (*z__ > 0.f) {
        /*<         IRET=2 >*/
        *iret = 2;
        /*<         GO TO 99 >*/
        // goto L99;
        printf(
            "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n");
        fprintf(stderr, "Error, positive z is specified!");
        printf(
            "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n");
        exit(EXIT_FAILURE);
        /*<       ENDIF >*/
    }
    /* ----- */
    /*<       DO 111 I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<         U(I)=F0 >*/
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
    dccon00_(&aalpha, &ddip);
    /* ====================================== */
    /* =====  REAL-SOURCE CONTRIBUTION  ===== */
    /* ====================================== */
    /*<       XX=X >*/
    xx = *x;
    /*<       YY=Y >*/
    yy = *y;
    /*<       ZZ=Z >*/
    zz = *z__;
    /*<       DD=DEPTH+Z >*/
    dd = *depth + *z__;
    /*<       CALL DCCON1(XX,YY,DD) >*/
    dccon1_(&xx, &yy, &dd);
    /*<       IF(R.EQ.F0) THEN >*/
    if (c1_1.r__ == f0) {
        /*<         IRET=1 >*/
        *iret = 1;
        /*<         GO TO 99 >*/
        goto L99;
        /*<       ENDIF >*/
    }
    /* ----- */
    /*<       PP1=POT1 >*/
    pp1 = *pot1;
    /*<       PP2=POT2 >*/
    pp2 = *pot2;
    /*<       PP3=POT3 >*/
    pp3 = *pot3;
    /*<       PP4=POT4 >*/
    pp4 = *pot4;
    /*<       CALL UA0(XX,YY,DD,PP1,PP2,PP3,PP4,DUA) >*/
    ua0_(&xx, &yy, &dd, &pp1, &pp2, &pp3, &pp4, dua);
    /* ----- */
    /*<       DO 222 I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<         IF(I.LT.10) U(I)=U(I)-DUA(I) >*/
        if (i__ < 10) {
            u[i__ - 1] -= dua[i__ - 1];
        }
        /*<         IF(I.GE.10) U(I)=U(I)+DUA(I) >*/
        if (i__ >= 10) {
            u[i__ - 1] += dua[i__ - 1];
        }
        /*<   222 CONTINUE >*/
        /* L222: */
    }
    /* ======================================= */
    /* =====  IMAGE-SOURCE CONTRIBUTION  ===== */
    /* ======================================= */
    /*<       DD=DEPTH-Z >*/
    dd = *depth - *z__;
    /*<       CALL DCCON1(XX,YY,DD) >*/
    dccon1_(&xx, &yy, &dd);
    /*<       CALL UA0(XX,YY,DD,PP1,PP2,PP3,PP4,DUA) >*/
    ua0_(&xx, &yy, &dd, &pp1, &pp2, &pp3, &pp4, dua);
    /*<       CALL UB0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,DUB) >*/
    ub0_(&xx, &yy, &dd, &zz, &pp1, &pp2, &pp3, &pp4, dub);
    /*<       CALL UC0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,DUC) >*/
    uc0_(&xx, &yy, &dd, &zz, &pp1, &pp2, &pp3, &pp4, duc);
    /* ----- */
    /*<       DO 333 I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<         DU=DUA(I)+DUB(I)+ZZ*DUC(I) >*/
        du = dua[i__ - 1] + dub[i__ - 1] + zz * duc[i__ - 1];
        /*<         IF(I.GE.10) DU=DU+DUC(I-9) >*/
        if (i__ >= 10) {
            du += duc[i__ - 10];
        }
        /*<         U(I)=U(I)+DU >*/
        u[i__ - 1] += du;
        /*<   333 CONTINUE >*/
        /* L333: */
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
/* ======================================= */
/* =====  IN CASE OF SINGULAR (R=0)  ===== */
/* ======================================= */
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
} /* dc3d0_ */

/*<       SUBROUTINE  UA0(X,Y,D,POT1,POT2,POT3,POT4,U) >*/
/* Subroutine */ int ua0_(doublereal *x, doublereal *y, doublereal *d__, doublereal *pot1,
                          doublereal *pot2, doublereal *pot3, doublereal *pot4, doublereal *u) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f3 = 3.;
    static doublereal pi2 = 6.283185307179586;

    static integer i__;
    static doublereal du[12];

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
    /*<       DIMENSION U(12),DU(12) >*/

    /* ******************************************************************** */
    /* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             ***** */
    /* *****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   ***** */
    /* ******************************************************************** */

    /* ***** INPUT */
    /* *****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM */
    /* *****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY */
    /* ***** OUTPUT */
    /* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /*<       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D >*/
    /*<        >*/
    /*<       DATA F0,F1,F3/0.D0,1.D0,3.D0/ >*/
    /* Parameter adjustments */
    --u;

    /* Function Body */
    /*<       DATA PI2/6.283185307179586D0/ >*/
    /* ----- */
    /*<       DO 111  I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<   111 U(I)=F0 >*/
        /* L111: */
        u[i__] = f0;
    }
    /* ====================================== */
    /* =====  STRIKE-SLIP CONTRIBUTION  ===== */
    /* ====================================== */
    /*<       IF(POT1.NE.F0) THEN >*/
    if (*pot1 != f0) {
        /*<         DU( 1)= ALP1*Q/R3    +ALP2*X2*QR >*/
        du[0] = c0_1.alp1 * c1_2.q / c1_2.r3 + c0_1.alp2 * c1_2.x2 * c1_2.qr;
        /*<         DU( 2)= ALP1*X/R3*SD +ALP2*XY*QR >*/
        du[1] = c0_1.alp1 * *x / c1_2.r3 * c0_1.sd + c0_1.alp2 * c1_2.xy * c1_2.qr;
        /*<         DU( 3)=-ALP1*X/R3*CD +ALP2*X*D*QR >*/
        du[2] = -c0_1.alp1 * *x / c1_2.r3 * c0_1.cd + c0_1.alp2 * *x * *d__ * c1_2.qr;
        /*<         DU( 4)= X*QR*(-ALP1 +ALP2*(F1+A5) ) >*/
        du[3] = *x * c1_2.qr * (-c0_1.alp1 + c0_1.alp2 * (f1 + c1_2.a5));
        /*<         DU( 5)= ALP1*A3/R3*SD +ALP2*Y*QR*A5 >*/
        du[4] = c0_1.alp1 * c1_2.a3 / c1_2.r3 * c0_1.sd + c0_1.alp2 * *y * c1_2.qr * c1_2.a5;
        /*<         DU( 6)=-ALP1*A3/R3*CD +ALP2*D*QR*A5 >*/
        du[5] = -c0_1.alp1 * c1_2.a3 / c1_2.r3 * c0_1.cd + c0_1.alp2 * *d__ * c1_2.qr * c1_2.a5;
        /*<         DU( 7)= ALP1*(SD/R3-Y*QR) +ALP2*F3*X2/R5*UY >*/
        du[6] = c0_1.alp1 * (c0_1.sd / c1_2.r3 - *y * c1_2.qr) +
                c0_1.alp2 * f3 * c1_2.x2 / c1_2.r5 * c1_2.uy;
        /*<         DU( 8)= F3*X/R5*(-ALP1*Y*SD +ALP2*(Y*UY+Q) ) >*/
        du[7] =
            f3 * *x / c1_2.r5 * (-c0_1.alp1 * *y * c0_1.sd + c0_1.alp2 * (*y * c1_2.uy + c1_2.q));
        /*<         DU( 9)= F3*X/R5*( ALP1*Y*CD +ALP2*D*UY ) >*/
        du[8] = f3 * *x / c1_2.r5 * (c0_1.alp1 * *y * c0_1.cd + c0_1.alp2 * *d__ * c1_2.uy);
        /*<         DU(10)= ALP1*(CD/R3+D*QR) +ALP2*F3*X2/R5*UZ >*/
        du[9] = c0_1.alp1 * (c0_1.cd / c1_2.r3 + *d__ * c1_2.qr) +
                c0_1.alp2 * f3 * c1_2.x2 / c1_2.r5 * c1_2.uz;
        /*<         DU(11)= F3*X/R5*( ALP1*D*SD +ALP2*Y*UZ ) >*/
        du[10] = f3 * *x / c1_2.r5 * (c0_1.alp1 * *d__ * c0_1.sd + c0_1.alp2 * *y * c1_2.uz);
        /*<         DU(12)= F3*X/R5*(-ALP1*D*CD +ALP2*(D*UZ-Q) ) >*/
        du[11] = f3 * *x / c1_2.r5 *
                 (-c0_1.alp1 * *d__ * c0_1.cd + c0_1.alp2 * (*d__ * c1_2.uz - c1_2.q));
        /*<         DO 222 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   222   U(I)=U(I)+POT1/PI2*DU(I) >*/
            /* L222: */
            u[i__] += *pot1 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* =================================== */
    /* =====  DIP-SLIP CONTRIBUTION  ===== */
    /* =================================== */
    /*<       IF(POT2.NE.F0) THEN >*/
    if (*pot2 != f0) {
        /*<         DU( 1)=            ALP2*X*P*QR >*/
        du[0] = c0_1.alp2 * *x * c1_2.p * c1_2.qr;
        /*<         DU( 2)= ALP1*S/R3 +ALP2*Y*P*QR >*/
        du[1] = c0_1.alp1 * c1_2.s / c1_2.r3 + c0_1.alp2 * *y * c1_2.p * c1_2.qr;
        /*<         DU( 3)=-ALP1*T/R3 +ALP2*D*P*QR >*/
        du[2] = -c0_1.alp1 * c1_2.t / c1_2.r3 + c0_1.alp2 * *d__ * c1_2.p * c1_2.qr;
        /*<         DU( 4)=                 ALP2*P*QR*A5 >*/
        du[3] = c0_1.alp2 * c1_2.p * c1_2.qr * c1_2.a5;
        /*<         DU( 5)=-ALP1*F3*X*S/R5 -ALP2*Y*P*QRX >*/
        du[4] = -c0_1.alp1 * f3 * *x * c1_2.s / c1_2.r5 - c0_1.alp2 * *y * c1_2.p * c1_2.qrx;
        /*<         DU( 6)= ALP1*F3*X*T/R5 -ALP2*D*P*QRX >*/
        du[5] = c0_1.alp1 * f3 * *x * c1_2.t / c1_2.r5 - c0_1.alp2 * *d__ * c1_2.p * c1_2.qrx;
        /*<         DU( 7)=                          ALP2*F3*X/R5*VY >*/
        du[6] = c0_1.alp2 * f3 * *x / c1_2.r5 * c1_2.vy;
        /*<         DU( 8)= ALP1*(S2D/R3-F3*Y*S/R5) +ALP2*(F3*Y/R5*VY+P*QR) >*/
        du[7] = c0_1.alp1 * (c0_1.s2d / c1_2.r3 - f3 * *y * c1_2.s / c1_2.r5) +
                c0_1.alp2 * (f3 * *y / c1_2.r5 * c1_2.vy + c1_2.p * c1_2.qr);
        /*<         DU( 9)=-ALP1*(C2D/R3-F3*Y*T/R5) +ALP2*F3*D/R5*VY >*/
        du[8] = -c0_1.alp1 * (c0_1.c2d / c1_2.r3 - f3 * *y * c1_2.t / c1_2.r5) +
                c0_1.alp2 * f3 * *d__ / c1_2.r5 * c1_2.vy;
        /*<         DU(10)=                          ALP2*F3*X/R5*VZ >*/
        du[9] = c0_1.alp2 * f3 * *x / c1_2.r5 * c1_2.vz;
        /*<         DU(11)= ALP1*(C2D/R3+F3*D*S/R5) +ALP2*F3*Y/R5*VZ >*/
        du[10] = c0_1.alp1 * (c0_1.c2d / c1_2.r3 + f3 * *d__ * c1_2.s / c1_2.r5) +
                 c0_1.alp2 * f3 * *y / c1_2.r5 * c1_2.vz;
        /*<         DU(12)= ALP1*(S2D/R3-F3*D*T/R5) +ALP2*(F3*D/R5*VZ-P*QR) >*/
        du[11] = c0_1.alp1 * (c0_1.s2d / c1_2.r3 - f3 * *d__ * c1_2.t / c1_2.r5) +
                 c0_1.alp2 * (f3 * *d__ / c1_2.r5 * c1_2.vz - c1_2.p * c1_2.qr);
        /*<         DO 333 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   333   U(I)=U(I)+POT2/PI2*DU(I) >*/
            /* L333: */
            u[i__] += *pot2 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ======================================== */
    /* =====  TENSILE-FAULT CONTRIBUTION  ===== */
    /* ======================================== */
    /*<       IF(POT3.NE.F0) THEN >*/
    if (*pot3 != f0) {
        /*<         DU( 1)= ALP1*X/R3 -ALP2*X*Q*QR >*/
        du[0] = c0_1.alp1 * *x / c1_2.r3 - c0_1.alp2 * *x * c1_2.q * c1_2.qr;
        /*<         DU( 2)= ALP1*T/R3 -ALP2*Y*Q*QR >*/
        du[1] = c0_1.alp1 * c1_2.t / c1_2.r3 - c0_1.alp2 * *y * c1_2.q * c1_2.qr;
        /*<         DU( 3)= ALP1*S/R3 -ALP2*D*Q*QR >*/
        du[2] = c0_1.alp1 * c1_2.s / c1_2.r3 - c0_1.alp2 * *d__ * c1_2.q * c1_2.qr;
        /*<         DU( 4)= ALP1*A3/R3     -ALP2*Q*QR*A5 >*/
        du[3] = c0_1.alp1 * c1_2.a3 / c1_2.r3 - c0_1.alp2 * c1_2.q * c1_2.qr * c1_2.a5;
        /*<         DU( 5)=-ALP1*F3*X*T/R5 +ALP2*Y*Q*QRX >*/
        du[4] = -c0_1.alp1 * f3 * *x * c1_2.t / c1_2.r5 + c0_1.alp2 * *y * c1_2.q * c1_2.qrx;
        /*<         DU( 6)=-ALP1*F3*X*S/R5 +ALP2*D*Q*QRX >*/
        du[5] = -c0_1.alp1 * f3 * *x * c1_2.s / c1_2.r5 + c0_1.alp2 * *d__ * c1_2.q * c1_2.qrx;
        /*<         DU( 7)=-ALP1*F3*XY/R5           -ALP2*X*QR*WY >*/
        du[6] = -c0_1.alp1 * f3 * c1_2.xy / c1_2.r5 - c0_1.alp2 * *x * c1_2.qr * c1_2.wy;
        /*<         DU( 8)= ALP1*(C2D/R3-F3*Y*T/R5) -ALP2*(Y*WY+Q)*QR >*/
        du[7] = c0_1.alp1 * (c0_1.c2d / c1_2.r3 - f3 * *y * c1_2.t / c1_2.r5) -
                c0_1.alp2 * (*y * c1_2.wy + c1_2.q) * c1_2.qr;
        /*<         DU( 9)= ALP1*(S2D/R3-F3*Y*S/R5) -ALP2*D*QR*WY >*/
        du[8] = c0_1.alp1 * (c0_1.s2d / c1_2.r3 - f3 * *y * c1_2.s / c1_2.r5) -
                c0_1.alp2 * *d__ * c1_2.qr * c1_2.wy;
        /*<         DU(10)= ALP1*F3*X*D/R5          -ALP2*X*QR*WZ >*/
        du[9] = c0_1.alp1 * f3 * *x * *d__ / c1_2.r5 - c0_1.alp2 * *x * c1_2.qr * c1_2.wz;
        /*<         DU(11)=-ALP1*(S2D/R3-F3*D*T/R5) -ALP2*Y*QR*WZ >*/
        du[10] = -c0_1.alp1 * (c0_1.s2d / c1_2.r3 - f3 * *d__ * c1_2.t / c1_2.r5) -
                 c0_1.alp2 * *y * c1_2.qr * c1_2.wz;
        /*<         DU(12)= ALP1*(C2D/R3+F3*D*S/R5) -ALP2*(D*WZ-Q)*QR >*/
        du[11] = c0_1.alp1 * (c0_1.c2d / c1_2.r3 + f3 * *d__ * c1_2.s / c1_2.r5) -
                 c0_1.alp2 * (*d__ * c1_2.wz - c1_2.q) * c1_2.qr;
        /*<         DO 444 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   444   U(I)=U(I)+POT3/PI2*DU(I) >*/
            /* L444: */
            u[i__] += *pot3 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ========================================= */
    /* =====  INFLATE SOURCE CONTRIBUTION  ===== */
    /* ========================================= */
    /*<       IF(POT4.NE.F0) THEN >*/
    if (*pot4 != f0) {
        /*<         DU( 1)=-ALP1*X/R3 >*/
        du[0] = -c0_1.alp1 * *x / c1_2.r3;
        /*<         DU( 2)=-ALP1*Y/R3 >*/
        du[1] = -c0_1.alp1 * *y / c1_2.r3;
        /*<         DU( 3)=-ALP1*D/R3 >*/
        du[2] = -c0_1.alp1 * *d__ / c1_2.r3;
        /*<         DU( 4)=-ALP1*A3/R3 >*/
        du[3] = -c0_1.alp1 * c1_2.a3 / c1_2.r3;
        /*<         DU( 5)= ALP1*F3*XY/R5 >*/
        du[4] = c0_1.alp1 * f3 * c1_2.xy / c1_2.r5;
        /*<         DU( 6)= ALP1*F3*X*D/R5 >*/
        du[5] = c0_1.alp1 * f3 * *x * *d__ / c1_2.r5;
        /*<         DU( 7)= DU(5) >*/
        du[6] = du[4];
        /*<         DU( 8)=-ALP1*B3/R3 >*/
        du[7] = -c0_1.alp1 * c1_2.b3 / c1_2.r3;
        /*<         DU( 9)= ALP1*F3*Y*D/R5 >*/
        du[8] = c0_1.alp1 * f3 * *y * *d__ / c1_2.r5;
        /*<         DU(10)=-DU(6) >*/
        du[9] = -du[5];
        /*<         DU(11)=-DU(9) >*/
        du[10] = -du[8];
        /*<         DU(12)= ALP1*C3/R3 >*/
        du[11] = c0_1.alp1 * c1_2.c3 / c1_2.r3;
        /*<         DO 555 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   555   U(I)=U(I)+POT4/PI2*DU(I) >*/
            /* L555: */
            u[i__] += *pot4 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* ua0_ */

/*<       SUBROUTINE  UB0(X,Y,D,Z,POT1,POT2,POT3,POT4,U) >*/
/* Subroutine */ int ub0_(doublereal *x, doublereal *y, doublereal *d__, doublereal *z__,
                          doublereal *pot1, doublereal *pot2, doublereal *pot3, doublereal *pot4,
                          doublereal *u) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal f3 = 3.;
    static doublereal f4 = 4.;
    static doublereal f5 = 5.;
    static doublereal f8 = 8.;
    static doublereal f9 = 9.;
    static doublereal pi2 = 6.283185307179586;

    static doublereal c__;
    static integer i__;
    static doublereal d12, d32, d33, d53, d54, rd, du[12], fi1, fi2, fi3, fi4, fi5, fj1, fj2, fj3,
        fj4, fk1, fk2, fk3;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
    /*<       DIMENSION U(12),DU(12) >*/

    /* ******************************************************************** */
    /* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             ***** */
    /* *****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   ***** */
    /* ******************************************************************** */

    /* ***** INPUT */
    /* *****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM */
    /* *****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY */
    /* ***** OUTPUT */
    /* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /*<       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D >*/
    /*<        >*/
    /*<        >*/
    /* Parameter adjustments */
    --u;

    /* Function Body */
    /*<       DATA PI2/6.283185307179586D0/ >*/
    /* ----- */
    /*<       C=D+Z >*/
    c__ = *d__ + *z__;
    /*<       RD=R+D >*/
    rd = c1_2.r__ + *d__;
    /*<       D12=F1/(R*RD*RD) >*/
    d12 = f1 / (c1_2.r__ * rd * rd);
    /*<       D32=D12*(F2*R+D)/R2 >*/
    d32 = d12 * (f2 * c1_2.r__ + *d__) / c1_2.r2;
    /*<       D33=D12*(F3*R+D)/(R2*RD) >*/
    d33 = d12 * (f3 * c1_2.r__ + *d__) / (c1_2.r2 * rd);
    /*<       D53=D12*(F8*R2+F9*R*D+F3*D2)/(R2*R2*RD) >*/
    d53 = d12 * (f8 * c1_2.r2 + f9 * c1_2.r__ * *d__ + f3 * c1_2.d2) / (c1_2.r2 * c1_2.r2 * rd);
    /*<       D54=D12*(F5*R2+F4*R*D+D2)/R3*D12 >*/
    d54 = d12 * (f5 * c1_2.r2 + f4 * c1_2.r__ * *d__ + c1_2.d2) / c1_2.r3 * d12;
    /* ----- */
    /*<       FI1= Y*(D12-X2*D33) >*/
    fi1 = *y * (d12 - c1_2.x2 * d33);
    /*<       FI2= X*(D12-Y2*D33) >*/
    fi2 = *x * (d12 - c1_2.y2 * d33);
    /*<       FI3= X/R3-FI2 >*/
    fi3 = *x / c1_2.r3 - fi2;
    /*<       FI4=-XY*D32 >*/
    fi4 = -c1_2.xy * d32;
    /*<       FI5= F1/(R*RD)-X2*D32 >*/
    fi5 = f1 / (c1_2.r__ * rd) - c1_2.x2 * d32;
    /*<       FJ1=-F3*XY*(D33-X2*D54) >*/
    fj1 = -f3 * c1_2.xy * (d33 - c1_2.x2 * d54);
    /*<       FJ2= F1/R3-F3*D12+F3*X2*Y2*D54 >*/
    fj2 = f1 / c1_2.r3 - f3 * d12 + f3 * c1_2.x2 * c1_2.y2 * d54;
    /*<       FJ3= A3/R3-FJ2 >*/
    fj3 = c1_2.a3 / c1_2.r3 - fj2;
    /*<       FJ4=-F3*XY/R5-FJ1 >*/
    fj4 = -f3 * c1_2.xy / c1_2.r5 - fj1;
    /*<       FK1=-Y*(D32-X2*D53) >*/
    fk1 = -(*y) * (d32 - c1_2.x2 * d53);
    /*<       FK2=-X*(D32-Y2*D53) >*/
    fk2 = -(*x) * (d32 - c1_2.y2 * d53);
    /*<       FK3=-F3*X*D/R5-FK2 >*/
    fk3 = -f3 * *x * *d__ / c1_2.r5 - fk2;
    /* ----- */
    /*<       DO 111  I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<   111 U(I)=F0 >*/
        /* L111: */
        u[i__] = f0;
    }
    /* ====================================== */
    /* =====  STRIKE-SLIP CONTRIBUTION  ===== */
    /* ====================================== */
    /*<       IF(POT1.NE.F0) THEN >*/
    if (*pot1 != f0) {
        /*<         DU( 1)=-X2*QR  -ALP3*FI1*SD >*/
        du[0] = -c1_2.x2 * c1_2.qr - c0_1.alp3 * fi1 * c0_1.sd;
        /*<         DU( 2)=-XY*QR  -ALP3*FI2*SD >*/
        du[1] = -c1_2.xy * c1_2.qr - c0_1.alp3 * fi2 * c0_1.sd;
        /*<         DU( 3)=-C*X*QR -ALP3*FI4*SD >*/
        du[2] = -c__ * *x * c1_2.qr - c0_1.alp3 * fi4 * c0_1.sd;
        /*<         DU( 4)=-X*QR*(F1+A5) -ALP3*FJ1*SD >*/
        du[3] = -(*x) * c1_2.qr * (f1 + c1_2.a5) - c0_1.alp3 * fj1 * c0_1.sd;
        /*<         DU( 5)=-Y*QR*A5      -ALP3*FJ2*SD >*/
        du[4] = -(*y) * c1_2.qr * c1_2.a5 - c0_1.alp3 * fj2 * c0_1.sd;
        /*<         DU( 6)=-C*QR*A5      -ALP3*FK1*SD >*/
        du[5] = -c__ * c1_2.qr * c1_2.a5 - c0_1.alp3 * fk1 * c0_1.sd;
        /*<         DU( 7)=-F3*X2/R5*UY      -ALP3*FJ2*SD >*/
        du[6] = -f3 * c1_2.x2 / c1_2.r5 * c1_2.uy - c0_1.alp3 * fj2 * c0_1.sd;
        /*<         DU( 8)=-F3*XY/R5*UY-X*QR -ALP3*FJ4*SD >*/
        du[7] = -f3 * c1_2.xy / c1_2.r5 * c1_2.uy - *x * c1_2.qr - c0_1.alp3 * fj4 * c0_1.sd;
        /*<         DU( 9)=-F3*C*X/R5*UY     -ALP3*FK2*SD >*/
        du[8] = -f3 * c__ * *x / c1_2.r5 * c1_2.uy - c0_1.alp3 * fk2 * c0_1.sd;
        /*<         DU(10)=-F3*X2/R5*UZ  +ALP3*FK1*SD >*/
        du[9] = -f3 * c1_2.x2 / c1_2.r5 * c1_2.uz + c0_1.alp3 * fk1 * c0_1.sd;
        /*<         DU(11)=-F3*XY/R5*UZ  +ALP3*FK2*SD >*/
        du[10] = -f3 * c1_2.xy / c1_2.r5 * c1_2.uz + c0_1.alp3 * fk2 * c0_1.sd;
        /*<         DU(12)= F3*X/R5*(-C*UZ +ALP3*Y*SD) >*/
        du[11] = f3 * *x / c1_2.r5 * (-c__ * c1_2.uz + c0_1.alp3 * *y * c0_1.sd);
        /*<         DO 222 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   222   U(I)=U(I)+POT1/PI2*DU(I) >*/
            /* L222: */
            u[i__] += *pot1 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* =================================== */
    /* =====  DIP-SLIP CONTRIBUTION  ===== */
    /* =================================== */
    /*<       IF(POT2.NE.F0) THEN >*/
    if (*pot2 != f0) {
        /*<         DU( 1)=-X*P*QR +ALP3*FI3*SDCD >*/
        du[0] = -(*x) * c1_2.p * c1_2.qr + c0_1.alp3 * fi3 * c0_1.sdcd;
        /*<         DU( 2)=-Y*P*QR +ALP3*FI1*SDCD >*/
        du[1] = -(*y) * c1_2.p * c1_2.qr + c0_1.alp3 * fi1 * c0_1.sdcd;
        /*<         DU( 3)=-C*P*QR +ALP3*FI5*SDCD >*/
        du[2] = -c__ * c1_2.p * c1_2.qr + c0_1.alp3 * fi5 * c0_1.sdcd;
        /*<         DU( 4)=-P*QR*A5 +ALP3*FJ3*SDCD >*/
        du[3] = -c1_2.p * c1_2.qr * c1_2.a5 + c0_1.alp3 * fj3 * c0_1.sdcd;
        /*<         DU( 5)= Y*P*QRX +ALP3*FJ1*SDCD >*/
        du[4] = *y * c1_2.p * c1_2.qrx + c0_1.alp3 * fj1 * c0_1.sdcd;
        /*<         DU( 6)= C*P*QRX +ALP3*FK3*SDCD >*/
        du[5] = c__ * c1_2.p * c1_2.qrx + c0_1.alp3 * fk3 * c0_1.sdcd;
        /*<         DU( 7)=-F3*X/R5*VY      +ALP3*FJ1*SDCD >*/
        du[6] = -f3 * *x / c1_2.r5 * c1_2.vy + c0_1.alp3 * fj1 * c0_1.sdcd;
        /*<         DU( 8)=-F3*Y/R5*VY-P*QR +ALP3*FJ2*SDCD >*/
        du[7] = -f3 * *y / c1_2.r5 * c1_2.vy - c1_2.p * c1_2.qr + c0_1.alp3 * fj2 * c0_1.sdcd;
        /*<         DU( 9)=-F3*C/R5*VY      +ALP3*FK1*SDCD >*/
        du[8] = -f3 * c__ / c1_2.r5 * c1_2.vy + c0_1.alp3 * fk1 * c0_1.sdcd;
        /*<         DU(10)=-F3*X/R5*VZ -ALP3*FK3*SDCD >*/
        du[9] = -f3 * *x / c1_2.r5 * c1_2.vz - c0_1.alp3 * fk3 * c0_1.sdcd;
        /*<         DU(11)=-F3*Y/R5*VZ -ALP3*FK1*SDCD >*/
        du[10] = -f3 * *y / c1_2.r5 * c1_2.vz - c0_1.alp3 * fk1 * c0_1.sdcd;
        /*<         DU(12)=-F3*C/R5*VZ +ALP3*A3/R3*SDCD >*/
        du[11] = -f3 * c__ / c1_2.r5 * c1_2.vz + c0_1.alp3 * c1_2.a3 / c1_2.r3 * c0_1.sdcd;
        /*<         DO 333 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   333   U(I)=U(I)+POT2/PI2*DU(I) >*/
            /* L333: */
            u[i__] += *pot2 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ======================================== */
    /* =====  TENSILE-FAULT CONTRIBUTION  ===== */
    /* ======================================== */
    /*<       IF(POT3.NE.F0) THEN >*/
    if (*pot3 != f0) {
        /*<         DU( 1)= X*Q*QR -ALP3*FI3*SDSD >*/
        du[0] = *x * c1_2.q * c1_2.qr - c0_1.alp3 * fi3 * c0_1.sdsd;
        /*<         DU( 2)= Y*Q*QR -ALP3*FI1*SDSD >*/
        du[1] = *y * c1_2.q * c1_2.qr - c0_1.alp3 * fi1 * c0_1.sdsd;
        /*<         DU( 3)= C*Q*QR -ALP3*FI5*SDSD >*/
        du[2] = c__ * c1_2.q * c1_2.qr - c0_1.alp3 * fi5 * c0_1.sdsd;
        /*<         DU( 4)= Q*QR*A5 -ALP3*FJ3*SDSD >*/
        du[3] = c1_2.q * c1_2.qr * c1_2.a5 - c0_1.alp3 * fj3 * c0_1.sdsd;
        /*<         DU( 5)=-Y*Q*QRX -ALP3*FJ1*SDSD >*/
        du[4] = -(*y) * c1_2.q * c1_2.qrx - c0_1.alp3 * fj1 * c0_1.sdsd;
        /*<         DU( 6)=-C*Q*QRX -ALP3*FK3*SDSD >*/
        du[5] = -c__ * c1_2.q * c1_2.qrx - c0_1.alp3 * fk3 * c0_1.sdsd;
        /*<         DU( 7)= X*QR*WY     -ALP3*FJ1*SDSD >*/
        du[6] = *x * c1_2.qr * c1_2.wy - c0_1.alp3 * fj1 * c0_1.sdsd;
        /*<         DU( 8)= QR*(Y*WY+Q) -ALP3*FJ2*SDSD >*/
        du[7] = c1_2.qr * (*y * c1_2.wy + c1_2.q) - c0_1.alp3 * fj2 * c0_1.sdsd;
        /*<         DU( 9)= C*QR*WY     -ALP3*FK1*SDSD >*/
        du[8] = c__ * c1_2.qr * c1_2.wy - c0_1.alp3 * fk1 * c0_1.sdsd;
        /*<         DU(10)= X*QR*WZ +ALP3*FK3*SDSD >*/
        du[9] = *x * c1_2.qr * c1_2.wz + c0_1.alp3 * fk3 * c0_1.sdsd;
        /*<         DU(11)= Y*QR*WZ +ALP3*FK1*SDSD >*/
        du[10] = *y * c1_2.qr * c1_2.wz + c0_1.alp3 * fk1 * c0_1.sdsd;
        /*<         DU(12)= C*QR*WZ -ALP3*A3/R3*SDSD >*/
        du[11] = c__ * c1_2.qr * c1_2.wz - c0_1.alp3 * c1_2.a3 / c1_2.r3 * c0_1.sdsd;
        /*<         DO 444 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   444   U(I)=U(I)+POT3/PI2*DU(I) >*/
            /* L444: */
            u[i__] += *pot3 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ========================================= */
    /* =====  INFLATE SOURCE CONTRIBUTION  ===== */
    /* ========================================= */
    /*<       IF(POT4.NE.F0) THEN >*/
    if (*pot4 != f0) {
        /*<         DU( 1)= ALP3*X/R3 >*/
        du[0] = c0_1.alp3 * *x / c1_2.r3;
        /*<         DU( 2)= ALP3*Y/R3 >*/
        du[1] = c0_1.alp3 * *y / c1_2.r3;
        /*<         DU( 3)= ALP3*D/R3 >*/
        du[2] = c0_1.alp3 * *d__ / c1_2.r3;
        /*<         DU( 4)= ALP3*A3/R3 >*/
        du[3] = c0_1.alp3 * c1_2.a3 / c1_2.r3;
        /*<         DU( 5)=-ALP3*F3*XY/R5 >*/
        du[4] = -c0_1.alp3 * f3 * c1_2.xy / c1_2.r5;
        /*<         DU( 6)=-ALP3*F3*X*D/R5 >*/
        du[5] = -c0_1.alp3 * f3 * *x * *d__ / c1_2.r5;
        /*<         DU( 7)= DU(5) >*/
        du[6] = du[4];
        /*<         DU( 8)= ALP3*B3/R3 >*/
        du[7] = c0_1.alp3 * c1_2.b3 / c1_2.r3;
        /*<         DU( 9)=-ALP3*F3*Y*D/R5 >*/
        du[8] = -c0_1.alp3 * f3 * *y * *d__ / c1_2.r5;
        /*<         DU(10)=-DU(6) >*/
        du[9] = -du[5];
        /*<         DU(11)=-DU(9) >*/
        du[10] = -du[8];
        /*<         DU(12)=-ALP3*C3/R3 >*/
        du[11] = -c0_1.alp3 * c1_2.c3 / c1_2.r3;
        /*<         DO 555 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   555   U(I)=U(I)+POT4/PI2*DU(I) >*/
            /* L555: */
            u[i__] += *pot4 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* ub0_ */

/*<       SUBROUTINE  UC0(X,Y,D,Z,POT1,POT2,POT3,POT4,U) >*/
/* Subroutine */ int uc0_(doublereal *x, doublereal *y, doublereal *d__, doublereal *z__,
                          doublereal *pot1, doublereal *pot2, doublereal *pot3, doublereal *pot4,
                          doublereal *u) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal f3 = 3.;
    static doublereal f5 = 5.;
    static doublereal f7 = 7.;
    static doublereal f10 = 10.;
    static doublereal f15 = 15.;
    static doublereal pi2 = 6.283185307179586;

    static doublereal c__;
    static integer i__;
    static doublereal a7, b5, b7, c5, q2, c7, d7, r7, du[12], dr5, qr5, qr7;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/
    /*<       DIMENSION U(12),DU(12) >*/

    /* ******************************************************************** */
    /* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             ***** */
    /* *****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   ***** */
    /* ******************************************************************** */

    /* ***** INPUT */
    /* *****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM */
    /* *****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY */
    /* ***** OUTPUT */
    /* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES */

    /*<       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D >*/
    /*<       COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3 >*/
    /*<        >*/
    /* Parameter adjustments */
    --u;

    /* Function Body */
    /*<       DATA PI2/6.283185307179586D0/ >*/
    /* ----- */
    /*<       C=D+Z >*/
    c__ = *d__ + *z__;
    /*<       Q2=Q*Q >*/
    q2 = c1_3.q * c1_3.q;
    /*<       R7=R5*R2 >*/
    r7 = c1_3.r5 * c1_3.r2;
    /*<       A7=F1-F7*X2/R2 >*/
    a7 = f1 - f7 * c1_3.x2 / c1_3.r2;
    /*<       B5=F1-F5*Y2/R2 >*/
    b5 = f1 - f5 * c1_3.y2 / c1_3.r2;
    /*<       B7=F1-F7*Y2/R2 >*/
    b7 = f1 - f7 * c1_3.y2 / c1_3.r2;
    /*<       C5=F1-F5*D2/R2 >*/
    c5 = f1 - f5 * c1_3.d2 / c1_3.r2;
    /*<       C7=F1-F7*D2/R2 >*/
    c7 = f1 - f7 * c1_3.d2 / c1_3.r2;
    /*<       D7=F2-F7*Q2/R2 >*/
    d7 = f2 - f7 * q2 / c1_3.r2;
    /*<       QR5=F5*Q/R2 >*/
    qr5 = f5 * c1_3.q / c1_3.r2;
    /*<       QR7=F7*Q/R2 >*/
    qr7 = f7 * c1_3.q / c1_3.r2;
    /*<       DR5=F5*D/R2 >*/
    dr5 = f5 * *d__ / c1_3.r2;
    /* ----- */
    /*<       DO 111  I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
        /*<   111 U(I)=F0 >*/
        /* L111: */
        u[i__] = f0;
    }
    /* ====================================== */
    /* =====  STRIKE-SLIP CONTRIBUTION  ===== */
    /* ====================================== */
    /*<       IF(POT1.NE.F0) THEN >*/
    if (*pot1 != f0) {
        /*<         DU( 1)=-ALP4*A3/R3*CD  +ALP5*C*QR*A5 >*/
        du[0] = -c0_1.alp4 * c1_3.a3 / c1_3.r3 * c0_1.cd + c0_1.alp5 * c__ * c1_3.qr * c1_3.a5;
        /*<         DU( 2)= F3*X/R5*( ALP4*Y*CD +ALP5*C*(SD-Y*QR5) ) >*/
        du[1] =
            f3 * *x / c1_3.r5 * (c0_1.alp4 * *y * c0_1.cd + c0_1.alp5 * c__ * (c0_1.sd - *y * qr5));
        /*<         DU( 3)= F3*X/R5*(-ALP4*Y*SD +ALP5*C*(CD+D*QR5) ) >*/
        du[2] = f3 * *x / c1_3.r5 *
                (-c0_1.alp4 * *y * c0_1.sd + c0_1.alp5 * c__ * (c0_1.cd + *d__ * qr5));
        /*<         DU( 4)= ALP4*F3*X/R5*(F2+A5)*CD   -ALP5*C*QRX*(F2+A7) >*/
        du[3] = c0_1.alp4 * f3 * *x / c1_3.r5 * (f2 + c1_3.a5) * c0_1.cd -
                c0_1.alp5 * c__ * c1_3.qrx * (f2 + a7);
        /*<         DU( 5)= F3/R5*( ALP4*Y*A5*CD +ALP5*C*(A5*SD-Y*QR5*A7) ) >*/
        du[4] = f3 / c1_3.r5 *
                (c0_1.alp4 * *y * c1_3.a5 * c0_1.cd +
                 c0_1.alp5 * c__ * (c1_3.a5 * c0_1.sd - *y * qr5 * a7));
        /*<         DU( 6)= F3/R5*(-ALP4*Y*A5*SD +ALP5*C*(A5*CD+D*QR5*A7) ) >*/
        du[5] = f3 / c1_3.r5 *
                (-c0_1.alp4 * *y * c1_3.a5 * c0_1.sd +
                 c0_1.alp5 * c__ * (c1_3.a5 * c0_1.cd + *d__ * qr5 * a7));
        /*<         DU( 7)= DU(5) >*/
        du[6] = du[4];
        /*<         DU( 8)= F3*X/R5*( ALP4*B5*CD -ALP5*F5*C/R2*(F2*Y*SD+Q*B7) ) >*/
        du[7] = f3 * *x / c1_3.r5 *
                (c0_1.alp4 * b5 * c0_1.cd -
                 c0_1.alp5 * f5 * c__ / c1_3.r2 * (f2 * *y * c0_1.sd + c1_3.q * b7));
        /*<         DU( 9)= F3*X/R5*(-ALP4*B5*SD +ALP5*F5*C/R2*(D*B7*SD-Y*C7*CD) ) >*/
        du[8] = f3 * *x / c1_3.r5 *
                (-c0_1.alp4 * b5 * c0_1.sd +
                 c0_1.alp5 * f5 * c__ / c1_3.r2 * (*d__ * b7 * c0_1.sd - *y * c7 * c0_1.cd));
        /*<         DU(10)= F3/R5*   (-ALP4*D*A5*CD +ALP5*C*(A5*CD+D*QR5*A7) ) >*/
        du[9] = f3 / c1_3.r5 *
                (-c0_1.alp4 * *d__ * c1_3.a5 * c0_1.cd +
                 c0_1.alp5 * c__ * (c1_3.a5 * c0_1.cd + *d__ * qr5 * a7));
        /*<         DU(11)= F15*X/R7*( ALP4*Y*D*CD  +ALP5*C*(D*B7*SD-Y*C7*CD) ) >*/
        du[10] = f15 * *x / r7 *
                 (c0_1.alp4 * *y * *d__ * c0_1.cd +
                  c0_1.alp5 * c__ * (*d__ * b7 * c0_1.sd - *y * c7 * c0_1.cd));
        /*<         DU(12)= F15*X/R7*(-ALP4*Y*D*SD  +ALP5*C*(F2*D*CD-Q*C7) ) >*/
        du[11] = f15 * *x / r7 *
                 (-c0_1.alp4 * *y * *d__ * c0_1.sd +
                  c0_1.alp5 * c__ * (f2 * *d__ * c0_1.cd - c1_3.q * c7));
        /*<         DO 222 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   222   U(I)=U(I)+POT1/PI2*DU(I) >*/
            /* L222: */
            u[i__] += *pot1 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* =================================== */
    /* =====  DIP-SLIP CONTRIBUTION  ===== */
    /* =================================== */
    /*<       IF(POT2.NE.F0) THEN >*/
    if (*pot2 != f0) {
        /*<         DU( 1)= ALP4*F3*X*T/R5          -ALP5*C*P*QRX >*/
        du[0] = c0_1.alp4 * f3 * *x * c1_3.t / c1_3.r5 - c0_1.alp5 * c__ * c1_3.p * c1_3.qrx;
        /*<         DU( 2)=-ALP4/R3*(C2D-F3*Y*T/R2) +ALP5*F3*C/R5*(S-Y*P*QR5) >*/
        du[1] = -c0_1.alp4 / c1_3.r3 * (c0_1.c2d - f3 * *y * c1_3.t / c1_3.r2) +
                c0_1.alp5 * f3 * c__ / c1_3.r5 * (c1_3.s - *y * c1_3.p * qr5);
        /*<         DU( 3)=-ALP4*A3/R3*SDCD         +ALP5*F3*C/R5*(T+D*P*QR5) >*/
        du[2] = -c0_1.alp4 * c1_3.a3 / c1_3.r3 * c0_1.sdcd +
                c0_1.alp5 * f3 * c__ / c1_3.r5 * (c1_3.t + *d__ * c1_3.p * qr5);
        /*<         DU( 4)= ALP4*F3*T/R5*A5              -ALP5*F5*C*P*QR/R2*A7 >*/
        du[3] = c0_1.alp4 * f3 * c1_3.t / c1_3.r5 * c1_3.a5 -
                c0_1.alp5 * f5 * c__ * c1_3.p * c1_3.qr / c1_3.r2 * a7;
        /*<         DU( 5)= F3*X/R5*(ALP4*(C2D-F5*Y*T/R2)-ALP5*F5*C/R2*(S-Y*P*QR7)) >*/
        du[4] = f3 * *x / c1_3.r5 *
                (c0_1.alp4 * (c0_1.c2d - f5 * *y * c1_3.t / c1_3.r2) -
                 c0_1.alp5 * f5 * c__ / c1_3.r2 * (c1_3.s - *y * c1_3.p * qr7));
        /*<         DU( 6)= F3*X/R5*(ALP4*(F2+A5)*SDCD   -ALP5*F5*C/R2*(T+D*P*QR7)) >*/
        du[5] = f3 * *x / c1_3.r5 *
                (c0_1.alp4 * (f2 + c1_3.a5) * c0_1.sdcd -
                 c0_1.alp5 * f5 * c__ / c1_3.r2 * (c1_3.t + *d__ * c1_3.p * qr7));
        /*<         DU( 7)= DU(5) >*/
        du[6] = du[4];
        /*<        >*/
        du[7] = f3 / c1_3.r5 *
                (c0_1.alp4 * (f2 * *y * c0_1.c2d + c1_3.t * b5) +
                 c0_1.alp5 * c__ * (c0_1.s2d - f10 * *y * c1_3.s / c1_3.r2 - c1_3.p * qr5 * b7));
        /*<         DU( 9)= F3/R5*(ALP4*Y*A5*SDCD-ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7)) >*/
        du[8] = f3 / c1_3.r5 *
                (c0_1.alp4 * *y * c1_3.a5 * c0_1.sdcd -
                 c0_1.alp5 * c__ * ((f3 + c1_3.a5) * c0_1.c2d + *y * c1_3.p * dr5 * qr7));
        /*<         DU(10)= F3*X/R5*(-ALP4*(S2D-T*DR5) -ALP5*F5*C/R2*(T+D*P*QR7)) >*/
        du[9] = f3 * *x / c1_3.r5 *
                (-c0_1.alp4 * (c0_1.s2d - c1_3.t * dr5) -
                 c0_1.alp5 * f5 * c__ / c1_3.r2 * (c1_3.t + *d__ * c1_3.p * qr7));
        /*<        >*/
        du[10] = f3 / c1_3.r5 *
                 (-c0_1.alp4 * (*d__ * b5 * c0_1.c2d + *y * c5 * c0_1.s2d) -
                  c0_1.alp5 * c__ * ((f3 + c1_3.a5) * c0_1.c2d + *y * c1_3.p * dr5 * qr7));
        /*<         DU(12)= F3/R5*(-ALP4*D*A5*SDCD-ALP5*C*(S2D-F10*D*T/R2+P*QR5*C7)) >*/
        du[11] = f3 / c1_3.r5 *
                 (-c0_1.alp4 * *d__ * c1_3.a5 * c0_1.sdcd -
                  c0_1.alp5 * c__ * (c0_1.s2d - f10 * *d__ * c1_3.t / c1_3.r2 + c1_3.p * qr5 * c7));
        /*<         DO 333 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   333   U(I)=U(I)+POT2/PI2*DU(I) >*/
            /* L333: */
            u[i__] += *pot2 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ======================================== */
    /* =====  TENSILE-FAULT CONTRIBUTION  ===== */
    /* ======================================== */
    /*<       IF(POT3.NE.F0) THEN >*/
    if (*pot3 != f0) {
        /*<         DU( 1)= F3*X/R5*(-ALP4*S +ALP5*(C*Q*QR5-Z)) >*/
        du[0] = f3 * *x / c1_3.r5 * (-c0_1.alp4 * c1_3.s + c0_1.alp5 * (c__ * c1_3.q * qr5 - *z__));
        /*<         DU( 2)= ALP4/R3*(S2D-F3*Y*S/R2)+ALP5*F3/R5*(C*(T-Y+Y*Q*QR5)-Y*Z) >*/
        du[1] = c0_1.alp4 / c1_3.r3 * (c0_1.s2d - f3 * *y * c1_3.s / c1_3.r2) +
                c0_1.alp5 * f3 / c1_3.r5 * (c__ * (c1_3.t - *y + *y * c1_3.q * qr5) - *y * *z__);
        /*<         DU( 3)=-ALP4/R3*(F1-A3*SDSD)   -ALP5*F3/R5*(C*(S-D+D*Q*QR5)-D*Z) >*/
        du[2] =
            -c0_1.alp4 / c1_3.r3 * (f1 - c1_3.a3 * c0_1.sdsd) -
            c0_1.alp5 * f3 / c1_3.r5 * (c__ * (c1_3.s - *d__ + *d__ * c1_3.q * qr5) - *d__ * *z__);
        /*<         DU( 4)=-ALP4*F3*S/R5*A5 +ALP5*(C*QR*QR5*A7-F3*Z/R5*A5) >*/
        du[3] = -c0_1.alp4 * f3 * c1_3.s / c1_3.r5 * c1_3.a5 +
                c0_1.alp5 * (c__ * c1_3.qr * qr5 * a7 - f3 * *z__ / c1_3.r5 * c1_3.a5);
        /*<        >*/
        du[4] = f3 * *x / c1_3.r5 *
                (-c0_1.alp4 * (c0_1.s2d - f5 * *y * c1_3.s / c1_3.r2) -
                 c0_1.alp5 * f5 / c1_3.r2 * (c__ * (c1_3.t - *y + *y * c1_3.q * qr7) - *y * *z__));
        /*<        >*/
        du[5] = f3 * *x / c1_3.r5 *
                (c0_1.alp4 * (f1 - (f2 + c1_3.a5) * c0_1.sdsd) +
                 c0_1.alp5 * f5 / c1_3.r2 *
                     (c__ * (c1_3.s - *d__ + *d__ * c1_3.q * qr7) - *d__ * *z__));
        /*<         DU( 7)= DU(5) >*/
        du[6] = du[4];
        /*<        >*/
        du[7] =
            f3 / c1_3.r5 *
            (-c0_1.alp4 * (f2 * *y * c0_1.s2d + c1_3.s * b5) -
             c0_1.alp5 *
                 (c__ * (f2 * c0_1.sdsd + f10 * *y * (c1_3.t - *y) / c1_3.r2 - c1_3.q * qr5 * b7) +
                  *z__ * b5));
        /*<        >*/
        du[8] = f3 / c1_3.r5 *
                (c0_1.alp4 * *y * (f1 - c1_3.a5 * c0_1.sdsd) +
                 c0_1.alp5 * (c__ * (f3 + c1_3.a5) * c0_1.s2d - *y * dr5 * (c__ * d7 + *z__)));
        /*<        >*/
        du[9] = f3 * *x / c1_3.r5 *
                (-c0_1.alp4 * (c0_1.c2d + c1_3.s * dr5) +
                 c0_1.alp5 * (f5 * c__ / c1_3.r2 * (c1_3.s - *d__ + *d__ * c1_3.q * qr7) - f1 -
                              *z__ * dr5));
        /*<        >*/
        du[10] = f3 / c1_3.r5 *
                 (c0_1.alp4 * (*d__ * b5 * c0_1.s2d - *y * c5 * c0_1.c2d) +
                  c0_1.alp5 *
                      (c__ * ((f3 + c1_3.a5) * c0_1.s2d - *y * dr5 * d7) - *y * (f1 + *z__ * dr5)));
        /*<        >*/
        du[11] =
            f3 / c1_3.r5 *
            (-c0_1.alp4 * *d__ * (f1 - c1_3.a5 * c0_1.sdsd) -
             c0_1.alp5 *
                 (c__ * (c0_1.c2d + f10 * *d__ * (c1_3.s - *d__) / c1_3.r2 - c1_3.q * qr5 * c7) +
                  *z__ * (f1 + c5)));
        /*<         DO 444 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   444   U(I)=U(I)+POT3/PI2*DU(I) >*/
            /* L444: */
            u[i__] += *pot3 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /* ========================================= */
    /* =====  INFLATE SOURCE CONTRIBUTION  ===== */
    /* ========================================= */
    /*<       IF(POT4.NE.F0) THEN >*/
    if (*pot4 != f0) {
        /*<         DU( 1)= ALP4*F3*X*D/R5 >*/
        du[0] = c0_1.alp4 * f3 * *x * *d__ / c1_3.r5;
        /*<         DU( 2)= ALP4*F3*Y*D/R5 >*/
        du[1] = c0_1.alp4 * f3 * *y * *d__ / c1_3.r5;
        /*<         DU( 3)= ALP4*C3/R3 >*/
        du[2] = c0_1.alp4 * c1_3.c3 / c1_3.r3;
        /*<         DU( 4)= ALP4*F3*D/R5*A5 >*/
        du[3] = c0_1.alp4 * f3 * *d__ / c1_3.r5 * c1_3.a5;
        /*<         DU( 5)=-ALP4*F15*XY*D/R7 >*/
        du[4] = -c0_1.alp4 * f15 * c1_3.xy * *d__ / r7;
        /*<         DU( 6)=-ALP4*F3*X/R5*C5 >*/
        du[5] = -c0_1.alp4 * f3 * *x / c1_3.r5 * c5;
        /*<         DU( 7)= DU(5) >*/
        du[6] = du[4];
        /*<         DU( 8)= ALP4*F3*D/R5*B5 >*/
        du[7] = c0_1.alp4 * f3 * *d__ / c1_3.r5 * b5;
        /*<         DU( 9)=-ALP4*F3*Y/R5*C5 >*/
        du[8] = -c0_1.alp4 * f3 * *y / c1_3.r5 * c5;
        /*<         DU(10)= DU(6) >*/
        du[9] = du[5];
        /*<         DU(11)= DU(9) >*/
        du[10] = du[8];
        /*<         DU(12)= ALP4*F3*D/R5*(F2+C5) >*/
        du[11] = c0_1.alp4 * f3 * *d__ / c1_3.r5 * (f2 + c5);
        /*<         DO 555 I=1,12 >*/
        for (i__ = 1; i__ <= 12; ++i__) {
            /*<   555   U(I)=U(I)+POT4/PI2*DU(I) >*/
            /* L555: */
            u[i__] += *pot4 / pi2 * du[i__ - 1];
        }
        /*<       ENDIF >*/
    }
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* uc0_ */

/*<       SUBROUTINE  DCCON0(ALPHA,DIP) >*/
/* Subroutine */ int dccon00_(doublereal *alpha, doublereal *dip) {
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

    /*<       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D >*/
    /*<       DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/ >*/
    /*<       DATA EPS/1.D-6/ >*/
    /* ----- */
    /*<       ALP1=(F1-ALPHA)/F2 >*/
    c0_1.alp1 = (f1 - *alpha) / f2;
    /*<       ALP2= ALPHA/F2 >*/
    c0_1.alp2 = *alpha / f2;
    /*<       ALP3=(F1-ALPHA)/ALPHA >*/
    c0_1.alp3 = (f1 - *alpha) / *alpha;
    /*<       ALP4= F1-ALPHA >*/
    c0_1.alp4 = f1 - *alpha;
    /*<       ALP5= ALPHA >*/
    c0_1.alp5 = *alpha;
    /* ----- */
    /*<       P18=PI2/360.D0 >*/
    p18 = pi2 / 360.;
    /*<       SD=DSIN(DIP*P18) >*/
    c0_1.sd = sin(*dip * p18);
    /*<       CD=DCOS(DIP*P18) >*/
    c0_1.cd = cos(*dip * p18);
    /*<       IF(DABS(CD).LT.EPS) THEN >*/
    if (fabs(c0_1.cd) < eps) {
        /*<         CD=F0 >*/
        c0_1.cd = f0;
        /*<         IF(SD.GT.F0) SD= F1 >*/
        if (c0_1.sd > f0) {
            c0_1.sd = f1;
        }
        /*<         IF(SD.LT.F0) SD=-F1 >*/
        if (c0_1.sd < f0) {
            c0_1.sd = -f1;
        }
        /*<       ENDIF >*/
    }
    /*<       SDSD=SD*SD >*/
    c0_1.sdsd = c0_1.sd * c0_1.sd;
    /*<       CDCD=CD*CD >*/
    c0_1.cdcd = c0_1.cd * c0_1.cd;
    /*<       SDCD=SD*CD >*/
    c0_1.sdcd = c0_1.sd * c0_1.cd;
    /*<       S2D=F2*SDCD >*/
    c0_1.s2d = f2 * c0_1.sdcd;
    /*<       C2D=CDCD-SDSD >*/
    c0_1.c2d = c0_1.cdcd - c0_1.sdsd;
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* dccon00_ */

/*<       SUBROUTINE  DCCON1(X,Y,D) >*/
/* Subroutine */ int dccon1_(doublereal *x, doublereal *y, doublereal *d__) {
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f3 = 3.;
    static doublereal f5 = 5.;
    static doublereal eps = 1e-6;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal r7;

    /*<       IMPLICIT REAL*8 (A-H,O-Z) >*/

    /* ********************************************************************** */
    /* *****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    ***** */
    /* ********************************************************************** */

    /* ***** INPUT */
    /* *****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM */
    /* ### CAUTION ### IF X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO */

    /*<       COMMON /C0/DUMMY(5),SD,CD >*/
    /*<        >*/
    /*<       DATA  F0,F1,F3,F5,EPS/0.D0,1.D0,3.D0,5.D0,1.D-6/ >*/
    /* ----- */
    /*<       IF(DABS(X).LT.EPS) X=F0 >*/
    if (fabs(*x) < eps) {
        *x = f0;
    }
    /*<       IF(DABS(Y).LT.EPS) Y=F0 >*/
    if (fabs(*y) < eps) {
        *y = f0;
    }
    /*<       IF(DABS(D).LT.EPS) D=F0 >*/
    if (fabs(*d__) < eps) {
        *d__ = f0;
    }
    /*<       P=Y*CD+D*SD >*/
    c1_2.p = *y * c0_2.cd + *d__ * c0_2.sd;
    /*<       Q=Y*SD-D*CD >*/
    c1_2.q = *y * c0_2.sd - *d__ * c0_2.cd;
    /*<       S=P*SD+Q*CD >*/
    c1_2.s = c1_2.p * c0_2.sd + c1_2.q * c0_2.cd;
    /*<       T=P*CD-Q*SD >*/
    c1_2.t = c1_2.p * c0_2.cd - c1_2.q * c0_2.sd;
    /*<       XY=X*Y >*/
    c1_2.xy = *x * *y;
    /*<       X2=X*X >*/
    c1_2.x2 = *x * *x;
    /*<       Y2=Y*Y >*/
    c1_2.y2 = *y * *y;
    /*<       D2=D*D >*/
    c1_2.d2 = *d__ * *d__;
    /*<       R2=X2+Y2+D2 >*/
    c1_2.r2 = c1_2.x2 + c1_2.y2 + c1_2.d2;
    /*<       R =DSQRT(R2) >*/
    c1_2.r__ = sqrt(c1_2.r2);
    /*<       IF(R.EQ.F0) RETURN >*/
    if (c1_2.r__ == f0) {
        return 0;
    }
    /*<       R3=R *R2 >*/
    c1_2.r3 = c1_2.r__ * c1_2.r2;
    /*<       R5=R3*R2 >*/
    c1_2.r5 = c1_2.r3 * c1_2.r2;
    /*<       R7=R5*R2 >*/
    r7 = c1_2.r5 * c1_2.r2;
    /* ----- */
    /*<       A3=F1-F3*X2/R2 >*/
    c1_2.a3 = f1 - f3 * c1_2.x2 / c1_2.r2;
    /*<       A5=F1-F5*X2/R2 >*/
    c1_2.a5 = f1 - f5 * c1_2.x2 / c1_2.r2;
    /*<       B3=F1-F3*Y2/R2 >*/
    c1_2.b3 = f1 - f3 * c1_2.y2 / c1_2.r2;
    /*<       C3=F1-F3*D2/R2 >*/
    c1_2.c3 = f1 - f3 * c1_2.d2 / c1_2.r2;
    /* ----- */
    /*<       QR=F3*Q/R5 >*/
    c1_2.qr = f3 * c1_2.q / c1_2.r5;
    /*<       QRX=F5*QR*X/R2 >*/
    c1_2.qrx = f5 * c1_2.qr * *x / c1_2.r2;
    /* ----- */
    /*<       UY=SD-F5*Y*Q/R2 >*/
    c1_2.uy = c0_2.sd - f5 * *y * c1_2.q / c1_2.r2;
    /*<       UZ=CD+F5*D*Q/R2 >*/
    c1_2.uz = c0_2.cd + f5 * *d__ * c1_2.q / c1_2.r2;
    /*<       VY=S -F5*Y*P*Q/R2 >*/
    c1_2.vy = c1_2.s - f5 * *y * c1_2.p * c1_2.q / c1_2.r2;
    /*<       VZ=T +F5*D*P*Q/R2 >*/
    c1_2.vz = c1_2.t + f5 * *d__ * c1_2.p * c1_2.q / c1_2.r2;
    /*<       WY=UY+SD >*/
    c1_2.wy = c1_2.uy + c0_2.sd;
    /*<       WZ=UZ+CD >*/
    c1_2.wz = c1_2.uz + c0_2.cd;
    /*<       RETURN >*/
    return 0;
    /*<       END >*/
} /* dccon1_ */
