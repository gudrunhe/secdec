#define SECDEC_RESULT_IS_COMPLEX 1
#include "common_cpu.h"

#define SecDecInternalSignCheckErrorPositivePolynomial(id) {*presult = nan("U"); return 1; }
#define SecDecInternalSignCheckErrorContourDeformation(id) {*presult = nan("F"); return 2; }

extern "C" int
builtin__gauge( // sunset, nu=(1,2,3), realp=(q2, m1sq, m2sq, m3sq), sector=1, order=0
    result_t * restrict presult,
    const uint64_t lattice,
    const uint64_t index1,
    const uint64_t index2,
    const uint64_t * restrict genvec,
    const real_t * restrict shift,
    const real_t * restrict realp,
    const complex_t * restrict complexp,
    const real_t * restrict deformp
)
{
    const real_t q2 = realp[0];
    const real_t m1sq = realp[1];
    const real_t m2sq = realp[2];
    const real_t m3sq = realp[3];
    const real_t SecDecInternalLambda0 = deformp[0];
    const real_t SecDecInternalLambda1 = deformp[1];
    const real_t invlattice = 1.0/lattice;
    resultvec_t acc = RESULTVEC_ZERO;
    uint64_t index = index1;
    int_t li_x0 = mulmod(genvec[0], index, lattice);
    int_t li_x1 = mulmod(genvec[1], index, lattice);
    for (; index < index2; index += 4) {
        int_t li_x0_0 = li_x0; li_x0 = warponce_i(li_x0 + genvec[0], lattice);
        int_t li_x0_1 = li_x0; li_x0 = warponce_i(li_x0 + genvec[0], lattice);
        int_t li_x0_2 = li_x0; li_x0 = warponce_i(li_x0 + genvec[0], lattice);
        int_t li_x0_3 = li_x0; li_x0 = warponce_i(li_x0 + genvec[0], lattice);
        realvec_t x0 = {{ li_x0_0*invlattice, li_x0_1*invlattice, li_x0_2*invlattice, li_x0_3*invlattice }};
        x0 = warponce(x0 + shift[0], 1);
        int_t li_x1_0 = li_x1; li_x1 = warponce_i(li_x1 + genvec[1], lattice);
        int_t li_x1_1 = li_x1; li_x1 = warponce_i(li_x1 + genvec[1], lattice);
        int_t li_x1_2 = li_x1; li_x1 = warponce_i(li_x1 + genvec[1], lattice);
        int_t li_x1_3 = li_x1; li_x1 = warponce_i(li_x1 + genvec[1], lattice);
        realvec_t x1 = {{ li_x1_0*invlattice, li_x1_1*invlattice, li_x1_2*invlattice, li_x1_3*invlattice }};
        x1 = warponce(x1 + shift[1], 1);
        auto w_x0 = korobov3x3_w(x0);
        auto w_x1 = korobov3x3_w(x1);
        auto w = w_x0*w_x1;
        if (unlikely(index + 1 >= index2)) w.x[1] = 0;
        if (unlikely(index + 2 >= index2)) w.x[2] = 0;
        if (unlikely(index + 3 >= index2)) w.x[3] = 0;
        x0 = korobov3x3_f(x0);
        x1 = korobov3x3_f(x1);
        auto tmp1_1 = -q2 + m2sq + m3sq + m1sq;
        auto tmp1_2 = 2*m3sq;
        auto tmp1_3 = tmp1_2*x0;
        auto tmp1_4 = tmp1_1 + tmp1_3;
        auto tmp1_5 = x1*tmp1_4;
        auto tmp3_1 = m2sq + tmp1_5;
        auto tmp3_2 = tmp1_3 + m1sq;
        auto tmp1_6 = 2*x0;
        auto tmp1_7 = tmp1_6*m1sq;
        auto tmp1_8 = 2*x1;
        auto tmp3_3 = tmp3_2*tmp1_8;
        auto tmp3_4 = tmp3_3 + tmp1_4;
        auto tmp3_5 = tmp1_2*x1;
        auto tmp1_9 = tmp1_1*x0;
        auto tmp3_6 = tmp1_9 + m1sq;
        auto tmp1_10 = x1*tmp1_7;
        auto tmp3_7 = tmp1_10 + tmp3_6;
        auto tmp1_11 = m3sq + tmp3_5;
        auto tmp1_12 = 1 + x0;
        auto tmp3_8 = m2sq*tmp1_12;
        auto tmp3_9 = x1*tmp3_6;
        auto tmp3_10 = tmp3_9 + tmp3_8;
        auto tmp3_11 = x0*m1sq;
        auto tmp1_13 = x1*m3sq;
        auto tmp1_14 = x0*SecDecInternalLambda0;
        auto tmp1_15 = -SecDecInternalLambda0 + tmp1_14;
        auto tmp1_16 = x1*SecDecInternalLambda1;
        auto tmp1_17 = -SecDecInternalLambda1 + tmp1_16;
        auto tmp1_18 = -SecDecInternalLambda1+2*tmp1_16;
        auto tmp3_12 = -1 + tmp1_6;
        auto tmp3_13 = SecDecInternalLambda0*tmp3_12;
        auto __PowCall1 = x0*x0;
        auto __PowCall3 = x1*x1;
        auto tmp2_11 = m3sq*__PowCall3;
        auto tmp3_14 = tmp2_11 + tmp1_13;
        auto tmp3_15 = __PowCall1*tmp3_14;
        auto tmp2_12 = tmp3_11*__PowCall3;
        auto tmp3_16 = tmp3_10 + tmp2_12 + tmp3_15;
        auto tmp3_17 = tmp1_11*__PowCall1;
        auto tmp3_18 = tmp3_7 + tmp3_17;
        auto tmp2_13 = tmp1_2*__PowCall3;
        auto tmp3_19 = tmp2_13 + tmp3_5;
        auto tmp2_14 = tmp1_2*__PowCall1;
        auto tmp3_20 = tmp2_14 + tmp1_7;
        auto tmp2_15 = tmp3_2*__PowCall3;
        auto tmp3_21 = tmp2_15 + tmp3_1;
        auto __RealPartCall3 = SecDecInternalRealPart(tmp3_4);
        auto tmp3_22 = SecDecInternalI(__RealPartCall3);
        auto tmp3_23 = __PowCall1*SecDecInternalLambda0;
        auto tmp3_24 = tmp3_23-tmp1_14;
        auto tmp3_25 = tmp3_24*tmp3_22;
        auto tmp3_26 = __PowCall3*SecDecInternalLambda1;
        auto tmp3_27 = tmp3_26-tmp1_16;
        auto tmp3_28 = tmp3_27*tmp3_22;
        auto __RealPartCall1 = SecDecInternalRealPart(tmp3_18);
        auto __RealPartCall2 = SecDecInternalRealPart(tmp3_19);
        auto __RealPartCall4 = SecDecInternalRealPart(tmp3_20);
        auto __RealPartCall5 = SecDecInternalRealPart(tmp3_21);
        auto tmp3_29 = SecDecInternalI(__RealPartCall5);
        auto tmp3_30 = tmp1_15*tmp3_29;
        auto tmp3_31 = 1 + tmp3_30;
        auto tmp3_32 = SecDecInternalI(__RealPartCall1);
        auto tmp3_33 = tmp1_17*tmp3_32;
        auto tmp3_34 = 1 + tmp3_33;
        auto tmp3_35 = SecDecInternalLambda0*__PowCall1;
        auto tmp3_36 = tmp3_35-tmp1_14;
        auto tmp3_37 = tmp3_29*tmp3_36;
        auto tmp3_38 = x0 + tmp3_37;
        auto tmp3_39 = SecDecInternalI(tmp3_36*__RealPartCall2);
        auto tmp3_40 = tmp3_13*tmp3_29;
        auto tmp3_41 = tmp3_40+1 + tmp3_39;
        auto tmp3_42 = SecDecInternalLambda1*__PowCall3;
        auto tmp3_43 = tmp3_42-tmp1_16;
        auto tmp3_44 = tmp3_32*tmp3_43;
        auto tmp3_45 = x1 + tmp3_44;
        auto tmp3_46 = SecDecInternalI(tmp3_43*__RealPartCall4);
        auto tmp3_47 = tmp1_18*tmp3_32;
        auto tmp3_48 = tmp3_47+1 + tmp3_46;
        auto tmp3_49 = 1 + tmp3_45;
        auto tmp3_50 = tmp3_38*tmp3_49;
        auto tmp3_51 = 1 + tmp3_50;
        auto tmp3_52 = tmp3_48*tmp3_41;
        auto tmp3_53 = -tmp3_28*tmp3_25;
        auto tmp3_54 = tmp3_52 + tmp3_53;
        auto __PowCall2 = tmp3_31*tmp3_31;
        auto __PowCall4 = tmp3_38*tmp3_38;
        auto __PowCall5 = tmp3_45*tmp3_45;
        auto tmp3_55 = __PowCall4*m3sq;
        auto tmp3_56 = tmp3_38*tmp1_1;
        auto tmp3_57 = tmp3_56 + m1sq + tmp3_55;
        auto tmp3_58 = tmp3_45*tmp3_57;
        auto tmp3_59 = __PowCall5*tmp3_55;
        auto tmp3_60 = __PowCall5*m1sq;
        auto tmp3_61 = m2sq + tmp3_60;
        auto tmp3_62 = tmp3_38*tmp3_61;
        auto tmp3_63 = tmp3_58 + tmp3_62 + m2sq + tmp3_59;
        auto __PowCall6 = tmp3_63*tmp3_63;
        auto __DenominatorCall1 = SecDecInternalDenominator(__PowCall6);
        auto tmp3_64 = -tmp3_16 + tmp3_63;
        auto tmp3_65 = x1*tmp3_54*tmp3_34*__PowCall1*__PowCall2*__DenominatorCall1;
        auto _SignCheckExpression = SecDecInternalImagPart(tmp3_64);
        if (unlikely(_SignCheckExpression>0)) SecDecInternalSignCheckErrorContourDeformation(1);
        auto tmp3_66 = SecDecInternalRealPart(tmp3_51);
        if (unlikely(tmp3_66<0)) SecDecInternalSignCheckErrorPositivePolynomial(1);
        acc = acc + w*(tmp3_65);
    }
    *presult = componentsum(acc);
    return 0;
}
