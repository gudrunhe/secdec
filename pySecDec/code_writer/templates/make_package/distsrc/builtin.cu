#define SECDEC_RESULT_IS_COMPLEX 1
#include "common_cuda.h"

// Replacing BlockReduce with WarpReduce (while decreasing the
// block size) doesn't seem to work on CUDA11 + A100: the sum
// gets wrong randomly. Inserting __threadfence() helps, but
// __syncthreads() doesn't -- and I don't know why. Lets use
// BlockReduce both here and in the integrands, so that there
// would be only one thing to debug, not two.

#define sum_kernel(name, value_t) \
    extern "C" __global__ void \
    name(value_t *dst, value_t *src, uint64_t n) \
    { \
        uint64_t bid = blockIdx.x; \
        uint64_t tid = threadIdx.x; \
        uint64_t idx = (bid*128 + tid)*8; \
        value_t val1 = (idx+0 < n) ? src[idx+0] : value_t(0); \
        value_t val2 = (idx+1 < n) ? src[idx+1] : value_t(0); \
        value_t val3 = (idx+2 < n) ? src[idx+2] : value_t(0); \
        value_t val4 = (idx+3 < n) ? src[idx+3] : value_t(0); \
        value_t val5 = (idx+4 < n) ? src[idx+4] : value_t(0); \
        value_t val6 = (idx+5 < n) ? src[idx+5] : value_t(0); \
        value_t val7 = (idx+6 < n) ? src[idx+6] : value_t(0); \
        value_t val8 = (idx+7 < n) ? src[idx+7] : value_t(0); \
        value_t val = ((val1 + val2) + (val3 + val4)) + ((val5 + val6) + (val7 + val8)); \
        typedef cub::BlockReduce<value_t, 128, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY> Reduce; \
        __shared__ typename Reduce::TempStorage shared; \
        value_t sum = Reduce(shared).Sum(val); \
        if (tid == 0) dst[bid] = sum; \
    }

sum_kernel(sum_d_b128_x1024, real_t)
sum_kernel(sum_c_b128_x1024, complex_t)

#define SecDecInternalSignCheckErrorPositivePolynomial(id) {val = nan("U"); break;}
#define SecDecInternalSignCheckErrorContourDeformation(id) {val = nan("F"); break;}

extern "C" __global__ void
builtin__gauge( // sunset, nu=(1,2,3), realp=(q2, m1sq, m2sq, m3sq), sector=1, order=0
    result_t * __restrict__ result,
    const uint64_t lattice,
    const uint64_t index1,
    const uint64_t index2,
    const uint64_t * __restrict__ genvec,
    const real_t * __restrict__ shift,
    const real_t * __restrict__ realp,
    const complex_t * __restrict__ complexp,
    const real_t * __restrict__ deformp
)
{
    // assert(blockDim.x == 128);
    const uint64_t bid = blockIdx.x;
    const uint64_t tid = threadIdx.x;
    const real_t q2 = realp[0];
    const real_t m1sq = realp[1];
    const real_t m2sq = realp[2];
    const real_t m3sq = realp[3];
    const real_t SecDecInternalLambda0 = deformp[0];
    const real_t SecDecInternalLambda1 = deformp[1];
    const real_t invlattice = 1.0/lattice;
    result_t val = 0.0;
    uint64_t index = index1 + (bid*128 + tid)*8;
    uint64_t li_x0 = mulmod(index, genvec[0], lattice);
    uint64_t li_x1 = mulmod(index, genvec[1], lattice);
    for (uint64_t i = 0; (i < 8) && (index < index2); i++, index++) {
        real_t x0 = warponce(li_x0*invlattice + shift[0], 1.0);
        li_x0 = warponce_i(li_x0 + genvec[0], lattice);
        real_t x1 = warponce(li_x1*invlattice + shift[1], 1.0);
        li_x1 = warponce_i(li_x1 + genvec[1], lattice);
        real_t w_x0 = korobov3x3_w(x0);
        real_t w_x1 = korobov3x3_w(x1);
        real_t w = w_x0*w_x1;
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
        auto tmp1_7 = m1sq*tmp1_6;
        auto tmp1_8 = 2*x1;
        auto tmp1_9 = tmp3_2*tmp1_8;
        auto tmp3_3 = tmp1_9 + tmp1_4;
        auto tmp3_4 = tmp1_2*x1;
        auto tmp1_10 = x0*m1sq;
        auto tmp1_11 = tmp1_8*tmp1_10;
        auto tmp1_12 = tmp1_1*x0;
        auto tmp3_5 = tmp1_12 + m1sq;
        auto tmp3_6 = tmp1_11 + tmp3_5;
        auto tmp1_13 = m3sq + tmp3_4;
        auto tmp1_14 = 1 + x0;
        auto tmp3_7 = m2sq*tmp1_14;
        auto tmp3_8 = x1*tmp3_5;
        auto tmp3_9 = tmp3_8 + tmp3_7;
        auto tmp3_10 = x1*m3sq;
        auto tmp1_15 = 3*x1;
        auto tmp1_16 = x0*SecDecInternalLambda0;
        auto tmp1_17 = -SecDecInternalLambda0 + tmp1_16;
        auto tmp1_18 = x1*SecDecInternalLambda1;
        auto tmp1_19 = -SecDecInternalLambda1 + tmp1_18;
        auto tmp1_20 = -1 + tmp1_8;
        auto tmp3_11 = SecDecInternalLambda1*tmp1_20;
        auto tmp3_12 = -1 + tmp1_6;
        auto tmp3_13 = SecDecInternalLambda0*tmp3_12;
        auto __PowCall1 = x0*x0;
        auto __PowCall5 = x1*x1;
        auto __DenominatorCall2 = SecDecInternalDenominator(x0);
        auto __DenominatorCall3 = SecDecInternalDenominator(x1);
        auto tmp2_11 = m3sq*__PowCall5;
        auto tmp3_14 = tmp2_11 + tmp3_10;
        auto tmp3_15 = __PowCall1*tmp3_14;
        auto tmp2_12 = tmp1_10*__PowCall5;
        auto tmp3_16 = tmp2_12 + tmp3_9 + tmp3_15;
        auto tmp3_17 = tmp1_13*__PowCall1;
        auto tmp3_18 = tmp3_6 + tmp3_17;
        auto tmp2_13 = tmp1_2*__PowCall5;
        auto tmp3_19 = tmp2_13 + tmp3_4;
        auto tmp2_14 = tmp1_2*__PowCall1;
        auto tmp3_20 = tmp2_14 + tmp1_7;
        auto tmp2_15 = tmp3_2*__PowCall5;
        auto tmp3_21 = tmp2_15 + tmp3_1;
        auto __PowCall3 = __DenominatorCall2*__DenominatorCall2;
        auto __PowCall4 = __DenominatorCall3*__DenominatorCall3;
        auto __RealPartCall3 = SecDecInternalRealPart(tmp3_3);
        auto tmp3_22 = __DenominatorCall3*__PowCall3;
        auto tmp3_23 = __PowCall4*__DenominatorCall2*__PowCall3;
        auto tmp3_24 = SecDecInternalLambda0*__PowCall1;
        auto tmp3_25 = -tmp1_16 + tmp3_24;
        auto tmp3_26 = SecDecInternalI(__RealPartCall3);
        auto tmp3_27 = tmp3_26*tmp3_25;
        auto tmp3_28 = SecDecInternalLambda1*__PowCall5;
        auto tmp3_29 = -tmp1_18 + tmp3_28;
        auto tmp3_30 = tmp3_26*tmp3_29;
        auto _logCall1 = SecDecInternalLog(tmp3_22);
        auto _logCall2 = SecDecInternalLog(tmp3_23);
        auto __RealPartCall1 = SecDecInternalRealPart(tmp3_18);
        auto __RealPartCall2 = SecDecInternalRealPart(tmp3_19);
        auto __RealPartCall4 = SecDecInternalRealPart(tmp3_20);
        auto __RealPartCall5 = SecDecInternalRealPart(tmp3_21);
        auto tmp3_31 = SecDecInternalI(__RealPartCall5);
        auto tmp3_32 = tmp1_17*tmp3_31;
        auto tmp3_33 = 1 + tmp3_32;
        auto tmp3_34 = SecDecInternalI(__RealPartCall1);
        auto tmp3_35 = tmp1_19*tmp3_34;
        auto tmp3_36 = 1 + tmp3_35;
        auto tmp3_37 = __PowCall1*SecDecInternalLambda0;
        auto tmp3_38 = tmp3_37-tmp1_16;
        auto tmp3_39 = tmp3_38*tmp3_31;
        auto tmp3_40 = x0 + tmp3_39;
        auto tmp3_41 = __RealPartCall2*tmp3_38;
        auto tmp3_42 = __RealPartCall5*tmp3_13;
        auto tmp3_43 = tmp3_42 + tmp3_41;
        auto tmp3_44 = SecDecInternalI(tmp3_43);
        auto tmp3_45 = 1 + tmp3_44;
        auto tmp3_46 = __PowCall5*SecDecInternalLambda1;
        auto tmp3_47 = tmp3_46-tmp1_18;
        auto tmp3_48 = tmp3_47*tmp3_34;
        auto tmp3_49 = x1 + tmp3_48;
        auto tmp3_50 = __RealPartCall4*tmp3_47;
        auto tmp3_51 = __RealPartCall1*tmp3_11;
        auto tmp3_52 = tmp3_51 + tmp3_50;
        auto tmp3_53 = SecDecInternalI(tmp3_52);
        auto tmp3_54 = 1 + tmp3_53;
        auto tmp3_55 = tmp3_49+1;
        auto tmp3_56 = tmp3_40*tmp3_55;
        auto tmp3_57 = 1 + tmp3_56;
        auto tmp3_58 = -tmp3_27*tmp3_30;
        auto tmp3_59 = tmp3_45*tmp3_54;
        auto tmp3_60 = tmp3_58 + tmp3_59;
        auto _logCall5 = SecDecInternalLog(tmp3_36);
        auto __PowCall2 = tmp3_33*tmp3_33;
        auto __PowCall6 = tmp3_40*tmp3_40;
        auto __PowCall7 = tmp3_49*tmp3_49;
        auto tmp3_61 = __PowCall6*m3sq;
        auto tmp3_62 = tmp3_40*tmp1_1;
        auto tmp3_63 = tmp3_62 + m1sq + tmp3_61;
        auto tmp3_64 = tmp3_49*tmp3_63;
        auto tmp3_65 = __PowCall7*tmp3_61;
        auto tmp3_66 = __PowCall7*m1sq;
        auto tmp3_67 = m2sq + tmp3_66;
        auto tmp3_68 = tmp3_40*tmp3_67;
        auto tmp3_69 = tmp3_64 + tmp3_68 + m2sq + tmp3_65;
        auto _logCall4 = SecDecInternalLog(tmp3_57);
        auto _logCall3 = SecDecInternalLog(tmp3_69);
        auto __PowCall8 = tmp3_69*tmp3_69;
        auto __DenominatorCall1 = SecDecInternalDenominator(__PowCall8);
        auto tmp3_70 = -tmp3_16 + tmp3_69;
        auto tmp3_71 = -_logCall2-_logCall3;
        auto tmp3_72 = tmp1_8*tmp3_71;
        auto tmp3_73 = _logCall1 + _logCall4;
        auto tmp3_74 = tmp1_15*tmp3_73;
        auto tmp3_75 = x1*_logCall5;
        auto tmp3_76 = tmp3_75 + tmp3_74 + tmp3_72;
        auto tmp3_77 = __DenominatorCall1*tmp3_76*tmp3_60*tmp3_36*__PowCall1*__PowCall2;
        auto _SignCheckExpression = SecDecInternalImagPart(tmp3_70);
        if (unlikely(_SignCheckExpression>0)) SecDecInternalSignCheckErrorContourDeformation(1);
        auto tmp3_78 = SecDecInternalRealPart(tmp3_57);
        if (unlikely(tmp3_78<0)) SecDecInternalSignCheckErrorPositivePolynomial(1);
        val += w*(tmp3_77);
    }
    // Sum up 128*8=1024 values across 4 warps.
    typedef cub::BlockReduce<result_t, 128, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY> Reduce;
    __shared__ typename Reduce::TempStorage shared;
    result_t sum = Reduce(shared).Sum(val);
    if (tid == 0) result[bid] = sum;
}
