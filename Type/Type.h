//
// Created by convi on 2020/3/22.
//

#ifndef NUMERICALANALYSIS_TYPE_H
#define NUMERICALANALYSIS_TYPE_H

#include <cstdint>
#include <vector>

namespace convivae {

#define U1_MAX UINT8_MAX
#define U2_MAX UINT16_MAX
#define U4_MAX UINT32_MAX
#define U8_MAX UINT64_MAX
#define U4_INF 0x3f3f3f3f

    using u1 = std::uint8_t;
    using u2 = std::uint16_t;
    using u4 = std::uint32_t;
    using u8 = std::uint64_t;

    using i1 = std::int8_t;
    using i2 = std::int16_t;
    using i4 = std::int32_t;
    using i8 = std::int64_t;

    using f4 = float;
    using f8 = double;

    using double_t = f8;
}

#endif //NUMERICALANALYSIS_TYPE_H
