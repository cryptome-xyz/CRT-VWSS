#pragma once
#include <NTL/ZZ.h>
#include <vector>

namespace three_square {
    bool decompose(NTL::ZZ n, std::vector<NTL::ZZ>& decompose);
    bool verify(const NTL::ZZ& n, const std::vector<NTL::ZZ>& decompose);
}