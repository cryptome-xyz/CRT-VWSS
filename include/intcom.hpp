// The integer commitment based on class groups of imaginary quadratic orders.
// It features a transparent setup (no trapdoor)
#pragma once

extern "C" {
#include <pari/pari.h>
}

#include <NTL/ZZ.h>
#include <string>
#include <sstream>

class IntCom {
public:
    struct Com {
        std::string c_x;
        NTL::ZZ r;   
    };

    IntCom();
    ~IntCom();

    // disable copy and move
    IntCom(const IntCom&) = delete;
    IntCom& operator=(const IntCom&) = delete;
    IntCom(IntCom&&) = delete;
    IntCom& operator=(IntCom&&) = delete;

    // commit to x with randomness r from [0,r_range)
    Com commit(const NTL::ZZ& x, const long r_range) const;
    bool open(const std::string& c_x, const NTL::ZZ& x, const NTL::ZZ& r) const;
    long U() const { return U_; }
    long commitment_bitlength(const std::string& c_x) const;
    std::string commit_elem(const NTL::ZZ& x, const NTL::ZZ& y) const;
    std::string pow(const std::string& a, const NTL::ZZ& x) const;

private:
    long U_; // the bit length upper bound of the group order
    GEN delta_;
    GEN g_;
    GEN h_;

    static GEN zzToGEN(const NTL::ZZ& x);
    static std::string genToString(GEN x);
    static GEN makeQfb(const std::string& a, const std::string& b, const std::string& c);

    static std::string serializeForm(GEN f);
    static GEN deserializeForm(const std::string& s);
};
