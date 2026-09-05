#pragma once

extern "C" {
#include <pari/pari.h>
}

#include <NTL/ZZ.h>
#include <string>
#include <sstream>
#include <vector>

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
    std::string g_str() const;
    std::string h_str() const;
    std::string pow(const std::string& a, const NTL::ZZ& x) const;
    std::string mul(const std::string& a, const std::string& b) const;
    std::string multi_pow_mul(const std::vector<std::string>& bases, const std::vector<NTL::ZZ>& exps) const;
    std::string g_pow(const NTL::ZZ& x) const;
    std::string h_pow(const NTL::ZZ& y) const;
    // Precompute the fixed-base g/h tables before a timed operation.
    void warmup(long max_bits) const;

private:
    static constexpr int WINDOW_BITS = 4;
    static constexpr int WINDOW_DIGITS = (1 << WINDOW_BITS) - 1; 
    static constexpr int MULTI_WINDOW_BITS = 4;
    static constexpr int MULTI_WINDOW_DIGITS = (1 << MULTI_WINDOW_BITS) - 1;

    long U_; // the bit length upper bound of the group order
    GEN delta_;
    GEN g_;
    GEN h_;
    mutable std::vector<std::vector<GEN>> g_window_;
    mutable std::vector<std::vector<GEN>> h_window_;

    static GEN zzToGEN(const NTL::ZZ& x);
    static std::string genToString(GEN x);
    static GEN makeQfb(const std::string& a, const std::string& b, const std::string& c);

    static std::string serializeForm(GEN f);
    static GEN deserializeForm(const std::string& s);

    static void ensure_window_table(std::vector<std::vector<GEN>>& table, GEN base, long num_windows_needed);
    static GEN window_pow(std::vector<std::vector<GEN>>& table, GEN base, const NTL::ZZ& e);
};
