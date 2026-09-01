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
    // open the commitment by revealing x and r, return true if the commitment can be recomputed from the opening
    bool open(const std::string& c_x, const NTL::ZZ& x, const NTL::ZZ& r) const;
    // get upper bound on the group order (in bits)
    long U() const { return U_; }
    // get the bit length of the commitment   
    long commitment_bitlength(const std::string& c_x) const;  
    // compute commitment to x with randomness y
    std::string commit_elem(const NTL::ZZ& x, const NTL::ZZ& y) const;
    // serialized string form of the base generators g and h, for folding g^x/h^y terms
    // into a multi_pow_mul call alongside other bases
    std::string g_str() const;
    std::string h_str() const;
    // compute a^x where a is a group element serialized as a string
    std::string pow(const std::string& a, const NTL::ZZ& x) const;
    // multiply two group elements, each serialized as a string
    std::string mul(const std::string& a, const std::string& b) const;
    // compute prod_i bases[i]^exps[i] via windowed simultaneous (Straus's) exponentiation: each
    // base gets its own small precomputed table of multiples (base_i, 2*base_i, ..., (2^W-1)*
    // base_i in additive/multiplicative notation), then one shared pass advances the shared
    // accumulator W bits at a time, multiplying in each base's table entry for its W-bit digit.
    // Squaring count is unchanged from the plain (W=1) version -- advancing the accumulator by
    // the same total number of bits needs the same number of squarings regardless of window
    // size -- but multiplication count drops roughly W-fold, since each base contributes at
    // most one multiplication per W-bit window instead of up to one per single bit.
    std::string multi_pow_mul(const std::vector<std::string>& bases, const std::vector<NTL::ZZ>& exps) const;
    // compute g^x / h^y using a lazily-extended precomputed windowed table: for window size
    // W, table[k][j-1] = base^(j * 2^(W*k)) for j = 1..2^W-1. Since g and h are fixed for the
    // lifetime of this object, the squaring chain needed to reach any exponent's bit length,
    // plus each window's digit values, are computed at most once and reused by every
    // subsequent call: exponentiation becomes ~bits/W table-lookup multiplications instead of
    // ~bits/2 (plain binary table) or ~bits squarings + ~bits/2 multiplications (no table).
    std::string g_pow(const NTL::ZZ& x) const;
    std::string h_pow(const NTL::ZZ& y) const;
    // extend both the g and h windowed tables to cover exponents up to max_bits bits, up front.
    // Optional: g_pow/h_pow/commit/commit_elem already extend the tables lazily on first use;
    // call this beforehand only if you want that one-time cost paid outside of a timed region,
    // rather than absorbed by whichever call happens to need a larger exponent first.
    void warmup(long max_bits) const;
private:
    static constexpr int WINDOW_BITS = 4;
    static constexpr int WINDOW_DIGITS = (1 << WINDOW_BITS) - 1; // 15
    // Window size for multi_pow_mul's per-call tables (distinct from WINDOW_BITS above: that
    // one amortizes its setup cost across every call over an IntCom's whole lifetime, since g
    // and h never change; this one pays its setup cost fresh on every single call, since the
    // bases differ each time, so the optimal size is smaller).
    static constexpr int MULTI_WINDOW_BITS = 4;
    static constexpr int MULTI_WINDOW_DIGITS = (1 << MULTI_WINDOW_BITS) - 1;

    long U_; // the bit length upper bound of the group order
    GEN delta_;
    GEN g_;
    GEN h_;
    // g_window_[k][j-1] = g^(j * 2^(WINDOW_BITS*k)), gclone'd, lazily extended one window at a time
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
