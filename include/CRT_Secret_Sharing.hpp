#pragma once

#include "class_group.hpp"
#include <map>
#include <set>
#include <string>
#include <vector>
#include <stdexcept>

using namespace NTL;

class CRT_Secret_Sharing {
public:
    struct PartyInfo {
        long id;
        long weight;
        NTL::ZZ modulus;   // p_i
    };

    struct Share {
        long party_id;
        NTL::ZZ value;     // S mod p_i
    };

    struct Params {
        long lambda = 128;         // security parameter
        long privacy_t = 0;        // unauthorized if total weight <= t
        long recon_T = 0;          // authorized if total weight >= T
        long scale_c = 1;          // scale factor for weighted construction
        NTL::ZZ p0;                // secret field modulus
        NTL::ZZ L;                 // lifting range
        std::vector<PartyInfo> parties;
    };

public:
    CRT_Secret_Sharing();
    explicit CRT_Secret_Sharing(const Params& params);
    virtual ~CRT_Secret_Sharing() = default;

    // ---------- Setup ----------
    virtual void set_params(const Params& params);
    virtual const Params& get_params() const;

    // Utility setup for weighted ramp CRT sharing.
    // This follows the paper's idea: choose p_i with bitlength about c * w_i,
    // and choose L around 2^(c*t + lambda). Distinct primes are used for coprimality.
    virtual void setup_weighted(
        long lambda,
        const std::vector<long>& weights,
        long privacy_t,
        long recon_T,
        long scale_c = 1,
        long min_modulus_bits = 16
    );

    // ---------- Core operations ----------
    virtual std::vector<Share> share_secret(const NTL::ZZ& secret) const;
    virtual NTL::ZZ reconstruct_secret(const std::vector<Share>& subset_shares) const;
    virtual NTL::ZZ reconstruct_lifted_value(const std::vector<Share>& subset_shares) const;

    // ---------- Access structure helpers ----------
    virtual bool is_authorized(const std::vector<long>& party_ids) const;
    virtual bool is_unauthorized(const std::vector<long>& party_ids) const;
    virtual long total_weight(const std::vector<long>& party_ids) const;

    // ---------- Validation ----------
    virtual bool validate_basic_coprimality() const;
    virtual bool validate_correctness_bound() const;
    virtual bool validate_secret(const NTL::ZZ& secret) const;

    // ---------- Convenience ----------
    virtual PartyInfo get_party(long party_id) const;
    virtual std::vector<long> all_party_ids() const;
    virtual std::map<long, Share> index_shares_by_party(const std::vector<Share>& shares) const;

protected:
    Params params_;
    std::map<long, PartyInfo> party_map_;

protected:
    virtual void rebuild_party_map();

    // Arithmetic helpers
    virtual NTL::ZZ mod_pos(const NTL::ZZ& x, const NTL::ZZ& m) const;
    virtual NTL::ZZ random_below(const NTL::ZZ& bound) const;
    virtual NTL::ZZ product_of_moduli(const std::vector<long>& party_ids) const;
    virtual NTL::ZZ compute_lambda_i(long target_party_id,
                                     const std::vector<long>& subset_party_ids) const;
    virtual std::vector<long> extract_party_ids(const std::vector<Share>& shares) const;
    virtual void ensure_distinct_party_ids(const std::vector<long>& party_ids) const;
    virtual void ensure_subset_is_known(const std::vector<long>& party_ids) const;

    // Parameter generation helpers
    virtual NTL::ZZ generate_distinct_prime(long bits,
                                            const std::set<std::string>& forbidden_values) const;
    virtual NTL::ZZ power_of_two(long exp) const;

    // Internal consistency checks
    virtual void ensure_initialized() const;
};
