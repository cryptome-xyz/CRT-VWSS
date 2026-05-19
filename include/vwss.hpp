// TODO: add verification of VWSS
// TODO: add range proof
// TODO: revise the code for hash.cpp
#pragma once

#include "intcom.hpp"
#include "hash.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <vector>
#include <set>
#include <algorithm>
class VWSS { 
public: 

struct Party {
        long id;
        long weight;    // weight after scaling
        NTL::ZZ modulus;
        NTL::ZZ share;
        NTL::ZZ global_CRT_coef;
    };

    struct Dealer {
        NTL::ZZ_p secret;
        NTL::ZZ lifted_secret;
    };

    struct Params {
        long lambda = 128;         // security parameter
        long t = 0;                 // privacy threhshold
        long T = 0;                //reconstruction threshold
        NTL::ZZ p0;                // secret field modulus
        long L;                     // lifting bit length (normally L = t + lambda = T - lambda)
        std::vector<long> weights;  // sorted weight vector after scaling
    };

    struct Parties {
        NTL::ZZ P;  // the product of all parties' moduli
        std::vector<Party> users; // all parties in the system
        std::vector<Party> A;  // Set A of parties with weightes > T
        std::vector<Party> B; // Set B of parties not in A 
    };

    struct RangeProof {
        NTL::ZZ z_1;
        NTL::ZZ z_2;
        NTL::ZZ z_3;
        NTL::ZZ t_1;
        NTL::ZZ t_2;
        NTL::ZZ t_3;
        std::string c1;
        std::string c2;
        std::string c3;
        NTL::ZZ tau;
    };

    struct Broadcast {
        // depending on whether A contains weight < lambda, 
        // if so, then ck, Rk, Rrk are defined
        std::string cs;
        std::vector<std::string> cv;
        std::string ck; 
        NTL::ZZ Rs;
        NTL::ZZ Rrs;
        NTL::ZZ Rk;
        NTL::ZZ Rrk;
        NTL::ZZ Rv;
        NTL::ZZ Rrv;
        hash::Hash rt;
        hash::Hash gamma;
        VWSS::RangeProof RP;
    };

    struct MSG_A {
        long id;
        NTL::ZZ vi;
        NTL::ZZ rvi;
    };

    struct MSG_B1 {
        long id;
        NTL::ZZ vj;
        NTL::ZZ Deltaj;
        hash::AuthPath pathj;
    };
    struct MSG_B2 {
        long id;
        NTL::ZZ vj;
        NTL::ZZ Deltaj;
        NTL::ZZ Rkj;
        NTL::ZZ Rrkj;
        std::string ckj;
        hash::AuthPath pathj;
    };

    struct Msg {
       Broadcast broadcast;
       std::vector<MSG_A> msg_A;
       std::vector<MSG_B1> msg_B1;
       std::vector<MSG_B2> msg_B2;
    };

    // Constructors and destructors
    VWSS();
    explicit VWSS(const Params& params);
    ~VWSS() = default;

    // Setup parameters
    void set_params(const Params& params);

    // getters
    const Params& get_params() const;
    const Dealer& get_dealer() const;
    const Party& get_party(long id) const;
    const Parties& get_parties() const;
    // const Broadcast& get_broadcast() const;

    void share() const;
    void verify_broadcast(const Broadcast& broadcast) const;

private:
    // system parameters
    Params params_;
    // dealer information
    Dealer dealer_;
    // parties information
    Parties parties_;
    void setup_dealer();
    void setup_parties();
    NTL::ZZ generate_party_modulus(long bits, const std::set<NTL::ZZ>& used_primes) const;
    void build_AB_small();
    void build_AB_large();
    const std::vector<long> broadcast_size(const Broadcast& broadcast) const ;
    const std::vector<long> msg_A_size(const std::vector<MSG_A>& msg_A) const ;
    const std::vector<long> msg_B1_size(const std::vector<MSG_B1>& msg_B1) const ;
    const std::vector<long> msg_B2_size(const std::vector<MSG_B2>& msg_B2) const ;

};

