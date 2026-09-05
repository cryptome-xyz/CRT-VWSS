#pragma once

#include "intcom.hpp"
#include "hash.hpp"
#include "three_square.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include <openssl/rand.h>
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
        long lambda = 128;         
        long t = 0;                 
        long T = 0;                
        NTL::ZZ p0;                
        long L;                     
        std::vector<long> weights; 
    };

    struct Parties {
        NTL::ZZ P;  
        long bit_P = 0;    
        NTL::ZZ PA; 
        std::vector<Party> users; 
        std::vector<Party> A;  
        std::vector<Party> B;
        long total_weights_in_A;
    };

    struct RangeProof {
        std::vector<std::string> commitments; 
        hash::Hash Delta;
        std::vector<std::pair<NTL::ZZ, NTL::ZZ>> responses; 
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

    // Messages sent to the parties in A
    struct MSG_A {
        long id;
        std::size_t msg_id;
        NTL::ZZ vi;
        NTL::ZZ rvi;
    };


    // Messages sent to the parties in B with weight >= lambda
    struct MSG_B1 {
        long id;
        std::size_t msg_id;
        NTL::ZZ vj;
        hash::AuthPath pathj;
        unsigned char rho[32];
    };

    // Messages sent to the parties in B with weight < lambda
    struct MSG_B2 {
        long id;
        std::size_t msg_id;
        NTL::ZZ vj;
        NTL::ZZ Rkj;
        NTL::ZZ Rrkj;
        std::string ckj;
        hash::AuthPath pathj;
        unsigned char rho[32];
    };

    struct Msg {
       Broadcast broadcast;
       std::vector<MSG_A> msg_A;
       std::vector<MSG_B1> msg_B1;
       std::vector<MSG_B2> msg_B2;
    };

    // Constructors and destructors
    VWSS();
    explicit VWSS(const Params& params, bool method1);
    ~VWSS() = default;

    // on input a set of parameters, set up the dealer and parties
    void setup_dealer_and_parties(const Params& params, bool method1);

    // Call after setup to precompute the fixed-base tables outside benchmark timings.
    void warmup() const;

    // getters
    const Params& get_params() const;
    const Dealer& get_dealer() const;
    const Party& get_party(long id) const;
    const Parties& get_parties() const;

    // Distribute shares and compute proof
    VWSS::Msg share() const;
   

    // verify the proof of broadcast message
    bool verify_broadcast(const Broadcast& broadcast) const;
     
    bool verify_party_in_A(const MSG_A& msg, const Broadcast& broadcast) const;

    bool verify_party_in_B1(const MSG_B1& msg, const NTL::ZZ& lambda_j, const NTL::ZZ& modulus, const Broadcast& broadcast) const;

    bool verify_party_in_B2(const MSG_B2& msg, const NTL::ZZ& lambda_j, const NTL::ZZ& modulus, const Broadcast& broadcast) const;

    // size calculation functions
    const std::vector<long> broadcast_size(const Broadcast& broadcast) const;
    const std::vector<long> msg_A_size(const std::vector<MSG_A>& msg_A) const ;
    const std::vector<long> msg_B1_size(const std::vector<MSG_B1>& msg_B1) const ;
    const std::vector<long> msg_B2_size(const std::vector<MSG_B2>& msg_B2) const ;
private:
    struct PendingB1 {
        size_t msg_index;
        size_t leaf_index;
    };

    struct PendingB2 {
        size_t msg_index;
        NTL::ZZ kj;
        NTL::ZZ bar_kj;
        IntCom::Com c_kj;
        IntCom::Com c_bar_kj;
        size_t leaf_index;
    };
    
    Params params_;
    Dealer dealer_; 
    Parties parties_;   
    IntCom intcom_; 
    void setup_dealer();
    void setup_parties(bool method1);
    NTL::ZZ generate_party_modulus(long bits, const std::set<NTL::ZZ>& used_primes) const;
    void build_AB_small();  // This is method 1 in the paper
    void build_AB_large();  // This is nmethod 2 in the paper


   

};

