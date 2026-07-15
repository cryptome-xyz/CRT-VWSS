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
        long L;                     // lifting bit length (normally L = t + 2 * lambda = T - lambda)
        std::vector<long> weights;  // sorted weight vector after scaling
    };

    struct Parties {
        NTL::ZZ P;  // the product of all parties' moduli
        long bit_P = 0;    // the bit length of P
        std::vector<Party> users; // all parties in the system
        std::vector<Party> A;  // Set A of parties with weights >= T
        std::vector<Party> B; // Set B of parties not in A 
    };

    struct RangeProof {
        std::vector<std::string> commitments; // commitments c_1,c_2,c_3 to the square decomposition of 4*x0*(B-x0)+1
        hash::Hash Delta;
        std::vector<std::pair<NTL::ZZ, NTL::ZZ>> responses; // (z_i, t_i) for i = 1..3
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

    // on input a set of parameters, set up the dealer and parties
    void setup_dealer_and_parties(const Params& params);

    // getters
    const Params& get_params() const;
    const Dealer& get_dealer() const;
    const Party& get_party(long id) const;
    const Parties& get_parties() const;
    // const Broadcast& get_broadcast() const;

    // Distribute shares and compute proof
    VWSS::Msg share() const;
   
     
    bool verify_party_in_A(const MSG_A& msg, const Broadcast& broadcast) const;

    bool verify_party_in_B1(const MSG_B1& msg, const NTL::ZZ& lambda_j, const NTL::ZZ& modulus, const Broadcast& broadcast) const;

    bool verify_party_in_B2(const MSG_B2& msg, const NTL::ZZ& lambda_j, const NTL::ZZ& modulus, const Broadcast& broadcast) const;

    // size calculation functions
    const std::vector<long> broadcast_size(const Broadcast& broadcast) const;
    const std::vector<long> msg_A_size(const std::vector<MSG_A>& msg_A) const ;
    const std::vector<long> msg_B1_size(const std::vector<MSG_B1>& msg_B1) const ;
    const std::vector<long> msg_B2_size(const std::vector<MSG_B2>& msg_B2) const ;
private:
    // for recording Merkle tree leaf of party in B with large weight
    struct PendingB1 {
        size_t msg_index;
        size_t leaf_index;
    };

    // for recording Merkle tree leaf and mask values of party in B with small weight
    struct PendingB2 {
        size_t msg_index;
        NTL::ZZ kj;
        NTL::ZZ bar_kj;
        IntCom::Com c_kj;
        IntCom::Com c_bar_kj;
        size_t leaf_index;
    };
    
    Params params_; // system parameters
    Dealer dealer_; // dealer
    Parties parties_;   // parties
    void setup_dealer();
    void setup_parties();
    NTL::ZZ generate_party_modulus(long bits, const std::set<NTL::ZZ>& used_primes) const;
    void build_AB_small();  // This is method 1 in the paper
    void build_AB_large();  // This is nmethod 2 in the paper

    // verify the proof of broadcast message
    bool verify_broadcast(const Broadcast& broadcast) const;

   

};

