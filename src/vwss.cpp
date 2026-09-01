#include "vwss.hpp"
#include <chrono>
#include <iostream>


// Constructors
VWSS::VWSS() = default;

VWSS::VWSS(const Params& params, bool method1){
    setup_dealer_and_parties(params, method1);
}

void VWSS::setup_dealer_and_parties(const Params& params, bool method1) {
    params_ = params;
    setup_dealer();
    setup_parties(method1);
}

void VWSS::setup_dealer() {
    NTL::ZZ_p::init(params_.p0);
    dealer_.secret = NTL::random_ZZ_p();
    NTL::ZZ u = NTL::RandomBits_ZZ(params_.L);
    // lifted secret = s + u * p0
    dealer_.lifted_secret = NTL::rep(dealer_.secret) + u * params_.p0;
}

void VWSS::setup_parties(bool method1) {
    // In CRT-based SS, primes must be distinct
    std::set<NTL::ZZ> used_primes;
    used_primes.insert(params_.p0);
    parties_.P = 1;
    parties_.users.clear();
    parties_.users.reserve(params_.weights.size());


    for(std::size_t i = 0; i < params_.weights.size(); ++i) {
        Party p;
        p.id = static_cast<long>(i);    
        p.weight = params_.weights[i];  
        // generate party modulus with p.weight bits
        p.modulus = generate_party_modulus(p.weight, used_primes);  
        parties_.P *= p.modulus;    
        parties_.bit_P += p.weight;
        p.share = dealer_.lifted_secret % p.modulus;    
        used_primes.insert(p.modulus);
        parties_.users.push_back(p);    
    }

    // compute global_CRT_coefficients for every party
    for(auto& p: parties_.users) {
        NTL::ZZ Qi = parties_.P / p.modulus;
        NTL::ZZ inv = NTL::InvMod(Qi % p.modulus, p.modulus);
        NTL::ZZ temp = NTL::MulMod(Qi, inv, parties_.P);
        // reduce temp to the range (-P/2, P/2) 
        if (temp > parties_.P / 2) {
            temp -= parties_.P;
        }
        p.global_CRT_coef = temp;
    }

    if (method1) {
        build_AB_small(); // method 1
    } else {
        build_AB_large(); // method 2
    }
}

const VWSS::Parties& VWSS::get_parties() const {
    return parties_;
}

const VWSS::Party& VWSS::get_party(long id) const {
    if (id < 0 || static_cast<std::size_t>(id) >= parties_.users.size()) {
        throw std::out_of_range("no party with the given id");
    }
    return parties_.users[id];
}

VWSS::Msg VWSS::share() const {
    using ShareClock = std::chrono::steady_clock;
    const auto t_start = ShareClock::now();
    // setup
    NTL::ZZ v(0); 
    NTL::ZZ rv(0); 
    NTL::ZZ bound_rv(0);
    VWSS::Msg msgs;

    const IntCom& intcom = intcom_;
    const long rnd_bit_length = intcom.U() + params_.lambda;
    const long E = parties_.bit_P - 1 + parties_.A.back().weight + ceil(log2(parties_.A.size())); // upper bound (in bit) on v
    const long Sv = parties_.bit_P - 1 + rnd_bit_length + ceil(log2(parties_.A.size())); // upper bound (in bit) on rv

    // commit to dealer's lifted secret s
    IntCom::Com c_s = intcom.commit(dealer_.lifted_secret, rnd_bit_length);
    msgs.broadcast.cs = c_s.c_x;

    // commit to parties in A's shares
    msgs.broadcast.cv.reserve(parties_.A.size());
    msgs.msg_A.reserve(parties_.A.size());


    for(const auto&p :parties_.A) {
        // compute v = \sum_{i \in A} \lambda_i * v_i
        v += p.global_CRT_coef * p.share;
        IntCom::Com c_vi = intcom.commit(p.share, rnd_bit_length);
        msgs.broadcast.cv.push_back(c_vi.c_x);
        // compute rv = sum_{ i\in A} \lambda_i * r_vi
        rv += p.global_CRT_coef * c_vi.r;

        VWSS::MSG_A msg;
        msg.vi = p.share;
        msg.rvi = c_vi.r;
        msg.id = p.id;
        msg.msg_id = msgs.msg_A.size();  
        msgs.msg_A.push_back(msg);
    }

    // commit to bar_s
    // first compute bar_s bit length upper bound (B')
    const long B_prime = std::max(params_.T, parties_.total_weights_in_A + params_.t);
    NTL::ZZ bar_s = NTL::RandomBits_ZZ(B_prime + 2 * params_.lambda);
    IntCom::Com c_bar_s = intcom.commit(bar_s, rnd_bit_length + 2 * params_.lambda);

    // commit to bar_v
    NTL::ZZ bar_v = NTL::RandomBits_ZZ(E + 2 * params_.lambda);
    IntCom::Com bar_c_v = intcom.commit(bar_v, Sv + 2 * params_.lambda);

    const NTL::ZZ bar_s_minus_bar_v = bar_s - bar_v;
    const NTL::ZZ dealer_secret_minus_v = dealer_.lifted_secret - v;
    std::vector<hash::Hash> hash_list;
    std::vector<PendingB1> pending_B1;
    std::vector<PendingB2> pending_B2;

    hash_list.reserve(parties_.B.size());

    for(const auto&p : parties_.B) {
        if(p.weight >= params_.lambda) {
            VWSS::MSG_B1 msg_j;
            msg_j.id = p.id;
            NTL::ZZ Deltaj = bar_s_minus_bar_v % p.modulus;
            msg_j.vj = p.share;
            msg_j.msg_id = msgs.msg_B1.size();
            msgs.msg_B1.push_back(msg_j);
            // hash the Delta_j
            hash_list.push_back(hash::hash(Deltaj)); 
            pending_B1.push_back({msgs.msg_B1.size() - 1, hash_list.size() - 1});
        } else {
            VWSS::MSG_B2 msg_j;
            msg_j.id = p.id;
            msg_j.vj = p.share;
            msg_j.msg_id = msgs.msg_B2.size();
            // generate kj 
            NTL::ZZ kj = (dealer_secret_minus_v - p.global_CRT_coef * p.share) / p.modulus;
            // commit to kj
            IntCom::Com c_kj = intcom.commit(kj, rnd_bit_length);
            msg_j.ckj = c_kj.c_x;
            // compute bar_kj 
            const long bar_kj_bit_length = 1 + std::max(parties_.bit_P - 1, std::max(params_.T, E) - p.weight) + 2 * params_.lambda;
            NTL::ZZ bar_kj = NTL::RandomBits_ZZ(bar_kj_bit_length);
            IntCom::Com c_bar_kj = intcom.commit(bar_kj, rnd_bit_length + 2 * params_.lambda);
            
            NTL::ZZ Deltaj = bar_s_minus_bar_v - p.modulus * bar_kj;
            msgs.msg_B2.push_back(msg_j);
            hash_list.push_back(hash::hashConcat(c_kj.c_x, c_bar_kj.c_x, Deltaj));

            pending_B2.push_back({msgs.msg_B2.size() - 1, kj, bar_kj, c_kj, c_bar_kj, hash_list.size() - 1});
        }
        
    }
    bool has_Merkle_tree = !hash_list.empty();
    std::unique_ptr<hash::MerkleTree> tree;

    if (has_Merkle_tree) {
        tree = std::make_unique<hash::MerkleTree>(hash_list);
        msgs.broadcast.rt = tree->root();
    } else {
        msgs.broadcast.rt = hash::Hash{};   // empty hash
    }

    const auto t_before_three_square = ShareClock::now();
    const NTL::ZZ B_bound = NTL::power2_ZZ(params_.T);
    const NTL::ZZ& x0 = dealer_.lifted_secret;
    std::vector<NTL::ZZ> rp_x(3); // x_1, x_2, x_3
    three_square::decompose(4 * x0 * (B_bound - x0) + 1, rp_x);
    const auto t_after_three_square = ShareClock::now();

    const long BCL_bit_length = params_.T + 2 * params_.lambda;  
    const long CLS_bit_length = intcom.U() + 3 * params_.lambda;
    std::vector<NTL::ZZ> rp_r(3), rp_m(3), rp_s(3);
    std::vector<std::string> rp_c(3), rp_d(3);
    for (int i = 0; i < 3; ++i) {
        IntCom::Com c_i = intcom.commit(rp_x[i], rnd_bit_length);
        rp_r[i] = c_i.r;
        rp_c[i] = c_i.c_x;

        rp_m[i] = NTL::RandomBits_ZZ(BCL_bit_length);
        IntCom::Com d_i = intcom.commit(rp_m[i], CLS_bit_length);
        rp_s[i] = d_i.r;
        rp_d[i] = d_i.c_x;
    }

    const long sigma_bit_length = params_.T + intcom.U() + 3 * params_.lambda + 2;
    const NTL::ZZ rp_sigma = NTL::RandomBits_ZZ(sigma_bit_length);
    std::string rp_d_elem = intcom.mul(
        intcom.h_pow(rp_sigma),
        intcom.multi_pow_mul(
            { c_s.c_x, rp_c[0], rp_c[1], rp_c[2] },
            { -4 * bar_s, -rp_m[0], -rp_m[1], -rp_m[2] }));

    std::ostringstream oss;
    std::ostringstream rp_oss;
    rp_oss << c_bar_s.c_x; 
    for (int i = 0; i < 3; i++) {
        rp_oss << rp_d[i];
    }
    rp_oss << rp_d_elem;
    const hash::Hash rp_Delta = hash::hash(rp_oss.str());

    msgs.broadcast.RP.commitments = rp_c;
    const auto t_after_rp_init = ShareClock::now();

    NTL::ZZ gamma_ZZ;
    NTL::ZZ bar_k(0);
    NTL::ZZ k(0);
    IntCom::Com c_k;
    IntCom::Com bar_c_k;

    NTL::ZZ delta;

    for(const auto& c : msgs.broadcast.cv) {
        oss << c;
    }

    for (const auto& c : rp_c) {
        oss << c;
    }
    oss.write(reinterpret_cast<const char*>(rp_Delta.data()), rp_Delta.size());
    

    if(parties_.A[0].weight < params_.lambda) {
        k = (dealer_.lifted_secret - v) / parties_.PA;
        c_k = intcom.commit(k, rnd_bit_length);
        msgs.broadcast.ck = c_k.c_x;

        // commit to bar_k
        const long bar_k_bit_length = std::max(params_.T, E) + 1 - parties_.total_weights_in_A + 2 * params_.lambda;
        bar_k = NTL::RandomBits_ZZ(bar_k_bit_length);
        bar_c_k = intcom.commit(bar_k, rnd_bit_length + 2 * params_.lambda);

        delta = bar_s_minus_bar_v - parties_.PA * bar_k;

        oss << c_s.c_x << c_bar_s.c_x << bar_c_v.c_x << c_k.c_x << bar_c_k.c_x;
    } else {
        msgs.broadcast.ck = ""; // not used

        delta = bar_s_minus_bar_v % parties_.PA;
        oss << c_s.c_x << c_bar_s.c_x << bar_c_v.c_x;
    }

    msgs.broadcast.gamma = hash::hashConcat(oss.str(), delta, msgs.broadcast.rt);
    // compute challenge (gamma); also serves as the range proof's challenge gamma in [0, 2^lambda)
    gamma_ZZ = hash::hashToZZ128(msgs.broadcast.gamma);

    // Responses
    msgs.broadcast.Rs = bar_s + dealer_.lifted_secret * gamma_ZZ;
    msgs.broadcast.Rrs = c_bar_s.r + c_s.r * gamma_ZZ;
    msgs.broadcast.Rv = bar_v + v * gamma_ZZ;
    msgs.broadcast.Rrv = bar_c_v.r + rv * gamma_ZZ;

    if(parties_.A[0].weight < params_.lambda){
        msgs.broadcast.Rk = bar_k + k * gamma_ZZ;
        msgs.broadcast.Rrk = bar_c_k.r + c_k.r * gamma_ZZ;
    } else {
        msgs.broadcast.Rk = 0; // not used
        msgs.broadcast.Rrk = 0; // not used
    }

    const auto t_before_rp_resp = ShareClock::now();
    // Range proof responses: z_i = m_i + gamma * x_i, t_i = s_i + gamma * r_i, i = 0..2
    msgs.broadcast.RP.responses.resize(3);
    for (int i = 0; i <= 2; ++i) {
        msgs.broadcast.RP.responses[i] = { rp_m[i] + gamma_ZZ * rp_x[i], rp_s[i] + gamma_ZZ * rp_r[i] };
    }
    // tau = sigma + gamma * (sum_{i=1..3} x_i * r_i - 4*(B - x_0) * r_0)
    msgs.broadcast.RP.tau = rp_sigma - gamma_ZZ * (4 * (B_bound - dealer_.lifted_secret) * c_s.r);
    for (int i = 0; i <= 2; ++i) {
        msgs.broadcast.RP.tau += gamma_ZZ * rp_x[i] * rp_r[i];
    }
    const auto t_after_rp_resp = ShareClock::now();

    if(has_Merkle_tree) {
        // compute Merkle tree authentication paths for parties in B1
        for (const auto&p : pending_B1) {
            msgs.msg_B1[p.msg_index].pathj = tree->authenticationPath(p.leaf_index);
        }

        // compute Merkle tree path and responses for parties in B with weight < lambda
        for (const auto&p : pending_B2) {
            msgs.msg_B2[p.msg_index].pathj = tree->authenticationPath(p.leaf_index);
            msgs.msg_B2[p.msg_index].Rkj =  p.bar_kj + p.kj* gamma_ZZ;
            msgs.msg_B2[p.msg_index].Rrkj = p.c_bar_kj.r + p.c_kj.r * gamma_ZZ;
        }
    }

    const auto t_end = ShareClock::now();
    const auto to_ms = [](auto d) {
        return std::chrono::duration<double, std::milli>(d).count();
    };
    const auto three_square_time = t_after_three_square - t_before_three_square;
    const auto range_proof_compute_time = (t_after_rp_init - t_after_three_square) + (t_after_rp_resp - t_before_rp_resp);
    const auto range_proof_total_time = three_square_time + range_proof_compute_time;
    const auto total_time = t_end - t_start;
    const auto proof_generation_time = total_time - range_proof_total_time;

    std::cout << "[share] proof generation (excl. range proof): " << to_ms(proof_generation_time) << " ms\n";
    std::cout << "[share] range proof - three-square decomposition: " << to_ms(three_square_time) << " ms\n";
    std::cout << "[share] range proof - computation (excl. three-square): " << to_ms(range_proof_compute_time) << " ms\n";
    std::cout << "[share] range proof total: " << to_ms(range_proof_total_time) << " ms\n";

    return msgs;
}

bool VWSS::verify_party_in_A(const MSG_A& msg, const Broadcast& broadcast) const {

    // check commitment c_vi ?= g^{vi} * h^{r_vi}
    return intcom_.open(broadcast.cv[msg.msg_id], msg.vi, msg.rvi) ;

}

bool VWSS::verify_party_in_B1(const MSG_B1& msg, const NTL::ZZ& lambda_j, const NTL::ZZ& modulus,const Broadcast& broadcast) const {

    NTL::ZZ gamma_ZZ = hash::hashToZZ128(broadcast.gamma);

    // compute Delta_j' = (R_s - R_v - gamma * lambda_j * v_j) mod p_j
    NTL::ZZ Delta_j_prime = ((broadcast.Rs - broadcast.Rv) - gamma_ZZ * lambda_j * msg.vj) % modulus;

    hash::Hash leaf = hash::hash(Delta_j_prime);
    hash::Hash recomputed_root = hash::MerkleTree::verify(leaf, msg.pathj);

    return (recomputed_root == broadcast.rt);
}

bool VWSS::verify_party_in_B2(const MSG_B2& msg, const NTL::ZZ& lambda_j, const NTL::ZZ& modulus, const Broadcast& broadcast) const {

    NTL::ZZ gamma_ZZ = hash::hashToZZ128(broadcast.gamma);

    const IntCom& intcom = intcom_;
    // compute bar_c_k_j_prime = g^Rkj h^Rrkj c_k_j^-gamma
    std::string bar_c_k_j_prime = intcom.mul(
        intcom.mul(intcom.g_pow(msg.Rkj), intcom.h_pow(msg.Rrkj)),
        intcom.pow(msg.ckj, -gamma_ZZ));

   

    // compute Delta_j_prime = R_s - R_v - gamma * lambda_j v_j - p_j R_kj
    NTL::ZZ Delta_j_prime = (broadcast.Rs - broadcast.Rv) - gamma_ZZ * lambda_j * msg.vj - modulus * msg.Rkj;

    hash::Hash leaf = hash::hashConcat(msg.ckj, bar_c_k_j_prime, Delta_j_prime);
    hash::Hash recomputed_root = hash::MerkleTree::verify(leaf, msg.pathj);

    return (recomputed_root == broadcast.rt);

}

bool VWSS::verify_broadcast(const Broadcast& broadcast) const {
    // verify the UPoM proof for c_s and c_v
    NTL::ZZ gamma_ZZ = hash::hashToZZ128(broadcast.gamma);
    const IntCom& intcom = intcom_;

    // bar_c_s' = g^Rs * h^Rrs * c_s^{-gamma}. g^Rs/h^Rrs
    std::string bar_c_s_prime = intcom.mul(
        intcom.mul(intcom.g_pow(broadcast.Rs), intcom.h_pow(broadcast.Rrs)),
        intcom.pow(broadcast.cs, -gamma_ZZ));

    // bar_c_v' = g^Rv * h^Rrv * (prod_{i in A} (c_vi)^{lambda_i})^{-gamma} 
    std::vector<std::string> cv_bases;
    std::vector<NTL::ZZ> cv_exps;
    cv_bases.reserve(parties_.A.size());
    cv_exps.reserve(parties_.A.size());
    for (size_t i = 0; i < parties_.A.size(); ++i) {
        cv_bases.push_back(broadcast.cv[i]);
        cv_exps.push_back(-gamma_ZZ * parties_.A[i].global_CRT_coef);
    }
    std::string bar_c_v_prime = intcom.mul(
        intcom.mul(intcom.g_pow(broadcast.Rv), intcom.h_pow(broadcast.Rrv)),
        intcom.multi_pow_mul(cv_bases, cv_exps));

    const NTL::ZZ Rs_minus_Rv = broadcast.Rs - broadcast.Rv;


    // f_i = g^{z_i} h^{t_i} * c_i^{-gamma}, i = 1..3
    std::vector<std::string> rp_f(3);
    for (std::size_t i = 0; i < 3; ++i) {
        rp_f[i] = intcom.mul(
            intcom.mul(intcom.g_pow(broadcast.RP.responses[i].first), intcom.h_pow(broadcast.RP.responses[i].second)),
            intcom.pow(broadcast.RP.commitments[i], -gamma_ZZ));
    }

    // f = h^tau * g^gamma * c_0^{4*(gamma*B - z_0)} * prod_{i=1..3} c_i^{-z_i}, B = 2^T, z_0 = Rs
    const NTL::ZZ B_bound = NTL::power2_ZZ(params_.T);
    std::string f = intcom.mul(
        intcom.mul(intcom.h_pow(broadcast.RP.tau), intcom.g_pow(gamma_ZZ)),
        intcom.multi_pow_mul(
            { broadcast.cs, broadcast.RP.commitments[0], broadcast.RP.commitments[1], broadcast.RP.commitments[2] },
            { 4 * (gamma_ZZ * B_bound - broadcast.Rs),
              -broadcast.RP.responses[0].first, -broadcast.RP.responses[1].first, -broadcast.RP.responses[2].first }));

    // Delta' = H(f_0, f_1, f_2, f_3, f)
    std::ostringstream rp_oss;
    rp_oss << bar_c_s_prime;
    for (std::size_t i = 0; i < 3; ++i) {
        rp_oss << rp_f[i];
    }
    rp_oss << f;
    const hash::Hash rp_Delta_prime = hash::hash(rp_oss.str());

    std::ostringstream oss;

    for (const auto& c : broadcast.cv) {
            oss << c;
    }
    for (const auto& c : broadcast.RP.commitments) {
            oss << c;
    }
    oss.write(reinterpret_cast<const char*>(rp_Delta_prime.data()), rp_Delta_prime.size());

    oss << broadcast.cs << bar_c_s_prime << bar_c_v_prime;

    hash::Hash gamma_prime;

    // check if A contains any party with weight < lambda
    if (parties_.A[0].weight >= params_.lambda) {
        // Delta_A' = R_s - R_v mod P_A
        NTL::ZZ Delta_A_prime = Rs_minus_Rv % parties_.PA;
        // gamma' ?= H(c_s, bar_c_s', bar_c_v', {c_v_i}, {RP.commitments}, Delta', Delta_A', rt)
        gamma_prime = hash::hashConcat(oss.str(), Delta_A_prime, broadcast.rt);
    } else {
        // bar_c_k' =  g^Rk * h^Rrk * c_k^{-gamma}
        std::string bar_c_k_prime = intcom.mul(
            intcom.mul(intcom.g_pow(broadcast.Rk), intcom.h_pow(broadcast.Rrk)),
            intcom.pow(broadcast.ck, -gamma_ZZ));

        // Delta_A' = R_s - R_v - P_A * R_k
        NTL::ZZ Delta_A_prime = Rs_minus_Rv - parties_.PA * broadcast.Rk;
        oss << broadcast.ck << bar_c_k_prime;
        // gamma' ?= H(c_s, bar_c_s', bar_c_v', {c_v_i}, {RP.commitments}, Delta', Delta_A', rt)
        gamma_prime = hash::hashConcat(oss.str(), Delta_A_prime, broadcast.rt);
    }

    return gamma_prime == broadcast.gamma;
    
    
}



NTL::ZZ VWSS::generate_party_modulus(long bits, const std::set<NTL::ZZ>& used_primes) const {
    while(true) {
        NTL::ZZ p = NTL::RandomPrime_ZZ(bits);
        if (used_primes.count(p) == 0) {
            return p;
        }
    }
}

/*
    Build the set A and B.
    A contains all parties with weight < lambda and enough parties with weight >= lambda to reach the reconstruction threshold T.
    B contains only parties with weight >= lambda.
*/
void VWSS::build_AB_small() {
    long weights_in_A = 0;
    NTL::ZZ PA(1);
    // absorb all parties with weight < lambda into A, and compute the total weight in A
    // NOTE: this assumes the weights are sorted in ascending order
    parties_.A.clear();
    parties_.B.clear();

    for (const auto& p: parties_.users) {
        if (p.weight < params_.lambda) {
            parties_.A.push_back(p);
            weights_in_A += p.weight;
            PA *= p.modulus;
        } else {
            parties_.B.push_back(p);
        }
    }

    // if the total weight in A is still smaller than the reconstruction threshold T, 
    // then we keep adding parties with weight >= lambda from B to A until reaching T
    // we reverse B to add parties with larger weight first, so that we can minimize the number of parties in A 
    while (weights_in_A < params_.T && !parties_.B.empty()) {
        parties_.A.push_back(parties_.B.back());
        weights_in_A += parties_.B.back().weight;
        PA *= parties_.B.back().modulus;
        parties_.B.pop_back();
    }

    parties_.total_weights_in_A = weights_in_A;
    parties_.PA = PA;

}

/*
    A contains parties with largest weight until the total weight in A reaches the reconstruction threshold T.
    If A contains any member with weight < lambda, then we will use UPoM with soundness boost. 
    B now may contain parties with weight < lambda, and we use UPoM with soundness boost for those parties. 
*/
void VWSS::build_AB_large() {
    long weights_in_A = 0;
    NTL::ZZ PA(1);

    parties_.A.clear();
    parties_.B.clear();

    // weights are assumed sorted in ascending order so iterate from the back to add parties with larger weight to A
    for (auto it = parties_.users.rbegin(); it != parties_.users.rend(); ++it) {
        const auto& user = *it; 
        if (weights_in_A < params_.T) {
            parties_.A.push_back(user);
            weights_in_A += user.weight;
            PA *= user.modulus;
        } else {
            parties_.B.push_back(user);
        }
    }
  
    // restore to ascending order for later use
    std::reverse(parties_.A.begin(), parties_.A.end());
    std::reverse(parties_.B.begin(), parties_.B.end());

    parties_.total_weights_in_A = weights_in_A;
    parties_.PA = PA;
}

const std::vector<long> VWSS::broadcast_size(const Broadcast& broadcast) const {
    const IntCom& intcom = intcom_;
    std::vector<long> bits_vec;
    long bits = 0;
    bits += NTL::NumBits(broadcast.Rs);
    bits += NTL::NumBits(broadcast.Rrs);
    bits += NTL::NumBits(broadcast.Rv);
    bits += NTL::NumBits(broadcast.Rrv);
    bits += intcom.commitment_bitlength(broadcast.cs);

    
    if (broadcast.Rk != 0) {
        bits += intcom.commitment_bitlength(broadcast.ck);
        bits += NTL::NumBits(broadcast.Rk);
        bits += NTL::NumBits(broadcast.Rrk);
    }

    // range proof: commitments c_1,c_2,c_3, response pairs (z_i,t_i), and tau
    // (RP.Delta is not counted: it's never transmitted, only recomputed by the verifier)
    for (const auto& c : broadcast.RP.commitments) {
        bits += intcom.commitment_bitlength(c);
    }
    for (const auto& resp : broadcast.RP.responses) {
        bits += NTL::NumBits(resp.first);
        bits += NTL::NumBits(resp.second);
    }
    bits += NTL::NumBits(broadcast.RP.tau);

    bits += 512; // these are 256 bits rt and 256 bits gamma
    bits_vec.push_back(bits);
    for (const auto& c : broadcast.cv) {
        bits += intcom.commitment_bitlength(c);
    }
    bits_vec.push_back(bits);
    return bits_vec;
}

const std::vector<long> VWSS::msg_A_size(const std::vector<MSG_A>& msg_A) const {
    std::vector<long> bits_vec;
    long bits = 0;
    for (const auto& msg : msg_A) {
        bits = 0;
        bits += NTL::NumBits(msg.vi);
        bits += NTL::NumBits(msg.rvi);
        bits_vec.push_back(bits);
    }
    return bits_vec;
}

const std::vector<long> VWSS::msg_B1_size(const std::vector<MSG_B1>& msg_B1) const {
    std::vector<long> bits_vec;
    for (const auto& msg : msg_B1) {
        long bits = 0;
        bits += NTL::NumBits(msg.vj);
        // each auth path step is a 256-bit sibling hash plus a 1-bit left/right flag
        bits += static_cast<long>(msg.pathj.size()) * (hash::HASH_SIZE * 8 + 1);
        bits_vec.push_back(bits);
    }
    return bits_vec;
}

const std::vector<long> VWSS::msg_B2_size(const std::vector<MSG_B2>& msg_B2) const {
    std::vector<long> bits_vec;
    const IntCom& intcom = intcom_;
    for (const auto& msg : msg_B2) {
        long bits = 0;
        bits += NTL::NumBits(msg.vj);
        bits += NTL::NumBits(msg.Rkj);
        bits += NTL::NumBits(msg.Rrkj);
        bits += intcom.commitment_bitlength(msg.ckj);
        // each auth path step is a 256-bit sibling hash plus a 1-bit left/right flag
        bits += static_cast<long>(msg.pathj.size()) * (hash::HASH_SIZE * 8 + 1);
        bits_vec.push_back(bits);
    }
    return bits_vec;
}
