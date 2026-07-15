#include "vwss.hpp"


// Constructors
VWSS::VWSS() = default;

VWSS::VWSS(const Params& params){
    setup_dealer_and_parties(params);
}

void VWSS::setup_dealer_and_parties(const Params& params) {
    params_ = params;
    setup_dealer();
    setup_parties();
}

void VWSS::setup_dealer() {
    NTL::ZZ_p::init(params_.p0);
    dealer_.secret = NTL::random_ZZ_p();
    NTL::ZZ u = NTL::RandomBits_ZZ(params_.L);
    // lifted secret = s + u * p0
    dealer_.lifted_secret = NTL::rep(dealer_.secret) + u * params_.p0;
}

void VWSS::setup_parties() {
    // In CRT-based SS, primes must be distinct
    std::set<NTL::ZZ> used_primes;
    // initialize the set of used primes with p0 
    used_primes.insert(params_.p0);
    // initialize the product of all parties' moduli as 1
    parties_.P = 1;
    parties_.users.clear();
    parties_.users.reserve(params_.weights.size());


    for(std::size_t i = 0; i < params_.weights.size(); ++i) {
        Party p;
        p.id = static_cast<long>(i);    // set party id
        p.weight = params_.weights[i];  // set party weight
        // generate party modulus with p.weight bits
        p.modulus = generate_party_modulus(p.weight, used_primes);  
        parties_.P *= p.modulus;    // update the product of all parties' moduli
        parties_.bit_P += p.weight;
        p.share = dealer_.lifted_secret % p.modulus;    // compute party share = lifted_secret mod party modulus
        used_primes.insert(p.modulus);
        parties_.users.push_back(p);    // add party to the list of users
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

    //build_AB_small(); // method 1
    build_AB_large(); // method 2
}

const VWSS::Parties& VWSS::get_parties() const {
    return parties_;
}

const VWSS::Party& VWSS::get_party(long id) const {
    // party ids are assigned as their index into parties_.users (see setup_parties),
    // so this can be a direct lookup instead of a linear scan.
    if (id < 0 || static_cast<std::size_t>(id) >= parties_.users.size()) {
        throw std::out_of_range("no party with the given id");
    }
    return parties_.users[id];
}


VWSS::Msg VWSS::share() const { 
    // setup
    NTL::ZZ v(0); // v = \sum_{i \in A} \lambda_i * v_i
    NTL::ZZ PA(1);  // PA = \prod_{i \in A} p_i
    NTL::ZZ rv(0); // rv = \sum_{i \in A} \lambda_i * r_vi
    NTL::ZZ bound_rv(0); // bound_rv = P/2 * 2^{U + \lambda} * |A|
    VWSS::Msg msgs;

    IntCom intcom;
    const long rnd_bit_length = intcom.U() + params_.lambda;
    const long E = parties_.bit_P - 1 + parties_.A.back().weight + ceil(log2(parties_.A.size())); // upper bound (in bit) on v
    const long Sv = parties_.bit_P - 1 + rnd_bit_length + ceil(log2(parties_.A.size())); // upper bound (in bit) on rv

    // commit to dealer's lifted secret s
    IntCom::Com c_s = intcom.commit(dealer_.lifted_secret, rnd_bit_length);
    msgs.broadcast.cs = c_s.c_x;

    // commit to parties in A's shares
    msgs.broadcast.cv.reserve(parties_.A.size());
    msgs.msg_A.reserve(parties_.A.size());



    // sum of weights in A
    long sum_weights_A = 0;
    for(const auto&p :parties_.A) {
        // compute v = \sum_{i \in A} \lambda_i * v_i
        v += p.global_CRT_coef * p.share;
        PA *= p.modulus;
        IntCom::Com c_vi = intcom.commit(p.share, rnd_bit_length);
        msgs.broadcast.cv.push_back(c_vi.c_x);
        // compute rv = sum_{ i\in A} \lambda_i * r_vi
        rv += p.global_CRT_coef * c_vi.r;
        sum_weights_A += p.weight;

        VWSS::MSG_A msg;
        msg.vi = p.share;
        msg.rvi = c_vi.r;
        msg.id = p.id;
        msg.msg_id = msgs.msg_A.size();  
        msgs.msg_A.push_back(msg);
    }

    // commit to bar_s
    // first compute bar_s bit length upper bound (B')
    const long B_prime = std::max(params_.T, sum_weights_A + params_.t);
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

    // Range proof (RP_CG): prove dealer_.lifted_secret (= x_0) lies in [0, 2^T] via
    // 4*x_0*(B - x_0) + 1 = x_1^2 + x_2^2 + x_3^2, B = 2^T.
    // Slot 0 (x_0) reuses the sigma-protocol pieces already computed above instead of
    // generating fresh ones: c_0 = c_s, r_0 = c_s.r, m_0 = bar_s, s_0 = c_bar_s.r,
    // d_0 = c_bar_s.c_x, and (z_0, t_0) = (Rs, Rrs) computed below.
    // Only slots 1..3 (the three square-root witnesses) need fresh commitments/masks.
    const NTL::ZZ B_bound = NTL::power2_ZZ(params_.T);
    const NTL::ZZ& x0 = dealer_.lifted_secret;
    std::vector<NTL::ZZ> rp_x(3); // x_1, x_2, x_3
    three_square::decompose(4 * x0 * (B_bound - x0) + 1, rp_x);

    const long BCL_bit_length = params_.T + 2 * params_.lambda;  // range of m_i, i = 1..3
    const long CLS_bit_length = intcom.U() + 3 * params_.lambda; // range of s_i, i = 1..3
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

        //std::cout << "d_" << i+1 << ": " << rp_d[i] << std::endl;
    }

    // sigma <- [0, 4*B*C*L*S), d = h^sigma * c_0^{-4*m_0} * prod_{i=1..3} c_i^{-m_i}, m_0 = bar_s
    const long sigma_bit_length = params_.T + intcom.U() + 3 * params_.lambda + 2;
    const NTL::ZZ rp_sigma = NTL::RandomBits_ZZ(sigma_bit_length);
    std::string rp_d_elem = intcom.commit_elem(NTL::ZZ(0), rp_sigma); // h^sigma
    rp_d_elem = intcom.mul(rp_d_elem, intcom.pow(c_s.c_x, -4 * bar_s));
    for (int i = 0; i < 3; ++i) {
        rp_d_elem = intcom.mul(rp_d_elem, intcom.pow(rp_c[i], -rp_m[i]));
    }

    std::ostringstream oss;
    std::ostringstream rp_oss;
    rp_oss << c_bar_s.c_x; // d_0, reusing the existing mask commitment to bar_s
    for (int i = 0; i < 3; i++) {
        rp_oss << rp_d[i];
    }
    rp_oss << rp_d_elem;
    const hash::Hash rp_Delta = hash::hash(rp_oss.str());

    msgs.broadcast.RP.commitments = rp_c;
    //msgs.broadcast.RP.Delta = rp_Delta;

    NTL::ZZ gamma_ZZ;
    NTL::ZZ bar_k(0);
    NTL::ZZ k(0);
    IntCom::Com c_k;
    IntCom::Com bar_c_k;

    NTL::ZZ delta;

    for(const auto& c : msgs.broadcast.cv) {
        oss << c;
    }

    // fold the range proof's first-move commitments into the same Fiat-Shamir transcript
    for (const auto& c : rp_c) {
        oss << c;
    }
    oss.write(reinterpret_cast<const char*>(rp_Delta.data()), rp_Delta.size());
    
    // since the weights in A is sorted in ascending order, we check if the smallest weight is smaller than lambda
    // if so, we will need to compute the UPoM with soundness boost for A
    if(parties_.A[0].weight < params_.lambda) {
        k = (dealer_.lifted_secret - v) / PA;
        c_k = intcom.commit(k, rnd_bit_length);
        msgs.broadcast.ck = c_k.c_x;

        // commit to bar_k
        const long bar_k_bit_length = std::max(params_.T, E) + 1 - sum_weights_A + 2 * params_.lambda;
        bar_k = NTL::RandomBits_ZZ(bar_k_bit_length);
        bar_c_k = intcom.commit(bar_k, rnd_bit_length);

        delta = bar_s_minus_bar_v - PA * bar_k;

        oss << c_s.c_x << c_bar_s.c_x << bar_c_v.c_x << c_k.c_x << bar_c_k.c_x;
    } else {
        msgs.broadcast.ck = ""; // not used

        delta = bar_s_minus_bar_v % PA;
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

    //std::cout << "bar_c_s:" << c_bar_s.c_x << std::endl;
    //std::cout << "bar_c_v:" << bar_c_v.c_x << std::endl;
    //std::cout << "Delta_A:" << delta << std::endl;


    return msgs;
}

bool VWSS::verify_party_in_A(const MSG_A& msg, const Broadcast& broadcast) const {
    // verify broadcast message
    bool broadcast_result = verify_broadcast(broadcast);

    IntCom intcom;
    // compute commitment c_vi' = g^{vi} * h^{r_vi}
    std::string c_vi_prime = intcom.commit_elem(msg.vi, msg.rvi);

    return (c_vi_prime == broadcast.cv[msg.msg_id] && broadcast_result);

}

bool VWSS::verify_party_in_B1(const MSG_B1& msg, const NTL::ZZ& lambda_j, const NTL::ZZ& modulus,const Broadcast& broadcast) const {
    // verify broadcast message
    bool broadcast_result = verify_broadcast(broadcast);
    NTL::ZZ gamma_ZZ = hash::hashToZZ128(broadcast.gamma);

    // compute Delta_j' = (R_s - R_v - gamma * lambda_j * v_j) mod p_j
    NTL::ZZ Delta_j_prime = ((broadcast.Rs - broadcast.Rv) - gamma_ZZ * lambda_j * msg.vj) % modulus;

    hash::Hash leaf = hash::hash(Delta_j_prime);
    hash::Hash recomputed_root = hash::MerkleTree::verify(leaf, msg.pathj);

    return broadcast_result && (recomputed_root == broadcast.rt);
}

bool VWSS::verify_party_in_B2(const MSG_B2& msg, const NTL::ZZ& lambda_j, const NTL::ZZ& modulus, const Broadcast& broadcast) const {
    // verify broadcast message
    bool broadcast_result = verify_broadcast(broadcast);
    NTL::ZZ gamma_ZZ = hash::hashToZZ128(broadcast.gamma);

    IntCom intcom;
    // compute bar_c_k_j_prime = g^Rkj h^Rrkj c_k_j^-gamma
    std::string bar_c_k_j_prime = intcom.mul(intcom.commit_elem(msg.Rkj, msg.Rrkj), intcom.pow(msg.ckj,-gamma_ZZ));

   

    // compute Delta_j_prime = R_s - R_v - gamma * lambda_j v_j - p_j R_kj
    NTL::ZZ Delta_j_prime = (broadcast.Rs - broadcast.Rv) - gamma_ZZ * lambda_j * msg.vj - modulus * msg.Rkj;

    hash::Hash leaf = hash::hashConcat(msg.ckj, bar_c_k_j_prime, Delta_j_prime);
    hash::Hash recomputed_root = hash::MerkleTree::verify(leaf, msg.pathj);

    return broadcast_result && (recomputed_root == broadcast.rt);

}

bool VWSS::verify_broadcast(const Broadcast& broadcast) const {    
    // verify the UPoM proof for c_s and c_v
    NTL::ZZ gamma_ZZ = hash::hashToZZ128(broadcast.gamma);
    // bar_c_s' = g^Rs * h^Rrs * c_s^{-gamma}
    IntCom intcom;
    std::string bar_c_s_prime = intcom.mul(intcom.commit_elem(broadcast.Rs, broadcast.Rrs),
                                            intcom.pow(broadcast.cs, -gamma_ZZ));


    //std::cout << "bar_c_s_prime: " << bar_c_s_prime << std::endl;

    NTL::ZZ PA(1);
    // bar_c_v' = g^Rv * h^Rrv * (prod_{i in A} (c_vi)^{lambda_i})^{-gamma}
    std::string cv_pow_prod = intcom.commit_elem(NTL::ZZ(0), NTL::ZZ(0)); // identity
    for (size_t i = 0; i < parties_.A.size(); ++i) {
        cv_pow_prod = intcom.mul(cv_pow_prod, intcom.pow(broadcast.cv[i], parties_.A[i].global_CRT_coef));
        PA *= parties_.A[i].modulus;
    }
    std::string cv_prod = intcom.pow(cv_pow_prod, -gamma_ZZ);
    std::string bar_c_v_prime = intcom.mul(intcom.commit_elem(broadcast.Rv, broadcast.Rrv), cv_prod);

    //std::cout << "bar_c_v_prime: " << bar_c_v_prime << std::endl;

    const NTL::ZZ Rs_minus_Rv = broadcast.Rs - broadcast.Rv;


    // f_i = g^{z_i} h^{t_i} * c_i^{-gamma}, i = 1..3
    std::vector<std::string> rp_f(3);
    for (std::size_t i = 0; i < 3; ++i) {
        rp_f[i] = intcom.mul(intcom.commit_elem(broadcast.RP.responses[i].first, broadcast.RP.responses[i].second),
                              intcom.pow(broadcast.RP.commitments[i], -gamma_ZZ));
    }

    // f = h^tau * g^gamma * c_0^{4*(gamma*B - z_0)} * prod_{i=1..3} c_i^{-z_i}, B = 2^T, z_0 = Rs
    const NTL::ZZ B_bound = NTL::power2_ZZ(params_.T);
    std::string f = intcom.commit_elem(gamma_ZZ, broadcast.RP.tau); // g^gamma * h^tau
    // c_0^{4*(gamma*B - z_0)} with c_0 = c_s and z_0 = Rs
    f = intcom.mul(f, intcom.pow(broadcast.cs, 4 * (gamma_ZZ * B_bound - broadcast.Rs)));

    // c_i^{-z_i} for i = 1..3
    for (std::size_t i = 0; i < 3; ++i) {
        f = intcom.mul(f, intcom.pow(broadcast.RP.commitments[i], -broadcast.RP.responses[i].first));
    }

    // Delta' = H(f_0, f_1, f_2, f_3, f), matching the order d_0,d_1,d_2,d_3,d used in share()
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
        NTL::ZZ Delta_A_prime = Rs_minus_Rv % PA;
        // gamma' ?= H(c_s, bar_c_s', bar_c_v', {c_v_i}, {RP.commitments}, Delta', Delta_A', rt)
        gamma_prime = hash::hashConcat(oss.str(), Delta_A_prime, broadcast.rt);
    } else {
        // bar_c_k' =  g^Rk * h^Rrk * c_k^{-gamma}
        std::string bar_c_k_prime = intcom.mul(intcom.commit_elem(broadcast.Rk, broadcast.Rrk), 
                                                intcom.pow(broadcast.ck, -gamma_ZZ));

        // Delta_A' = R_s - R_v - P_A * R_k
        NTL::ZZ Delta_A_prime = Rs_minus_Rv - PA * broadcast.Rk;
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
    // absorb all parties with weight < lambda into A, and compute the total weight in A
    // NOTE: this assumes the weights are sorted in ascending order
    parties_.A.clear();
    parties_.B.clear();

    for (const auto& p: parties_.users) {
        if (p.weight < params_.lambda) {
            parties_.A.push_back(p);
            weights_in_A += p.weight;
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
        parties_.B.pop_back();
    }

}

/*
    A contains parties with largest weight until the total weight in A reaches the reconstruction threshold T.
    If A contains any member with weight < lambda, then we will use UPoM with soundness boost. 
    B now may contain parties with weight < lambda, and we use UPoM with soundness boost for those parties. 
*/
void VWSS::build_AB_large() {
    long weights_in_A = 0;

    parties_.A.clear();
    parties_.B.clear();

    // weights are assumed sorted in ascending order so iterate from the back to add parties with larger weight to A
    for (auto it = parties_.users.rbegin(); it != parties_.users.rend(); ++it) {
        const auto& user = *it; 
        if (weights_in_A < params_.T) {
            parties_.A.push_back(user);
            weights_in_A += user.weight;
        } else {
            parties_.B.push_back(user);
        }
    }
  
    // restore to ascending order for later use
    std::reverse(parties_.A.begin(), parties_.A.end());
    std::reverse(parties_.B.begin(), parties_.B.end());
}

const std::vector<long> VWSS::broadcast_size(const Broadcast& broadcast) const {
    IntCom intcom;
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
    IntCom intcom;
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
    long bits = 0;
    for (const auto& msg : msg_B1) {
        bits = 0;
        bits += NTL::NumBits(msg.vj);
        bits += log2(parties_.B.size()) * 256; // for the proof path, we assume the hash is 256 bits, so the size of the proof is log2(|B|) * 256 bits
        bits_vec.push_back(bits);
    }
    return bits_vec;
}

const std::vector<long> VWSS::msg_B2_size(const std::vector<MSG_B2>& msg_B2) const {
    std::vector<long> bits_vec;
    long bits = 0;
    IntCom intcom;
    for (const auto& msg : msg_B2) {
        bits = 0;
        bits += NTL::NumBits(msg.vj);
        bits += NTL::NumBits(msg.Rkj);
        bits += NTL::NumBits(msg.Rrkj);
        bits += intcom.commitment_bitlength(msg.ckj);
        bits += log2(parties_.B.size()) * 256;
        bits_vec.push_back(bits);
    }
    return bits_vec;
}

