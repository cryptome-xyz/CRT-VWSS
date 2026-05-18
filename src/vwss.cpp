#include "vwss.hpp"


// Constructors
VWSS::VWSS() = default;

VWSS::VWSS(const Params& params){
    set_params(params);
}

void VWSS::set_params(const Params& params) {
    // TODO: add validator of params
    params_ = params;
    setup_dealer();
    setup_parties();
}

void VWSS::setup_dealer() {
    NTL::ZZ_p::init(params_.p0);
    dealer_.secret = NTL::random_ZZ_p();
    NTL::ZZ u = NTL::RandomBits_ZZ(params_.L);
    dealer_.lifted_secret = NTL::rep(dealer_.secret) + u * params_.p0;
}

void VWSS::setup_parties() {
    // In CRT-based SS, primes must be distinct
    std::set<NTL::ZZ> used_primes;
    // initialize the set of used primes with p0 
    used_primes.insert(params_.p0);
    // initialize the product of all parties' moduli as 1
    parties_.P = 1;


    for(std::size_t i = 0; i < params_.weights.size(); ++i) {
        Party p;
        p.id = static_cast<long>(i);
        p.weight = params_.weights[i];
        p.modulus = generate_party_modulus(p.weight, used_primes);
        parties_.P *= p.modulus;    // update the product of all parties' moduli
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

    //build_AB_small();
    build_AB_large();
}

const VWSS::Parties& VWSS::get_parties() const {
    return parties_;
}


void VWSS::share() const { 
    auto start = std::chrono::high_resolution_clock::now();

    
    // setup
    NTL::ZZ v(0); // v = \sum_{i \in A} \lambda_i * v_i
    NTL::ZZ PA(1);  // PA = \prod_{i \in A} p_i
    NTL::ZZ rv(0); // rv = \sum_{i \in A} \lambda_i * r_vi
    VWSS::Msg msgs;

    IntCom intcom;
    long rnd_bit_length = intcom.U() + params_.lambda;

    // commit to dealer's lifted secret s
    IntCom::Com c_s = intcom.commit(dealer_.lifted_secret, rnd_bit_length);
    msgs.broadcast.cs = c_s.c_x;

    // commit to parties in A's shares
    msgs.broadcast.cv.reserve(parties_.A.size());
    for(const auto&p :parties_.A) {
        // compute v = \sum_{i \in A} \lambda_i * v_i (mod P)
        v += p.global_CRT_coef * p.share;
        //v = NTL::AddMod(v, NTL::MulMod(p.global_CRT_coef, p.share, parties_.P), parties_.P);
        PA *= p.modulus;
        // NTL::ZZ cv_rand = NTL::RandomBnd(bound);
        // rM += p.global_CRT_coef * cv_rand;
        IntCom::Com c_vi = intcom.commit(p.share, rnd_bit_length);
        msgs.broadcast.cv.push_back(c_vi.c_x);
        // compute rv = sum_{ i\in A} \lambda_i * r_vi
        rv += p.global_CRT_coef * c_vi.r;
        
        VWSS::MSG_A msg;
        msg.vi = p.share;
        msg.rvi = c_vi.r;
        msg.id = p.id;
        msgs.msg_A.push_back(msg);
    }

    // commit to bar_s
    // here, we assume p0 has lambda bits, so s has length L + lmabda, and we add 2lambda bits for hiding
    long bar_s_bit_length = params_.L + 3 * params_.lambda ; 
    NTL::ZZ bar_s = NTL::RandomBits_ZZ(bar_s_bit_length);
    IntCom::Com c_bar_s = intcom.commit(bar_s, rnd_bit_length);

    // commit to bar_v
    long bar_v_bit_length = NTL::NumBits(parties_.P) + 2 * params_.lambda;
    NTL::ZZ bar_v = NTL::RandomBits_ZZ(bar_v_bit_length);
    long bar_v_rand_bit_length = rnd_bit_length * NTL::NumBits(parties_.P/2) * parties_.A.size(); 
    IntCom::Com bar_c_v = intcom.commit(bar_v, rnd_bit_length);

    NTL::ZZ delta = (bar_s - bar_v) % PA;

    long temp_range = std::max(NTL::NumBits(parties_.P), params_.L + params_.lambda);

    std::vector<hash::Hash> hash_list;


    for(const auto&p : parties_.B) {
        if(p.modulus >= params_.lambda) {
            VWSS::MSG_B1 msg_j;
            msg_j.id = p.id;
            NTL::ZZ Delta_j = (bar_s - bar_v) % p.modulus;
            msg_j.Deltaj = Delta_j;
            msg_j.vj = p.share;
            msgs.msg_B1.push_back(msg_j);
            // hash the Delta_j
            hash_list.push_back(hash::hash(Delta_j));

            
        } else {
            VWSS::MSG_B2 msg_j;
            msg_j.id = p.id;
            // generate kj 
            NTL::ZZ kj = (dealer_.lifted_secret - v - p.global_CRT_coef * p.share) % p.modulus;
            // commit to kj
            IntCom::Com c_kj = intcom.commit(kj, rnd_bit_length);
            msg_j.ckj = c_kj.c_x;
            // compute bar_kj 
            long bar_kj_bit_length = NTL::NumBits(temp_range - p.weight) + NTL::NumBits(parties_.P/2) + 2 * params_.lambda;
            NTL::ZZ bar_kj = NTL::RandomBits_ZZ(bar_kj_bit_length);
            IntCom::Com c_bar_kj = intcom.commit(bar_kj, rnd_bit_length);
            NTL::ZZ Delta_j = bar_s - bar_v - p.global_CRT_coef * bar_kj;
            msg_j.Deltaj = Delta_j;
            msgs.msg_B2.push_back(msg_j);
            hash_list.push_back(hash::hashConcat(c_kj.c_x, c_bar_kj.c_x, Delta_j));
        }
        
    }
    
    if(!hash_list.empty()) {
        hash::MerkleTree tree(hash_list);
        msgs.broadcast.rt = tree.root();
    }
    

   

    NTL::ZZ gamma_ZZ;
    NTL::ZZ bar_k;
    NTL::ZZ k;
    IntCom::Com c_k;
    IntCom::Com bar_c_k;

    // since the weights in A is sorted in ascending order, we check if the smallest weight is smaller than lambda
    // if so, we will need to compute the UPoM with soundness boost for A
    if(parties_.A[0].weight < params_.lambda) {
        delta = 0;
        k = (dealer_.lifted_secret - v) / PA;
        IntCom::Com c_k = intcom.commit(k, rnd_bit_length);
        msgs.broadcast.ck = c_k.c_x;

        // commit to bar_k
        long bar_k_bit_length = NTL::NumBits(temp_range) + 2 * params_.lambda;
        bar_k = NTL::RandomBits_ZZ(bar_k_bit_length);
        bar_c_k = intcom.commit(bar_k, rnd_bit_length);

        // overwrite delta for UPoM with soundness boost
        delta = bar_s - bar_v - PA * bar_k;

        std::ostringstream oss;
        oss << c_s.c_x << c_bar_s.c_x << bar_c_v.c_x << c_k.c_x << bar_c_k.c_x;
        for(const auto& c : msgs.broadcast.cv) {
            oss << c;
        }
        msgs.broadcast.gamma = hash::hashConcat(oss.str(), delta, msgs.broadcast.rt);
        gamma_ZZ = hash::hashToZZ128(msgs.broadcast.gamma);

    } else {
        msgs.broadcast.ck = ""; // not used

        std::ostringstream oss;
        oss << c_s.c_x << c_bar_s.c_x << bar_c_v.c_x;
        for(const auto& c : msgs.broadcast.cv) {
            oss << c;
        }
        msgs.broadcast.gamma = hash::hashConcat(oss.str(), delta, msgs.broadcast.rt);
        gamma_ZZ = hash::hashToZZ128(msgs.broadcast.gamma);

    }

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

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken for sharing: " << elapsed.count() << " ms\n";

    std::vector<long> broadcast_bit_length = broadcast_size(msgs.broadcast);

    std::cout << "Broadcast Message Bit Length (without user commitmenets): " << broadcast_bit_length[0] << std::endl;
    std::cout << "Broadcast Message Bit Length (with user commitmenets): " << broadcast_bit_length[1] << std::endl;

    std::cout << "A Size: " << parties_.A.size() << std::endl;
    std::cout << "B Size: " << parties_.B.size() << std::endl;
}


void VWSS::verify(const Broadcast& broadcast) const {
    auto start = std::chrono::high_resolution_clock::now();



    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken for verify: " << elapsed.count() << " ms\n";
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
    A contains parties with weight >= T and parties with weight < lambda
    B contains only parties with weight >= lambda.
    This is one way to build A and B for comparing the performance with another way.
*/

void VWSS::build_AB_small() {
    long weights_in_A = 0;
    // absorb all parties with weight < lambda into A, and compute the total weight in A
    // this assumes the weights are sorted in ascending order
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
    for (auto i = parties_.users.rbegin(); i != parties_.users.rend(); ++i) {
        if (weights_in_A < params_.T) {
            parties_.A.push_back(*i);
            weights_in_A += i->weight;
        } else {
            parties_.B.push_back(*i);
        }
    }

    // sort to ascending order for later use
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
        bits += NTL::NumBits(msg.vj);
        bits += NTL::NumBits(msg.Deltaj);
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
        bits += NTL::NumBits(msg.vj);
        bits += NTL::NumBits(msg.Deltaj);
        bits += NTL::NumBits(msg.Rkj);
        bits += NTL::NumBits(msg.Rrkj);
        bits += intcom.commitment_bitlength(msg.ckj);
        bits += log2(parties_.B.size()) * 256;
    }
    return bits_vec;
}

