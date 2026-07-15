#include "vwss.hpp"
#include <iostream>
#include <random>
#include "three_square.hpp"

std::vector<long> generate_random_weights(int num_parties, int min_weight, int max_weight) {
    std::random_device rd;              // seed
    std::mt19937 gen(rd());             // Mersenne Twister
    std::uniform_int_distribution<> dist(min_weight, max_weight);
    std::vector<long> weights(num_parties);
    
    for(int i = 0; i < num_parties; ++i) {
         weights[i] = dist(gen);
    }

    // Sort ascending
    std::sort(weights.begin(), weights.end());
    return weights;
}


int main() {
    // initialize PARI with 128 MB stack and 2 preallocated GENs
    pari_init(128 * 1024 * 1024, 2);
    
    VWSS::Params params;
    params.lambda = 128;
    params.p0 = NTL::GenPrime_ZZ(128);  
    
    // generate random weights for 128 parties in the range [10, 255]
    // params.weights = generate_random_weights(128, 10, 255);
    
    // For reproducibility, the test vectors used in the paper are recorded here.
    // test vector randomly generated for 4 parties
    // params.weights = {29, 43, 212, 240};
    // test vector randomly generated for 8 parties
    params.weights = {42, 47, 50, 61, 124, 166, 201, 245};
    //params.weights = {128, 128, 128, 128, 128, 166, 201, 245};
    // test vector randomly generated for 16 parties
    //params.weights = {18, 60, 72, 87, 106, 122, 164, 168, 177, 183, 200, 207, 216, 231, 235, 244};
    // test vector randomly generated for 32 parties
    //params.weights = {19,23,26,41,53,53,55,64,68,72,75,81,86,93,96,100,115,118,127,154,157,173,186,202,208,209,211,217,218,227,229,253};
    // test vector randomly generated for 64 parties
    // params.weights = {14, 19, 21, 23, 28, 31, 31, 31, 31, 33, 38, 39, 51, 52, 59, 62,
    //     64, 66, 80, 82, 82, 84, 84, 102, 105, 110, 110, 120, 132, 145, 150, 153, 155,
    //     158, 159, 167, 168, 174, 174, 177, 178, 179, 179, 189, 193, 202, 209, 209,
    //     212, 214, 219, 219, 223, 227, 228, 229, 230, 237, 239, 245, 248, 253, 254, 255};
    // test vector randomly generated for 128 parties
    // params.weights = {10, 11, 11, 13, 15, 20, 21, 26, 27, 32, 34, 35, 38, 39, 40, 
    //     42, 42, 43, 44, 45, 47, 48, 53, 59, 61, 63, 67, 69, 73, 73, 77, 78, 83, 90, 90,
    //      90, 94,97, 99, 104, 107, 107, 108, 113, 116, 119, 119, 120, 121, 123, 125,
    //       128, 129, 133, 133, 134, 135, 138, 138, 139, 140, 141, 144, 145, 146, 147, 
    //       149, 149, 150, 150, 153, 154 ,155, 159, 161, 162, 163, 166, 168, 168, 174, 
    //       176, 180, 181, 182, 185, 186, 189, 198 ,201, 203, 203, 204, 204, 206, 208, 
    //       208, 208, 209, 209, 209, 217, 219, 219, 222, 224, 226, 228, 232, 233, 233, 
    //       237, 237, 239, 240, 242, 242, 242, 243, 246, 248, 250, 250, 251, 253, 253, 254, 255};

    // A ETH-based test vector from Shehata et al.'s paper
    // params.weights = {10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12,12,13,13,13,13,14,14,14,15,15,16,17,18,19,19,
    // 20,22,23,26,30,33,35,38,40,48,49,54,65,67,73,73,102,109,123,150,155,169,178,208,749,1611};
    
    long total_weights(0);
    for (const auto i : params.weights) {
        total_weights += i;
    }
    std::cout << "Total Weights: " <<total_weights << std::endl;
    params.T = total_weights * 0.51;  // reconstruction threshold is set to more than half of total weights
    std::cout << "Reconstruction Threshold: " << params.T << std::endl;
    params.t = params.T - 3 * params.lambda;    // privacy threshold t
    std::cout << "Privacy Threshold: " << params.t << std::endl;
    params.L = params.t + 2 * params.lambda;    // L = t + 2 * lambda
    std::cout << "Lift Parameter: " << params.L << std::endl;


    VWSS vwss;
    vwss.setup_dealer_and_parties(params);
    auto start = std::chrono::high_resolution_clock::now();
    VWSS::Msg msgs = vwss.share();
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken for sharing: " << elapsed.count() << " ms\n";

    std::vector<long> broadcast_bit_length = vwss.broadcast_size(msgs.broadcast);

    std::cout << "Broadcast Message Bit Length (without user commitmenets): " << broadcast_bit_length[0] << std::endl;
    std::cout << "Broadcast Message Bit Length (with user commitmenets): " << broadcast_bit_length[1] << std::endl;

    long total_msg_bit_length = 0;
    std::vector<long> msg_A_bit_length = vwss.msg_A_size(msgs.msg_A);
    for (size_t i = 0; i < msg_A_bit_length.size(); i++) {
        // std::cout << "MSG_A Bit Length for party " << msgs.msg_A[i].id << ": " << msg_A_bit_length[i] << std::endl;
        total_msg_bit_length += msg_A_bit_length[i];
    } 
    

    std::vector<long> msg_B1_bit_length = vwss.msg_B1_size(msgs.msg_B1);
    for(size_t i = 0; i < msg_B1_bit_length.size(); i++) {
        //std::cout << "MSG_B1 Bit Length for party " << msgs.msg_B1[i].id << ": " << msg_B1_bit_length[i] << std::endl;
        total_msg_bit_length += msg_B1_bit_length[i];
    }

    std::vector<long> msg_B2_bit_length = vwss.msg_B2_size(msgs.msg_B2);
    for(size_t i = 0; i < msg_B2_bit_length.size(); i++) {
        //std::cout << "MSG_B2 Bit Length for party " << msgs.msg_B2[i].id << ": " << msg_B2_bit_length[i] << std::endl;
        total_msg_bit_length += msg_B2_bit_length[i];
    }

    std::cout << "Total Private Message Bit Length: " << total_msg_bit_length << std::endl;
    std::cout << "Average Private Message Bit Length: " << total_msg_bit_length / (msg_A_bit_length.size() + msg_B1_bit_length.size() + msg_B2_bit_length.size()) << std::endl;

    std::chrono::milliseconds total_verify_A_time{0};
    for (const auto& a : msgs.msg_A) {
        auto v_start = std::chrono::high_resolution_clock::now();
        bool result = vwss.verify_party_in_A(a, msgs.broadcast);
        auto v_end = std::chrono::high_resolution_clock::now();
        total_verify_A_time += std::chrono::duration_cast<std::chrono::milliseconds>(v_end - v_start);
        std::cout << "Message verification " << (result ? "successful" : "NOT successful")
                   << " for user" << a.id << " in A!" << std::endl;
    }

    std::chrono::milliseconds total_verify_B1_time{0};
    for (const auto& b1 : msgs.msg_B1) {
        const VWSS::Party& p = vwss.get_party(b1.id);
        auto v_start = std::chrono::high_resolution_clock::now();
        bool result = vwss.verify_party_in_B1(b1, p.global_CRT_coef, p.modulus, msgs.broadcast);
        auto v_end = std::chrono::high_resolution_clock::now();
        total_verify_B1_time += std::chrono::duration_cast<std::chrono::milliseconds>(v_end - v_start);
        std::cout << "Message verification " << (result ? "successful" : "NOT successful")
                   << " for user" << p.id << " in B1!" << std::endl;
    }

    std::chrono::milliseconds total_verify_B2_time{0};
    for (const auto& b2 : msgs.msg_B2) {
        const VWSS::Party& p = vwss.get_party(b2.id);
        auto v_start = std::chrono::high_resolution_clock::now();
        bool result = vwss.verify_party_in_B2(b2, p.global_CRT_coef, p.modulus, msgs.broadcast);
        auto v_end = std::chrono::high_resolution_clock::now();
        total_verify_B2_time += std::chrono::duration_cast<std::chrono::milliseconds>(v_end - v_start);
        std::cout << "Message verification " << (result ? "successful" : "NOT successful")
                   << " for user" << p.id << " in B2!" << std::endl;
    }

    if (!msgs.msg_A.empty()) {
        std::cout << "Average verification time for A: "
                   << total_verify_A_time.count() / static_cast<double>(msgs.msg_A.size()) << " ms\n";
    }
    if (!msgs.msg_B1.empty()) {
        std::cout << "Average verification time for B1: "
                   << total_verify_B1_time.count() / static_cast<double>(msgs.msg_B1.size()) << " ms\n";
    }
    if (!msgs.msg_B2.empty()) {
        std::cout << "Average verification time for B2: "
                   << total_verify_B2_time.count() / static_cast<double>(msgs.msg_B2.size()) << " ms\n";
    }

    long total_verify_count = static_cast<long>(msgs.msg_A.size() + msgs.msg_B1.size() + msgs.msg_B2.size());
    if (total_verify_count > 0) {
        auto total_verify_time = total_verify_A_time + total_verify_B1_time + total_verify_B2_time;
        std::cout << "Average verification time overall: "
                   << total_verify_time.count() / static_cast<double>(total_verify_count) << " ms\n";
    }

    pari_close();


    // NTL::ZZ n = NTL::power2_ZZ(18000) * 4 + 1;
    // std::vector<NTL::ZZ> decompose(3);
    // auto start = std::chrono::high_resolution_clock::now();
    // three_square::decompose(n, decompose);
    // auto end = std::chrono::high_resolution_clock::now();
    // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "Time taken for three-square decomposition: " << elapsed.count() << " ms\n";

    // if (three_square::verify(n, decompose)) {
    //     std::cout << "Decomposition successful!" << std::endl;
    //     std::cout << "x: " << decompose[0] << ", y: " << decompose[1] << ", z: " << decompose[2] << std::endl;
    // } else {
    //     std::cout << "Decomposition failed!" << std::endl;
    // }


    return 0;
}
