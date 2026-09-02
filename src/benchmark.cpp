#include "vwss.hpp"
#include <chrono>
#include <iostream>
#include <cstdlib>
#include <random>
#include "three_square.hpp"

using BenchmarkClock = std::chrono::steady_clock;
using FractionalMilliseconds = std::chrono::duration<double, std::milli>;

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

/*
 * Test vectors for timing reported in Table 2
 */
const VWSS::Params test_vector(int num_parties) {
    // initialize PARI with 128 MB stack and 2 preallocated GENs
    pari_init(128 * 1024 * 1024, 2);
    VWSS::Params params;
    params.lambda = 128;
    params.p0 = NTL::GenPrime_ZZ(128);
    switch (num_parties) {
        // this cannot be set to 51% of total weight since the privacy threshold t would be negative
        case 4 : {
            params.weights = {100, 132, 144, 148};
            long total_weights(0);
            for (const auto i : params.weights) {
                total_weights += i;
            }
            params.T = total_weights;
            // privacy threshold is T - 2lambda since the total weight is small in this case
            params.t = params.T - 2 * params.lambda;
            params.L = params.T - params.lambda;
            std::cout << "Total Weights: " << total_weights << std::endl;
            std::cout << "Reconstruction Threshold: " << params.T << std::endl;
            return params;
        }
        case 8 : {
            // test vector randomly generated for 8 parties (1,1,2,2,2,3,3,4), scaling factor 128
            params.weights = {35, 82, 91, 103, 129, 131, 158, 207};
            break;
        } 
        case 16 : {
            // test vector randomly generated for 16 parties
            params.weights = {18, 60, 72, 87, 106, 122, 164, 168, 177, 183, 200, 207, 216, 231, 235, 244};
            break;
        }
        case 32 : {
            // test vector randomly generated for 32 parties
            params.weights = {19,23,26,41,53,53,55,64,68,72,75,81,86,93,96,100,115,118,127,154,157,173,186,202,208,209,211,217,218,227,229,253};
            break;
        }
        case 64 : {
            // test vector randomly generated for 64 parties
            params.weights = {14, 19, 21, 23, 28, 31, 31, 31, 31, 33, 38, 39, 51, 52, 59, 62,
                64, 66, 80, 82, 82, 84, 84, 102, 105, 110, 110, 120, 132, 145, 150, 153, 155,
                158, 159, 167, 168, 174, 174, 177, 178, 179, 179, 189, 193, 202, 209, 209,
                212, 214, 219, 219, 223, 227, 228, 229, 230, 237, 239, 245, 248, 253, 254, 255};
            break;
        }
        case 128 : {
            // test vector randomly generated for 128 parties
            params.weights = {10, 11, 11, 13, 15, 20, 21, 26, 27, 32, 34, 35, 38, 39, 40, 
                42, 42, 43, 44, 45, 47, 48, 53, 59, 61, 63, 67, 69, 73, 73, 77, 78, 83, 90, 90,
                90, 94,97, 99, 104, 107, 107, 108, 113, 116, 119, 119, 120, 121, 123, 125,
                128, 129, 133, 133, 134, 135, 138, 138, 139, 140, 141, 144, 145, 146, 147, 
                149, 149, 150, 150, 153, 154 ,155, 159, 161, 162, 163, 166, 168, 168, 174, 
                176, 180, 181, 182, 185, 186, 189, 198 ,201, 203, 203, 204, 204, 206, 208, 
                208, 208, 209, 209, 209, 217, 219, 219, 222, 224, 226, 228, 232, 233, 233, 
                237, 237, 239, 240, 242, 242, 242, 243, 246, 248, 250, 250, 251, 253, 253, 254, 255};
            break;
        }
    }
    long total_weights(0);
    for (const auto i : params.weights) {
        total_weights += i;
    }
    params.T = static_cast<long>(total_weights * 0.51);
    params.t = params.T - 3 * params.lambda;
    params.L = params.t + 2 * params.lambda;
    std::cout << "Total Weights: " << total_weights << std::endl;
    std::cout << "Reconstruction Threshold: " << params.T << std::endl;
    return params;
    
}

/*
 * Test vectors for timing reported in Table 1
 * Total weightes = 4711
 * Reconstruction Threshold = 2402 (51% of total weights)
 * Privacy Threshold = 2402 - 3 * lambda = 2050
 */
const VWSS::Params test_vector_eth() {
    // initialize PARI with 128 MB stack and 2 preallocated GENs
    pari_init(128 * 1024 * 1024, 2);
    VWSS::Params params;
    params.lambda = 128;
    params.p0 = NTL::GenPrime_ZZ(128);
    params.weights = {10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12,12,
                      13,13,13,13,14,14,14,15,15,16,17,18,19,19,20,22,23,26,30,33,35,38,40,
                      48,49,54,65,67,73,73,102,109,123,150,155,169,178,208,749,1611};
    long total_weights(0);
    for (const auto i : params.weights) {
        total_weights += i;
    }
    params.T = static_cast<long>(total_weights * 0.51); 
    params.t = params.T - 3 * params.lambda;    
    params.L = params.t + 2 * params.lambda;
    std::cout << "Total Weights: " << total_weights << std::endl;
    std::cout << "Reconstruction Threshold: " << params.T << std::endl;
    return params;
}

int main() {

{
    // choose the test vector for benchmarking
    VWSS::Params params = test_vector(4);
    VWSS vwss;
    // true for VWSS-method 1 
    vwss.setup_dealer_and_parties(params, true);
    auto start = BenchmarkClock::now();
    VWSS::Msg msgs = vwss.share();
    auto end = BenchmarkClock::now();
    auto elapsed = FractionalMilliseconds(end - start);
    std::cout << "Time taken for sharing: " << elapsed.count() << " ms\n";

    std::vector<long> broadcast_bit_length = vwss.broadcast_size(msgs.broadcast);

    std::cout << "Broadcast Message Bit Length (without user commitmenets): " << broadcast_bit_length[0] << std::endl;
    std::cout << "Broadcast Message Bit Length (with user commitmenets): " << broadcast_bit_length[1] << std::endl;

    long total_msg_bit_length = 0;
    std::vector<long> msg_A_bit_length = vwss.msg_A_size(msgs.msg_A);
    for (size_t i = 0; i < msg_A_bit_length.size(); i++) {
        total_msg_bit_length += msg_A_bit_length[i];
    } 
    

    std::vector<long> msg_B1_bit_length = vwss.msg_B1_size(msgs.msg_B1);
    for(size_t i = 0; i < msg_B1_bit_length.size(); i++) {
        total_msg_bit_length += msg_B1_bit_length[i];
    }

    std::vector<long> msg_B2_bit_length = vwss.msg_B2_size(msgs.msg_B2);
    for(size_t i = 0; i < msg_B2_bit_length.size(); i++) {
        total_msg_bit_length += msg_B2_bit_length[i];
    }

    std::cout << "Total Private Message Bit Length: " << total_msg_bit_length << std::endl;
    std::cout << "Average Private Message Bit Length: " << total_msg_bit_length / (msg_A_bit_length.size() + msg_B1_bit_length.size() + msg_B2_bit_length.size()) << std::endl;

    start = BenchmarkClock::now();

    bool result = vwss.verify_broadcast(msgs.broadcast);

    end = BenchmarkClock::now();
    elapsed = FractionalMilliseconds(end - start);
    std::cout << "Time taken for verifying broadcast message: " << elapsed.count() << " ms\n";

    BenchmarkClock::duration total_verify_A_time{};
    for (const auto& a : msgs.msg_A) {
        auto v_start = BenchmarkClock::now();
        bool result = vwss.verify_party_in_A(a, msgs.broadcast);
        auto v_end = BenchmarkClock::now();
        total_verify_A_time += v_end - v_start;
        if (!result) {
            std::cerr << "Message verification FAILED for user" << a.id << " in A!" << std::endl;
            pari_close();
            std::exit(EXIT_FAILURE);
        }
    }

    BenchmarkClock::duration total_verify_B1_time{};
    for (const auto& b1 : msgs.msg_B1) {
        const VWSS::Party& p = vwss.get_party(b1.id);
        auto v_start = BenchmarkClock::now();
        bool result = vwss.verify_party_in_B1(b1, p.global_CRT_coef, p.modulus, msgs.broadcast);
        auto v_end = BenchmarkClock::now();
        total_verify_B1_time += v_end - v_start;
        if (!result) {
            std::cerr << "Message verification FAILED for user" << p.id << " in B1!" << std::endl;
            pari_close();
            std::exit(EXIT_FAILURE);
        }
    }

    BenchmarkClock::duration total_verify_B2_time{};
    for (const auto& b2 : msgs.msg_B2) {
        const VWSS::Party& p = vwss.get_party(b2.id);
        auto v_start = BenchmarkClock::now();
        bool result = vwss.verify_party_in_B2(b2, p.global_CRT_coef, p.modulus, msgs.broadcast);
        auto v_end = BenchmarkClock::now();
        total_verify_B2_time += v_end - v_start;
        if (!result) {
            std::cerr << "Message verification FAILED for user" << p.id << " in B2!" << std::endl;
            pari_close();
            std::exit(EXIT_FAILURE);
        }
    }

    if (!msgs.msg_A.empty()) {
        std::cout << "Average verification time for A: "
                   << FractionalMilliseconds(total_verify_A_time).count() /
                          static_cast<double>(msgs.msg_A.size())
                   << " ms\n";
    }
    if (!msgs.msg_B1.empty()) {
        std::cout << "Average verification time for B1: "
                   << FractionalMilliseconds(total_verify_B1_time).count() /
                          static_cast<double>(msgs.msg_B1.size())
                   << " ms\n";
    }
    if (!msgs.msg_B2.empty()) {
        std::cout << "Average verification time for B2: "
                   << FractionalMilliseconds(total_verify_B2_time).count() /
                          static_cast<double>(msgs.msg_B2.size())
                   << " ms\n";
    }

    long total_verify_count = static_cast<long>(msgs.msg_A.size() + msgs.msg_B1.size() + msgs.msg_B2.size());
    if (total_verify_count > 0) {
        auto total_verify_time = total_verify_A_time + total_verify_B1_time + total_verify_B2_time;
        std::cout << "Average verification time overall: "
                   << FractionalMilliseconds(total_verify_time).count() /
                          static_cast<double>(total_verify_count)
                   << " ms\n";
    }
}
    pari_close();

    return 0;
}
