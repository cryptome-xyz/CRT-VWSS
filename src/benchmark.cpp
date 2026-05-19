#include "vwss.hpp"
#include <iostream>
#include <random>
#include "three_square.hpp"

int main() {
    /*
            // initialize PARI with 128 MB stack and 2 preallocated GENs
    pari_init(128 * 1024 * 1024, 2);
    
    VWSS::Params params;
    params.lambda = 128;
    params.p0 = NTL::GenPrime_ZZ(128);
    // params.weights = {16,16,32,32,32,64,64,64,128,128,128,128,192,192,192,192};
    //params.weights = {16, 32, 48, 64};
    //params.weights = {16, 32, 48, 64, 128, 192, 256, 512};
    //params.weights = {16, 32, 48, 64, 128, 192, 256, 512, 1024, 2048, 4096, 8192};
    // Random number generator
    int a = 16; // minimum weight
    int b = 512; // maximum weight
    int t = 128; // length
    std::random_device rd;              // seed
    std::mt19937 gen(rd());             // Mersenne Twister
    std::uniform_int_distribution<> dist(a, b);
    std::vector<long> weights_temp(t);
    
    for(int i = 0; i < t; ++i) {
        weights_temp[i] = dist(gen);
    }

     // Sort ascending
    std::sort(weights_temp.begin(), weights_temp.end());

    params.weights = weights_temp;

    for (int x : params.weights) {
        std::cout << x << " ";
    }
    std::cout << std::endl;


    long total_weights(0);
    for (const auto i : params.weights) {
        total_weights += i;
    }
    std::cout << "Total Weights: " <<total_weights << std::endl;
    params.T = total_weights / 2;  // reconstruction threshold is set to half of total weights
    std::cout << "Reconstruction Threshold: " <<params.T << std::endl;
    params.t = params.T - 3 * params.lambda;    // privacy threshold t
    params.L = params.t + params.lambda;

    VWSS vwss;
    vwss.set_params(params);
    vwss.share();

    pari_close();
    */

    NTL::ZZ n = NTL::power2_ZZ(18000) * 4 + 1;
    std::vector<NTL::ZZ> decompose(3);
    auto start = std::chrono::high_resolution_clock::now();
    three_square::decompose(n, decompose);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken for three-square decomposition: " << elapsed.count() << " ms\n";

    if (three_square::verify(n, decompose)) {
        std::cout << "Decomposition successful!" << std::endl;
        std::cout << "x: " << decompose[0] << ", y: " << decompose[1] << ", z: " << decompose[2] << std::endl;
    } else {
        std::cout << "Decomposition failed!" << std::endl;
    }


    return 0;
}
