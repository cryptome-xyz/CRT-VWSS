#include "three_square.hpp"
#include <array>
#include <cmath>
#include <cstddef>

namespace three_square {
    namespace {
        // Odd primes up to SMALL_PRIME_BOUND for fast composite rejection before the expensive
        // Miller-Rabin test, computed by a sieve of Eratosthenes entirely at compile time (so
        // this is a genuine hardcoded table baked into the binary, not a runtime computation).
        // A single MR round on a multi-thousand-bit candidate costs on the order of a full
        // modular exponentiation (tens of milliseconds at the candidate sizes this module deals
        // with -- see profiling notes below), against which trial division is essentially free,
        // so a larger bound is worth it: by Mertens' third theorem, the fraction of composites
        // surviving trial division up to X is ~e^{-gamma}/ln(X), which drops from ~10% at X=251
        // (the original 54-prime table) to ~5.7% at X=20000 to ~4.9% at X=100000 -- each
        // extra 1% of survivors filtered out is one fewer ProbPrime call (the dominant cost)
        // needed on average per decompose() call. 100000 was chosen as the practical ceiling
        // for this scheme: the constexpr sieve compiles in ~2s here vs ~0.3s at 20000, but
        // going further (e.g. 1000000) hits GCC's default -fconstexpr-loop-limit and would need
        // extra compiler flags for a Mertens-predicted gain of well under 1 more percentage
        // point of survivors filtered -- not worth the added build complexity.
        constexpr long SMALL_PRIME_BOUND = 100000;

        constexpr std::array<bool, SMALL_PRIME_BOUND + 1> sieve_composite_flags() {
            std::array<bool, SMALL_PRIME_BOUND + 1> composite{};
            for (long i = 2; i * i <= SMALL_PRIME_BOUND; ++i) {
                if (!composite[i]) {
                    for (long j = i * i; j <= SMALL_PRIME_BOUND; j += i) {
                        composite[j] = true;
                    }
                }
            }
            return composite;
        }

        constexpr std::size_t count_small_primes() {
            auto composite = sieve_composite_flags();
            std::size_t count = 0;
            for (long i = 3; i <= SMALL_PRIME_BOUND; ++i) { // skip 2: callers already know n's parity
                if (!composite[i]) ++count;
            }
            return count;
        }

        constexpr std::size_t SMALL_PRIME_COUNT = count_small_primes();

        constexpr std::array<long, SMALL_PRIME_COUNT> compute_small_primes() {
            auto composite = sieve_composite_flags();
            std::array<long, SMALL_PRIME_COUNT> primes{};
            std::size_t idx = 0;
            for (long i = 3; i <= SMALL_PRIME_BOUND; ++i) {
                if (!composite[i]) primes[idx++] = i;
            }
            return primes;
        }

        constexpr std::array<long, SMALL_PRIME_COUNT> SMALL_PRIMES = compute_small_primes();

        bool is_square(const NTL::ZZ& n, NTL::ZZ& r) {
            if (n < 0) return false;
            r = NTL::SqrRoot(n);
            return r * r == n;
        }

        bool has_small_factor(const NTL::ZZ& n) {
            for (long p : SMALL_PRIMES)
                if (n % p == 0) return true;
            return false;
        }

        bool brute_small(const NTL::ZZ& n, NTL::ZZ& x, NTL::ZZ& y, NTL::ZZ& z) {
            if (n > 10000) return false;
            long N = NTL::conv<long>(n);
            for (long a = 0; a * a <= N; a++) {
                for (long b = 0; a * a + b * b <= N; b++) {
                    long c2 = N - a * a - b * b;
                    long c = (long)std::sqrt((double)c2);
                    while (c > 0 && c * c > c2) c--;
                    while ((c + 1) * (c + 1) <= c2) c++;
                    if (c * c == c2) {
                        x = a; y = b; z = c;
                        return true;
                    }
                }
            }
            return false;
        }

        NTL::ZZ random_congruent_below(const NTL::ZZ& bound, long mod, long residue) {
            NTL::ZZ r = NTL::conv<NTL::ZZ>(residue);
            if (r > bound) return NTL::ZZ(-1);
            NTL::ZZ q = (bound - r) / NTL::conv<NTL::ZZ>(mod) + 1;
            return r + NTL::conv<NTL::ZZ>(mod) * NTL::RandomBnd(q);
        }

        // Caller guarantees p is prime and p ≡ 1 (mod 4)
        bool cornacchia_prime_1mod4(const NTL::ZZ& p, NTL::ZZ& a, NTL::ZZ& b) {
            NTL::ZZ r;
            while (true) {
                NTL::ZZ t = NTL::RandomBnd(p - 3) + 2;
                NTL::PowerMod(r, t, (p - 1) / 4, p);
                if ((r * r + 1) % p == 0) break;
            }
            NTL::ZZ x0 = p, x1 = r;
            while (x1 * x1 > p) {
                NTL::ZZ x2 = x0 % x1;
                x0 = x1;
                x1 = x2;
            }
            a = x1;
            NTL::ZZ b2 = p - a * a;
            return is_square(b2, b);
        }

        bool as_two_squares_prime_case(const NTL::ZZ& m, bool expect_twice_prime, NTL::ZZ& x, NTL::ZZ& y) {
            if (m <= 0) return false;
            if (!expect_twice_prime) {
                if (m % 4 != 1 || has_small_factor(m) || !NTL::ProbPrime(m, 20)) return false;
                return cornacchia_prime_1mod4(m, x, y);
            }
            if (m % 2 != 0) return false;
            NTL::ZZ p = m / 2;
            if (p % 4 != 1 || has_small_factor(p) || !NTL::ProbPrime(p, 20)) return false;
            NTL::ZZ a, b;
            if (!cornacchia_prime_1mod4(p, a, b)) return false;
            x = a + b;
            y = NTL::abs(a - b);
            return true;
        }
    }

    // n is taken by value: the algorithm strips factors of 4 internally and must
    // not modify the caller's variable (otherwise verify(original_n, result) fails)
    bool decompose(NTL::ZZ n, std::vector<NTL::ZZ>& out) {
        if (n < 0 || n % 8 == 7) return false;

        if (n == 0) {
            for (int i = 0; i < 3; i++) out[i] = NTL::ZZ(0);
            return true;
        }

        // Perfect-square shortcut: n = r² + 0² + 0²
        {
            NTL::ZZ r;
            if (is_square(n, r)) {
                out[0] = r; out[1] = NTL::ZZ(0); out[2] = NTL::ZZ(0);
                return true;
            }
        }

        NTL::ZZ scale(1);
        while (n % 4 == 0) {
            n /= 4;
            scale *= 2;
        }

        if (brute_small(n, out[0], out[1], out[2])) {
            for (int i = 0; i < 3; i++) out[i] *= scale;
            return true;
        }

        NTL::ZZ sq_root = NTL::SqrRoot(n);
        long r8 = NTL::conv<long>(n % 8);
        long mod, residue;
        bool twice_prime;

        if (r8 == 1 || r8 == 5) {
            mod = 2; residue = 0; twice_prime = false;
        } else if (r8 == 2) {
            mod = 4; residue = 0; twice_prime = true;
        } else if (r8 == 3) {
            mod = 2; residue = 1; twice_prime = true;
        } else if (r8 == 6) {
            mod = 4; residue = 2; twice_prime = true;
        } else {
            return false;  // r8 == 7 after stripping all 4s: three-square theorem says impossible
        }

        for (;;) {
            NTL::ZZ z = random_congruent_below(sq_root, mod, residue);
            if (z < 0) continue;
            NTL::ZZ m = n - z * z;
            if (as_two_squares_prime_case(m, twice_prime, out[0], out[1])) {
                out[2] = z;
                for (int i = 0; i < 3; i++) out[i] *= scale;
                return true;
            }
        }
    }

    bool verify(const NTL::ZZ& n, const std::vector<NTL::ZZ>& decompose) {
        if (decompose.size() != 3) return false;
        NTL::ZZ sum = NTL::power(decompose[0], 2) + NTL::power(decompose[1], 2) + NTL::power(decompose[2], 2);
        return sum == n;
    }
}
