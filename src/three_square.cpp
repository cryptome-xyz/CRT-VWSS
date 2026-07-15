#include "three_square.hpp"
#include <cmath>

namespace three_square {
    namespace {
        // Small odd primes for fast composite rejection before Miller-Rabin
        static constexpr long SMALL_PRIMES[] = {
            3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
            59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
            127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
            191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251
        };

        bool is_square(const NTL::ZZ& n, NTL::ZZ& r) {
            if (n < 0) return false;
            r = NTL::SqrRoot(n);
            return r * r == n;
        }

        // Rejects ~87% of composites before the expensive ProbPrime call
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
