#include "three_square.hpp"

namespace three_square {
    namespace {
        bool is_square(const NTL::ZZ& n, NTL::ZZ& r) {
            if (n < 0) return false;
            r = NTL::SqrRoot(n);
            return r * r == n;
        }

        bool brute_small(const NTL::ZZ& n, NTL::ZZ& x, NTL::ZZ& y, NTL::ZZ& z) {
            if (n <= 10000) {
                long N = NTL::conv<long>(n);
                for (long a = 0; a * a <= N; a++) {
                    for (long b = 0; a + b * b <= N; b++) {
                        long c2 = N - a * a - b * b;
                        long c = NTL::conv<long>(NTL::SqrRoot(NTL::conv<NTL::ZZ>(c2)));
                        if (c * c == c2) {
                            x = NTL::conv<NTL::ZZ>(a);
                            y = NTL::conv<NTL::ZZ>(b);
                            z = NTL::conv<NTL::ZZ>(c);
                            return true;
                        }
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
        
        bool cornacchia_prime_1mod4(const NTL::ZZ& p, NTL::ZZ& a, NTL::ZZ& b) {
            if (p <= 2 || p % 4 != 1 || !NTL::ProbPrime(p, 30)) return false;

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
                if (m % 4 != 1 || !NTL::ProbPrime(m, 30)) return false;
                return cornacchia_prime_1mod4(m, x, y);
            }

            if (m % 2 != 0) return false;
            NTL::ZZ p = m / 2;
            if (p % 4 != 1 || !NTL::ProbPrime(p, 30)) return false;

            NTL::ZZ a, b;
            if (!cornacchia_prime_1mod4(p, a, b)) return false;

            x = a + b;
            y = NTL::abs(a - b);
            return true;
        }
    }

    bool decompose(NTL::ZZ& n, std::vector<NTL::ZZ>& decompose) {
        if (n < 0 || n % 8 == 7) {
            return false;
        }

        if (n == 0) {
            for(int i = 0; i < 3; i++) {
                decompose[i] = NTL::ZZ(0);
            }
            return true;
        }

        NTL::ZZ scale(1);
        while(n % 4 == 0) {
            n = n / 4;
            scale *= 2;
        }

        if(brute_small(n,decompose[0], decompose[1], decompose[2])) {
            for(int i = 0; i < 3; i++) {
                decompose[i] *= scale;
            }
            return true;
        }

        NTL::ZZ root = NTL::SqrRoot(n);
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
            return false;
        }

        for (;;) {
            NTL::ZZ z = random_congruent_below(root, mod, residue);
            if (z < 0) continue;

            NTL::ZZ m = n - z * z;
            decompose[2] = z;
            if (as_two_squares_prime_case(m, twice_prime, decompose[0], decompose[1])) {
                decompose[0] *= scale;
                decompose[1] *= scale;
                decompose[2] *= scale;
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
