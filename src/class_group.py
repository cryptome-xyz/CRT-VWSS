# Compute Class Group Parameters
from sage.all import *

# For prototype, 1600 bits discriminant should provide 120 bit security based on "Survey on IQ cryptography" 3.2 
# (https://github.com/ZenGo-X/class/blob/master/src/primitives/vdf.rs)
DISCRIMINANT_BITS = 1600

# Generate a negative discriminant Delta
def generate_discriminant():
    lower = Integer(1) << (DISCRIMINANT_BITS - 1)
    upper = (Integer(1) << DISCRIMINANT_BITS) - 1
    set_random_seed(8)  # for reproducibility
    while True:
        p = random_prime(upper, lbound=lower)
        if p % 4 == 3:
            Delta = -p  # Delta = 1 mod 4
            print("Delta = ", Delta)
            print("bits = ", abs(Delta).nbits())
            return Integer(Delta)


def find_split_prime(Delta, ell):
    if not is_prime(ell):
        raise ValueError("ell must be prime")  
    while True:
        if kronecker(Delta, ell) == 1:  # check if Delta is a quadratic residue mod ell 
            return ell
        ell = next_prime(ell)

def to_tuple(f):
    return tuple(Integer(x) for x in f)
         
if __name__ == "__main__":
    Delta = generate_discriminant()
    ell_g = find_split_prime(Delta, 3)
    ell_h = find_split_prime(Delta, next_prime(ell_g))

    g = pari.qfbprimeform(Delta, ell_g)
    h = pari.qfbprimeform(Delta, ell_h)

    g_tuple = to_tuple(g)
    h_tuple = to_tuple(h)

    print("g =", g_tuple)
    print("h =", h_tuple)

