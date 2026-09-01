#include "intcom.hpp"

// Initiate the parameters for class group
IntCom::IntCom() {
    
    // Hardcoded parameters for testing (run class_group.py to reproduce these values)
    static const char* DELTA_STR = "-23923022441763514261374255739786997872154045257289336511219878864586757005883847742286508461631820307302220187308685760906140702068851306644750434419403587277865557539680380941053142852567384143183014267815968038334067863306588533292757980804343431174589385259282220240810559571427552836110886157467881094886440802877753492881430021649955280632003021008024078620405216665734908127940106881093763222655030182552935493449216919743220992862329961978035087007194361418148522199034072339";
    
    static const char* G_A_STR = "3";
    static const char* G_B_STR = "1";
    static const char* G_C_STR = "1993585203480292855114521311648916489346170438107444709268323238715563083823653978523875705135985025608518348942390480075511725172404275553729202868283632273155463128306698411754428571047282011931917855651330669861172321942215711107729831733695285931215782104940185020067546630952296069675907179788990091240536733573146124406785835137496273386000251750668673218367101388811242343995008906757813601887919181879411291120768076645268416071860830164836257250599530118179043516586172695";

    static const char* H_A_STR = "5";
    static const char* H_B_STR = "1";
    static const char* H_C_STR = "1196151122088175713068712786989349893607702262864466825560993943229337850294192387114325423081591015365111009365434288045307035103442565332237521720970179363893277876984019047052657142628369207159150713390798401916703393165329426664637899040217171558729469262964111012040527978571377641805544307873394054744322040143887674644071501082497764031600151050401203931020260833286745406397005344054688161132751509127646774672460845987161049643116498098901754350359718070907426109951703617";

    U_ = 1600/2 + 20; // Hard code upper bound of the group order since Delta is 1600 bits.

    pari_sp av = avma;
    
    delta_ = gclone(strtoi(const_cast<char*>(DELTA_STR)));
    g_ = gclone(makeQfb(G_A_STR, G_B_STR, G_C_STR));
    h_ = gclone(makeQfb(H_A_STR, H_B_STR, H_C_STR));
    
    avma = av;
}

IntCom::~IntCom() {
    for (auto& window : g_window_) {
        for (GEN g : window) gunclone(g);
    }
    for (auto& window : h_window_) {
        for (GEN h : window) gunclone(h);
    }
    gunclone(delta_);
    gunclone(g_);
    gunclone(h_);
}


IntCom::Com IntCom::commit(const NTL::ZZ& x, const long r_range) const {
    if (r_range < 0) {
        throw std::invalid_argument("r_range must be nonnegative");
    }
    //Generate random r in [0, 2^r_range)
    NTL::ZZ r = NTL::RandomBits_ZZ(r_range);

    pari_sp av = avma;
    GEN gx = window_pow(g_window_, g_, x);    // g^x
    GEN hr = window_pow(h_window_, h_, r);    // h^r
    GEN c_x = gmul(gx, hr); // c_x = g^x * h^r

    Com com;
    com.c_x = serializeForm(c_x);
    com.r = r;

    avma = av;

    return com;

}


bool IntCom::open(const std::string& c_x, const NTL::ZZ& x, const NTL::ZZ& r) const {
        pari_sp av = avma;

        GEN c_claimed = deserializeForm(c_x);
        GEN gx = window_pow(g_window_, g_, x);    // g^x
        GEN hr = window_pow(h_window_, h_, r);    // h^r
        GEN c_recomputed = gmul(gx, hr); // c_recomputed = g^x * h^r

        bool ok = gequal(c_claimed, c_recomputed);

        avma = av;
        return ok;
}

long IntCom::commitment_bitlength(const std::string& c_x) const {
    pari_sp av = avma;

    GEN c_gen = deserializeForm(c_x);

    long bits = expi(gel(c_gen, 1)) + expi(gel(c_gen, 2)) + expi(gel(c_gen, 3)) + 1; // +1 for the sign bit

    avma = av;
    return bits;
}

std::string IntCom::commit_elem(const NTL::ZZ& x, const NTL::ZZ& y) const {
    pari_sp av = avma;

    GEN gx = window_pow(g_window_, g_, x);
    GEN hy = window_pow(h_window_, h_, y);
    GEN c_x = gmul(gx, hy);

    std::string s = serializeForm(c_x);
    avma = av;
    return s;
}

std::string IntCom::g_str() const {
    pari_sp av = avma;
    std::string s = serializeForm(g_);
    avma = av;
    return s;
}

std::string IntCom::h_str() const {
    pari_sp av = avma;
    std::string s = serializeForm(h_);
    avma = av;
    return s;
}

std::string IntCom::pow(const std::string& a, const NTL::ZZ& e) const {
    pari_sp av = avma;
    GEN a_gen = deserializeForm(a);
    GEN e_gen = zzToGEN(e);
    GEN out = gpow(a_gen, e_gen, 0);

    std::string s = serializeForm(out);
    avma = av;
    return s;
}

std::string IntCom::mul(const std::string& a, const std::string& b) const {
    pari_sp av = avma;
    GEN a_gen = deserializeForm(a);
    GEN b_gen = deserializeForm(b);
    GEN out = gmul(a_gen, b_gen);

    std::string s = serializeForm(out);
    avma = av;
    return s;
}

std::string IntCom::multi_pow_mul(const std::vector<std::string>& bases, const std::vector<NTL::ZZ>& exps) const {
    if (bases.empty() || bases.size() != exps.size()) {
        throw std::invalid_argument("multi_pow_mul: bases and exps must be non-empty and of equal length");
    }

    pari_sp av = avma;

    const std::size_t k = bases.size();
    std::vector<NTL::ZZ> mag(k); // |exps[i]|
    long max_bits = 0;

    // Per-base window table: table[i][j-1] = operand_i^j for j = 1..2^W-1, where operand_i is
    // bases[i] (or its inverse if exps[i] < 0). Built fresh each call -- unlike g_pow/h_pow's
    // table, which is amortized across every call over this object's lifetime since g and h
    // never change, these bases differ on every call, so there's no reuse to amortize into.
    std::vector<std::vector<GEN>> table(k);
    for (std::size_t i = 0; i < k; ++i) {
        GEN b = deserializeForm(bases[i]);
        mag[i] = NTL::abs(exps[i]);
        GEN operand = (exps[i] < 0) ? ginv(b) : b;
        max_bits = std::max(max_bits, NTL::NumBits(mag[i]));

        table[i].reserve(MULTI_WINDOW_DIGITS);
        table[i].push_back(operand); // j = 1
        for (int j = 2; j <= MULTI_WINDOW_DIGITS; ++j) {
            table[i].push_back(gmul(table[i].back(), operand));
        }
    }

    GEN acc = gpow(table[0][0], gen_0, 0); // identity element of the (shared) group
    long num_windows = (max_bits + MULTI_WINDOW_BITS - 1) / MULTI_WINDOW_BITS;

    // Reclaim each iteration's gsqr/gmul scratch immediately instead of letting it accumulate
    // for the whole loop: exponents can span thousands of bits (e.g. CRT coefficients scale
    // with the bit length of the product of all parties' moduli), and without this the
    // unreclaimed garbage can overflow the PARI stack long before the final result does.
    for (long w = num_windows - 1; w >= 0; --w) {
        pari_sp av_loop = avma;
        for (int s = 0; s < MULTI_WINDOW_BITS; ++s) {
            acc = gsqr(acc); // advance the shared accumulator by one window's worth of bits
        }
        for (std::size_t i = 0; i < k; ++i) {
            int digit = 0;
            for (int b = 0; b < MULTI_WINDOW_BITS; ++b) {
                if (NTL::bit(mag[i], MULTI_WINDOW_BITS * w + b)) {
                    digit |= (1 << b);
                }
            }
            if (digit != 0) {
                acc = gmul(acc, table[i][digit - 1]);
            }
        }
        acc = gerepileupto(av_loop, acc);
    }

    std::string s = serializeForm(acc);
    avma = av;
    return s;
}

std::string IntCom::g_pow(const NTL::ZZ& x) const {
    pari_sp av = avma;
    std::string s = serializeForm(window_pow(g_window_, g_, x));
    avma = av;
    return s;
}

std::string IntCom::h_pow(const NTL::ZZ& y) const {
    pari_sp av = avma;
    std::string s = serializeForm(window_pow(h_window_, h_, y));
    avma = av;
    return s;
}

void IntCom::warmup(long max_bits) const {
    long num_windows = (max_bits + WINDOW_BITS - 1) / WINDOW_BITS;
    ensure_window_table(g_window_, g_, num_windows);
    ensure_window_table(h_window_, h_, num_windows);
}

// extend table so that table[k] (k = 0 .. num_windows_needed-1) holds
// [base^(1*2^(W*k)), base^(2*2^(W*k)), ..., base^((2^W-1)*2^(W*k))].
// Built one window at a time: window k's j=1 entry ("boundary") is either base itself (k=0)
// or the previous window's boundary squared W times; entries j=2..2^W-1 are then a simple
// running product by the boundary. Every new entry is cloned onto PARI's permanent heap so
// it survives avma resets elsewhere; avma is reset after each window so volatile stack usage
// stays bounded regardless of how far the table is extended.
void IntCom::ensure_window_table(std::vector<std::vector<GEN>>& table, GEN base, long num_windows_needed) {
    while (static_cast<long>(table.size()) < num_windows_needed) {
        pari_sp av = avma;

        GEN boundary; // base^(2^(W*k)) for the window k being built
        if (table.empty()) {
            boundary = base;
        } else {
            GEN cur = table.back()[0];
            for (int i = 0; i < WINDOW_BITS; ++i) {
                cur = gsqr(cur);
            }
            boundary = cur;
        }

        std::vector<GEN> window;
        window.reserve(WINDOW_DIGITS);
        GEN cur = boundary;
        window.push_back(gclone(cur)); // j = 1
        for (int j = 2; j <= WINDOW_DIGITS; ++j) {
            cur = gmul(cur, boundary);
            window.push_back(gclone(cur)); // j
        }

        avma = av;
        table.push_back(std::move(window));
    }
}

// compute base^e using the (lazily extended) windowed table: split |e| into WINDOW_BITS-sized
// windows, multiply in the precomputed value for each window's nonzero digit, then invert if
// e is negative. No squaring is ever done here directly — the squaring chain and per-window
// digit precomputation live entirely in ensure_window_table and are paid at most once per
// window position across the lifetime of this IntCom instance.
GEN IntCom::window_pow(std::vector<std::vector<GEN>>& table, GEN base, const NTL::ZZ& e) {
    NTL::ZZ mag = NTL::abs(e);
    long bits = NTL::NumBits(mag);

    if (bits == 0) {
        return gpow(base, gen_0, 0); // e == 0: identity element
    }

    long num_windows = (bits + WINDOW_BITS - 1) / WINDOW_BITS;
    ensure_window_table(table, base, num_windows);

    // Combine the (already-precomputed) per-window digit values. Reclaim each iteration's gmul
    // scratch immediately: for large exponents this can run thousands of times, and without
    // reclaiming, that unreclaimed garbage can overflow the PARI stack on its own.
    GEN acc = nullptr;
    for (long k = 0; k < num_windows; ++k) {
        int digit = 0;
        for (int b = 0; b < WINDOW_BITS; ++b) {
            if (NTL::bit(mag, WINDOW_BITS * k + b)) {
                digit |= (1 << b);
            }
        }
        if (digit != 0) {
            GEN term = table[k][digit - 1]; // gclone'd: lives on the permanent heap, not the regular stack
            if (!acc) {
                acc = term; // no stack garbage created yet, nothing to reclaim
            } else {
                pari_sp av_loop = avma;
                acc = gerepileupto(av_loop, gmul(acc, term));
            }
        }
    }

    if (!acc) {
        acc = gpow(base, gen_0, 0); // shouldn't happen since bits > 0, but guard anyway
    }

    return (e < 0) ? ginv(acc) : acc;
}

GEN IntCom::zzToGEN(const NTL::ZZ& x) {
    // strtoi() only parses unsigned decimal digits (it returns 0 on a leading '-'),
    // so the sign has to be handled separately.
    std::ostringstream oss;
    oss << NTL::abs(x);
    std::string s = oss.str();
    GEN g = strtoi(const_cast<char*>(s.c_str()));
    return (x < 0) ? gneg(g) : g;
}

std::string IntCom::genToString(GEN x) {
    char* raw = GENtostr(x);
    std::string out(raw);
    pari_free(raw);
    return out;
}

std::string IntCom::serializeForm(GEN f) {
    return genToString(gel(f, 1)) + "|" + genToString(gel(f, 2)) + "|" + genToString(gel(f, 3));
}

GEN IntCom::deserializeForm(const std::string& s) {
    size_t p1 = s.find('|');
    if (p1 == std::string::npos) {
        throw std::invalid_argument("bad serialized form");
    }
    
    size_t p2 = s.find('|', p1 + 1);
    if (p2 == std::string::npos) {
        throw std::invalid_argument("bad serialized form");
    }

    std::string a = s.substr(0, p1);
    std::string b = s.substr(p1 + 1, p2 - p1 - 1);
    std::string c = s.substr(p2 + 1);
    return makeQfb(a, b, c);
}

GEN IntCom::makeQfb(const std::string& a, const std::string& b, const std::string& c) {
    std::string expr = "Qfb(" + a + "," + b + "," + c + ")";
    return gp_read_str(const_cast<char*>(expr.c_str()));
}