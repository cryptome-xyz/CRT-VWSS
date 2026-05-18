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
    GEN x_gen = zzToGEN(x);
    GEN r_gen = zzToGEN(r);

    GEN gx = gpow(g_, x_gen, 0);    // g^x
    GEN hr = gpow(h_, r_gen, 0);    // h^r
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
        GEN x_gen = zzToGEN(x);
        GEN r_gen = zzToGEN(r);

        GEN gx = gpow(g_, x_gen, 0);    // g^x
        GEN hr = gpow(h_, r_gen, 0);    // h^r
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
    
    GEN x_gen = zzToGEN(x);
    GEN y_gen = zzToGEN(y);

    GEN gx = gpow(g_, x_gen, 0);
    GEN hy = gpow(h_, y_gen, 0);
    GEN c_x = gmul(gx, hy);

    avma = av;

    return serializeForm(c_x);

}

std::string IntCom::pow(const std::string& a, const NTL::ZZ& e) const {
    //if (e < 0) {
        //return pow(inv(a), -e);
    //}

    pari_sp av = avma;
    GEN a_gen = deserializeForm(a);
    GEN e_gen = zzToGEN(e);
    GEN out = gpow(a_gen, e_gen, 0);

    std::string s = serializeForm(out);
    avma = av;
    return s;
}

GEN IntCom::zzToGEN(const NTL::ZZ& x) {
    std::ostringstream oss;
    oss << x;
    std::string s = oss.str();
    return strtoi(const_cast<char*>(s.c_str()));
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