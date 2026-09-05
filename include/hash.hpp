#pragma once
#include <NTL/ZZ.h>
#include <array>
#include <vector>

namespace hash {

using namespace NTL;

// SHA-256 output size
constexpr size_t HASH_SIZE = 32;

using Hash = std::array<unsigned char, HASH_SIZE>;

struct AuthStep {
    Hash sibling;
    bool siblingIsLeft;
};

using AuthPath = std::vector<AuthStep>;

// ------------------------------
// Hash functions
// ------------------------------

// Hash hashLeaf(const ZZ& x);
// Hash hashInternal(const Hash& left, const Hash& right);
Hash hashBytes(const unsigned char* data, size_t len);
Hash hashBytes(const std::vector<unsigned char>& data);

// Typed hashing
Hash hash(const NTL::ZZ& x);
Hash hash(const std::string& s);
Hash hash(const char* s);
Hash hash(const Hash& h);

Hash hashConcat(const std::string& s, const NTL::ZZ& z, const Hash& h);
Hash hashConcat(const std::string& s, const std::string& t, const NTL::ZZ& z, const NTL::ZZ& w, const unsigned char bytes[32]);
Hash hashConcat(const NTL::ZZ& z, const NTL::ZZ& w, const unsigned char bytes[32]);

Hash hashCombine(const Hash& left, const Hash& right);

// Hash to ZZ
NTL::ZZ hashToZZ(const Hash& h);
NTL::ZZ hashToZZ128(const Hash& h);

// ------------------------------
// Merkle Tree
// ------------------------------

class MerkleTree {
public:

    explicit MerkleTree(const std::vector<Hash>& leaves);

    Hash root() const;

    AuthPath authenticationPath(size_t leafIndex) const;

    static Hash verify(const Hash& leaf,
                       const AuthPath& path);

private:
    std::vector<std::vector<Hash>> levels_;
};

}

