#include "hash.hpp"

#include <openssl/sha.h>
#include <NTL/ZZ.h>

#include <stdexcept>

namespace hash {

using namespace NTL;


void appendRaw(std::vector<unsigned char>& out,
               const unsigned char* data,
               size_t len) {
    out.insert(out.end(), data, data + len);
}

// 8-byte big-endian length prefix
void appendUint64(std::vector<unsigned char>& out, uint64_t x) {
    for (int i = 7; i >= 0; --i) {
        out.push_back(static_cast<unsigned char>((x >> (8 * i)) & 0xff));
    }
}

void appendStringEncoded(std::vector<unsigned char>& out, const std::string& s) {
    out.push_back(0x01); // type tag: string
    appendUint64(out, static_cast<uint64_t>(s.size()));
    appendRaw(out,
              reinterpret_cast<const unsigned char*>(s.data()),
              s.size());
}

void appendZZEncoded(std::vector<unsigned char>& out, const NTL::ZZ& z) {
    out.push_back(0x02); // type tag: ZZ

    long n = NTL::NumBytes(z);
    if (n < 0) {
        throw std::runtime_error("Invalid ZZ size");
    }

    appendUint64(out, static_cast<uint64_t>(n));

    if (n > 0) {
        size_t oldSize = out.size();
        out.resize(oldSize + static_cast<size_t>(n));
        NTL::BytesFromZZ(out.data() + oldSize, z, n);
    }
}

void appendHashEncoded(std::vector<unsigned char>& out, const Hash& h) {
    out.push_back(0x03); // type tag: Hash
    appendUint64(out, static_cast<uint64_t>(h.size()));
    appendRaw(out, h.data(), h.size());
}


static Hash sha256(const std::vector<unsigned char>& data) {
    Hash out{};
    SHA256(data.data(), data.size(), out.data());
    return out;
}

static std::vector<unsigned char> concatBytes(const Hash& a, const Hash& b) {
    std::vector<unsigned char> out;
    out.reserve(a.size() + b.size());
    out.insert(out.end(), a.begin(), a.end());
    out.insert(out.end(), b.begin(), b.end());
    return out;
}



Hash hashBytes(const unsigned char* data, size_t len) {
    Hash out{};
    SHA256(data, len, out.data());
    return out;
}

Hash hashBytes(const std::vector<unsigned char>& data) {
    return hashBytes(data.data(), data.size());
}

// Domain-separated hashing for ZZ
Hash hash(const NTL::ZZ& x) {
    long n = NTL::NumBytes(x);
    if (n < 0) {
        throw std::runtime_error("Invalid ZZ size");
    }

    std::vector<unsigned char> bytes;
    bytes.reserve(1 + static_cast<size_t>(n));
    bytes.push_back(0x00); // domain tag for ZZ

    if (n > 0) {
        size_t oldSize = bytes.size();
        bytes.resize(oldSize + static_cast<size_t>(n));
        NTL::BytesFromZZ(bytes.data() + oldSize, x, n);
    }

    return hashBytes(bytes);
}

// Domain-separated hashing for string
Hash hash(const std::string& s) {
    std::vector<unsigned char> bytes;
    bytes.reserve(1 + s.size());
    bytes.push_back(0x01); // domain tag for string
    appendRaw(bytes,
                reinterpret_cast<const unsigned char*>(s.data()),
                s.size());
    return hashBytes(bytes);
}

Hash hash(const char* s) {
    return hash(std::string(s));
}

// Domain-separated hashing for Hash itself
Hash hash(const Hash& h) {
    std::vector<unsigned char> bytes;
    bytes.reserve(1 + h.size());
    bytes.push_back(0x02); // domain tag for Hash
    appendRaw(bytes, h.data(), h.size());
    return hashBytes(bytes);
}

// Combine two hashes as an internal node
Hash hashCombine(const Hash& left, const Hash& right) {
    std::vector<unsigned char> bytes;
    bytes.reserve(1 + left.size() + right.size());
    bytes.push_back(0x03); // domain tag for internal/combine
    appendRaw(bytes, left.data(), left.size());
    appendRaw(bytes, right.data(), right.size());
    return hashBytes(bytes);
}

Hash hashConcat(const std::string& s, const NTL::ZZ& z, const Hash& h) {
    std::vector<unsigned char> buf;
    buf.reserve(1 + 8 + s.size() + 1 + 8 + static_cast<size_t>(NTL::NumBytes(z)) + 1 + 8 + h.size());

    appendStringEncoded(buf, s);
    appendZZEncoded(buf, z);
    appendHashEncoded(buf, h);

    return hashBytes(buf);
}

Hash hashConcat(const std::string& s, const std::string& t, const NTL::ZZ& z) {
    std::vector<unsigned char> buf;
    buf.reserve(1 + 8 + s.size() + 1 + 8 + t.size() + 1 + 8 + static_cast<size_t>(NTL::NumBytes(z)));

    appendStringEncoded(buf, s);
    appendStringEncoded(buf, t);
    appendZZEncoded(buf, z);

    return hashBytes(buf);
}

NTL::ZZ hashToZZ(const Hash& h) {
    NTL::ZZ z;
    NTL::ZZFromBytes(z, h.data(), h.size());
    return z;
}

NTL::ZZ hashToZZ128(const Hash& h) {
    NTL::ZZ z;
    NTL::ZZFromBytes(z, h.data(), 16);
    return z;
}

// ------------------------------
// MerkleTree implementation
// ------------------------------

MerkleTree::MerkleTree(const std::vector<Hash>& leaves) {
    if (leaves.empty()) {
        throw std::invalid_argument("Merkle tree requires at least one leaf");
    }

    std::vector<Hash> current = leaves;
    //current.reserve(leaves.size());

    //for (const auto& x : leaves) {
        //current.push_back(hash(x));
    //}

    levels_.push_back(current);

    while (current.size() > 1) {
        std::vector<Hash> next;
        next.reserve((current.size() + 1) / 2);

        for (size_t i = 0; i < current.size(); i += 2) {
            const Hash& left = current[i];
            const Hash& right =
                (i + 1 < current.size()) ? current[i + 1] : current[i];

            next.push_back(hashCombine(left, right));
        }

        levels_.push_back(next);
        current = next;
    }
}

Hash MerkleTree::root() const {
    return levels_.back().front();
}

AuthPath MerkleTree::authenticationPath(size_t leafIndex) const {
    if (leafIndex >= levels_.front().size()) {
        throw std::out_of_range("Leaf index out of range");
    }

    AuthPath path;
    size_t idx = leafIndex;

    for (size_t level = 0; level + 1 < levels_.size(); ++level) {
        const auto& nodes = levels_[level];

        size_t siblingIdx;
        bool siblingIsLeft;

        if (idx % 2 == 0) {
            siblingIdx = (idx + 1 < nodes.size()) ? idx + 1 : idx;
            siblingIsLeft = false;
        } else {
            siblingIdx = idx - 1;
            siblingIsLeft = true;
        }

        path.push_back({nodes[siblingIdx], siblingIsLeft});
        idx /= 2;
    }

    return path;
}

Hash MerkleTree::verify(const ZZ& leaf,
                        const AuthPath& path) {
    Hash current = hash(leaf);

    for (const auto& step : path) {
        if (step.siblingIsLeft) {
            current = hashCombine(step.sibling, current);
        } else {
            current = hashCombine(current, step.sibling);
        }
    }

    return current;
}

}
