#include "aether/persistence/FieldArchive.hpp"

#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace aether::persistence {

namespace {

constexpr char kMagic[4] = {'A', 'E', 'C', 'F'};
constexpr std::uint32_t kVersion = 1;

void writeRaw(std::ostream& out, const void* data, std::size_t bytes) {
    out.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
}

void readRaw(std::istream& in, void* data, std::size_t bytes) {
    in.read(static_cast<char*>(data), static_cast<std::streamsize>(bytes));
    if (!in) {
        throw std::runtime_error("FieldArchive::load: unexpected end of file (truncated archive)");
    }
}

void writeString(std::ostream& out, const std::string& s) {
    const std::uint64_t length = s.size();
    writeRaw(out, &length, sizeof(length));
    writeRaw(out, s.data(), length);
}

std::string readString(std::istream& in) {
    std::uint64_t length = 0;
    readRaw(in, &length, sizeof(length));
    std::string s(length, '\0');
    if (length > 0) {
        readRaw(in, s.data(), length);
    }
    return s;
}

} // namespace

void FieldArchive::setMetadata(const std::string& key, double value) { metadata_[key] = value; }

double FieldArchive::metadata(const std::string& key) const { return metadata_.at(key); }

bool FieldArchive::hasMetadata(const std::string& key) const { return metadata_.find(key) != metadata_.end(); }

void FieldArchive::setField(const std::string& name, std::vector<double> data) { fields_[name] = std::move(data); }

const std::vector<double>& FieldArchive::field(const std::string& name) const { return fields_.at(name); }

bool FieldArchive::hasField(const std::string& name) const { return fields_.find(name) != fields_.end(); }

void FieldArchive::save(const std::string& path) const {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("FieldArchive::save: could not open '" + path + "' for writing");
    }

    writeRaw(out, kMagic, sizeof(kMagic));
    writeRaw(out, &kVersion, sizeof(kVersion));

    const std::uint64_t metadataCount = metadata_.size();
    writeRaw(out, &metadataCount, sizeof(metadataCount));
    for (const auto& [key, value] : metadata_) {
        writeString(out, key);
        writeRaw(out, &value, sizeof(value));
    }

    const std::uint64_t fieldCount = fields_.size();
    writeRaw(out, &fieldCount, sizeof(fieldCount));
    for (const auto& [name, data] : fields_) {
        writeString(out, name);
        const std::uint64_t elementCount = data.size();
        writeRaw(out, &elementCount, sizeof(elementCount));
        if (elementCount > 0) {
            writeRaw(out, data.data(), elementCount * sizeof(double));
        }
    }

    if (!out) {
        throw std::runtime_error("FieldArchive::save: write failed for '" + path + "'");
    }
}

FieldArchive FieldArchive::load(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("FieldArchive::load: could not open '" + path + "' for reading");
    }

    char magic[4];
    readRaw(in, magic, sizeof(magic));
    if (magic[0] != kMagic[0] || magic[1] != kMagic[1] || magic[2] != kMagic[2] || magic[3] != kMagic[3]) {
        throw std::runtime_error("FieldArchive::load: '" + path + "' is not an Aether checkpoint file (bad magic)");
    }

    std::uint32_t version = 0;
    readRaw(in, &version, sizeof(version));
    if (version != kVersion) {
        throw std::runtime_error("FieldArchive::load: '" + path + "' has unsupported version " +
                                  std::to_string(version) + " (expected " + std::to_string(kVersion) + ")");
    }

    FieldArchive archive;

    std::uint64_t metadataCount = 0;
    readRaw(in, &metadataCount, sizeof(metadataCount));
    for (std::uint64_t i = 0; i < metadataCount; ++i) {
        const std::string key = readString(in);
        double value = 0.0;
        readRaw(in, &value, sizeof(value));
        archive.metadata_[key] = value;
    }

    std::uint64_t fieldCount = 0;
    readRaw(in, &fieldCount, sizeof(fieldCount));
    for (std::uint64_t i = 0; i < fieldCount; ++i) {
        const std::string name = readString(in);
        std::uint64_t elementCount = 0;
        readRaw(in, &elementCount, sizeof(elementCount));
        std::vector<double> data(elementCount);
        if (elementCount > 0) {
            readRaw(in, data.data(), elementCount * sizeof(double));
        }
        archive.fields_[name] = std::move(data);
    }

    return archive;
}

} // namespace aether::persistence
