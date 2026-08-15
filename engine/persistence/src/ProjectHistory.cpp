#include "aether/persistence/ProjectHistory.hpp"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace aether::persistence {

namespace {

constexpr char kMagic[4] = {'A', 'E', 'H', 'I'};
constexpr std::uint32_t kVersion = 1;

void writeRaw(std::ostream& out, const void* data, std::size_t bytes) {
    out.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
}

void readRaw(std::istream& in, void* data, std::size_t bytes) {
    in.read(static_cast<char*>(data), static_cast<std::streamsize>(bytes));
    if (!in) {
        throw std::runtime_error("ProjectHistory: unexpected end of index file (truncated)");
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

double currentTimestamp() {
    return std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
}

} // namespace

ProjectHistory::ProjectHistory(std::string directory) : directory_(std::move(directory)) {
    std::filesystem::create_directories(directory_);
}

std::string ProjectHistory::indexPath() const {
    return (std::filesystem::path(directory_) / "index.aehi").string();
}

std::vector<HistoryEntry> ProjectHistory::entries() const {
    const std::string path = indexPath();
    if (!std::filesystem::exists(path)) {
        return {};
    }

    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("ProjectHistory::entries: could not open '" + path + "' for reading");
    }

    char magic[4];
    readRaw(in, magic, sizeof(magic));
    if (magic[0] != kMagic[0] || magic[1] != kMagic[1] || magic[2] != kMagic[2] || magic[3] != kMagic[3]) {
        throw std::runtime_error("ProjectHistory::entries: '" + path + "' is not an Aether history index (bad magic)");
    }
    std::uint32_t version = 0;
    readRaw(in, &version, sizeof(version));
    if (version != kVersion) {
        throw std::runtime_error("ProjectHistory::entries: '" + path + "' has unsupported version " +
                                  std::to_string(version));
    }

    std::uint64_t count = 0;
    readRaw(in, &count, sizeof(count));

    std::vector<HistoryEntry> result;
    result.reserve(count);
    for (std::uint64_t i = 0; i < count; ++i) {
        HistoryEntry entry;
        entry.label = readString(in);
        entry.filename = readString(in);
        readRaw(in, &entry.timestamp, sizeof(entry.timestamp));
        result.push_back(std::move(entry));
    }
    return result;
}

void ProjectHistory::record(const std::string& label, const FieldArchive& archive) {
    std::vector<HistoryEntry> existing = entries();

    std::ostringstream filenameStream;
    filenameStream << "checkpoint_" << std::setw(4) << std::setfill('0') << existing.size() << ".aecf";

    HistoryEntry newEntry;
    newEntry.label = label;
    newEntry.filename = filenameStream.str();
    newEntry.timestamp = currentTimestamp();

    archive.save((std::filesystem::path(directory_) / newEntry.filename).string());

    existing.push_back(std::move(newEntry));

    const std::string path = indexPath();
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("ProjectHistory::record: could not open '" + path + "' for writing");
    }
    writeRaw(out, kMagic, sizeof(kMagic));
    writeRaw(out, &kVersion, sizeof(kVersion));
    const std::uint64_t count = existing.size();
    writeRaw(out, &count, sizeof(count));
    for (const HistoryEntry& entry : existing) {
        writeString(out, entry.label);
        writeString(out, entry.filename);
        writeRaw(out, &entry.timestamp, sizeof(entry.timestamp));
    }
    if (!out) {
        throw std::runtime_error("ProjectHistory::record: write failed for '" + path + "'");
    }
}

FieldArchive ProjectHistory::load(const HistoryEntry& entry) const {
    return FieldArchive::load((std::filesystem::path(directory_) / entry.filename).string());
}

} // namespace aether::persistence
