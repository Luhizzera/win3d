#pragma once

#include <map>
#include <string>
#include <vector>

namespace aether::persistence {

// Module 11 (persistence), first pass: a generic checkpoint container --
// a named bag of scalar metadata (doubles) plus named double arrays -- that
// can be saved to and loaded back from a single binary file.
//
// Deliberately solver-agnostic: this class knows nothing about staggered
// grids, turbulence closures, or any other solver-specific detail. It is
// up to the caller (a solver's own code, or the Python orchestration layer
// above it) to decide which of a solver's fields (u, v, w, p, k, omega,
// nut, ...) and which scalar parameters (nx, ny, nz, dx, dy, dz, time,
// viscosity, ...) to store, and under what names. That keeps this module
// reusable across every solver in engine/solver without needing to grow a
// case for each one.
//
// No external dependency (no SQLite, no serialization library): this
// project's persistence needs at this scale -- a single user, no
// concurrent access, no relational queries -- don't need a database
// engine, so a hand-rolled format is the dependency-free choice, the same
// default this project has used everywhere except Module 10's CUDA (which
// had no dependency-free equivalent for raw GPU throughput).
//
// File format (little-endian; this project only targets Windows x64, so no
// endianness handling is implemented -- would need one if that changed):
//   char[4]    magic       = "AECF" (Aether Checkpoint Format)
//   uint32_t   version     = 1
//   uint64_t   metadataCount
//   metadataCount times:
//     uint64_t   keyLength
//     char[keyLength] key      (UTF-8, not null-terminated)
//     double     value
//   uint64_t   fieldCount
//   fieldCount times:
//     uint64_t   nameLength
//     char[nameLength] name    (UTF-8, not null-terminated)
//     uint64_t   elementCount
//     double[elementCount] data
class FieldArchive {
public:
    void setMetadata(const std::string& key, double value);
    double metadata(const std::string& key) const;
    bool hasMetadata(const std::string& key) const;

    void setField(const std::string& name, std::vector<double> data);
    const std::vector<double>& field(const std::string& name) const;
    bool hasField(const std::string& name) const;

    // Throws std::runtime_error if the file cannot be opened for writing.
    void save(const std::string& path) const;

    // Throws std::runtime_error if the file cannot be opened, is truncated,
    // or does not start with the expected magic/version.
    static FieldArchive load(const std::string& path);

private:
    std::map<std::string, double> metadata_;
    std::map<std::string, std::vector<double>> fields_;
};

} // namespace aether::persistence
