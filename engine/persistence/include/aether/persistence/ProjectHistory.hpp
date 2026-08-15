#pragma once

#include "aether/persistence/FieldArchive.hpp"

#include <string>
#include <vector>

namespace aether::persistence {

struct HistoryEntry {
    std::string label;
    std::string filename; // relative to the ProjectHistory's own directory
    double timestamp;      // seconds since the Unix epoch (std::chrono::system_clock)
};

// Module 11, second piece: an ordered, disk-persisted list of checkpoints
// within one directory -- the "history" and "projects" half of the
// module's own name (FieldArchive alone only covers "results": a single
// checkpoint at a time, with no notion of a run's history).
//
// Deliberately thin: a ProjectHistory does not hold any checkpoint data in
// memory itself. record() writes the given FieldArchive to a generated
// filename inside its directory and appends one entry to an index file
// there; entries() re-reads that index file from disk every call rather
// than caching it, so two ProjectHistory instances pointed at the same
// directory (e.g. across separate process runs) always agree -- the index
// file on disk is the single source of truth, not any in-memory state.
//
// Index file format (little-endian, no endianness handling -- same
// rationale as FieldArchive: this project only targets Windows x64):
//   char[4]    magic   = "AEHI" (Aether History Index)
//   uint32_t   version = 1
//   uint64_t   entryCount
//   entryCount times:
//     uint64_t   labelLength
//     char[labelLength] label
//     uint64_t   filenameLength
//     char[filenameLength] filename
//     double     timestamp
class ProjectHistory {
public:
    // Creates `directory` if it does not already exist. Does not touch or
    // require an existing index file there -- the first record() call
    // creates one.
    explicit ProjectHistory(std::string directory);

    // Saves `archive` under a generated filename inside this history's
    // directory (checkpoint_0000.aecf, checkpoint_0001.aecf, ...) and
    // appends a new entry (label, that filename, the current wall-clock
    // time) to the index file. Throws std::runtime_error on any I/O
    // failure (mirrors FieldArchive::save()/load()).
    void record(const std::string& label, const FieldArchive& archive);

    // Reads the index file fresh from disk. Returns an empty vector if no
    // index file exists yet (a history with no recorded checkpoints is not
    // an error).
    std::vector<HistoryEntry> entries() const;

    // Convenience: loads the FieldArchive a given entry points to.
    FieldArchive load(const HistoryEntry& entry) const;

private:
    std::string indexPath() const;

    std::string directory_;
};

} // namespace aether::persistence
