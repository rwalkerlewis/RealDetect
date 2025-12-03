#pragma once

/**
 * CSS3.0 Database Interface
 * 
 * Provides SQLite-based storage following the CSS3.0 schema for
 * seismic event catalogs, arrivals, magnitudes, and station metadata.
 */

#include "css30_schema.hpp"
#include "../core/types.hpp"
#include "../core/event.hpp"
#include "../core/station.hpp"
#include <memory>
#include <vector>
#include <functional>
#include <mutex>

// Forward declare sqlite3 types
struct sqlite3;
struct sqlite3_stmt;

namespace seisproc {

/**
 * CSS30Database - SQLite database implementing CSS3.0 schema
 */
class CSS30Database {
public:
    CSS30Database();
    ~CSS30Database();
    
    // Prevent copying
    CSS30Database(const CSS30Database&) = delete;
    CSS30Database& operator=(const CSS30Database&) = delete;
    
    // Connection management
    bool open(const std::string& filename);
    bool isOpen() const { return db_ != nullptr; }
    void close();
    
    // Schema management
    bool createSchema();
    bool dropSchema();
    std::string schemaVersion() const;
    
    // Configuration
    void setAuthor(const std::string& auth) { author_ = auth; }
    const std::string& author() const { return author_; }
    
    // Event storage
    int64_t storeEvent(const Event& event);
    int64_t storeOrigin(int64_t evid, const Origin& origin);
    int64_t storeOriginError(int64_t orid, const Origin& origin);
    int64_t storeArrival(const Pick& pick, const StreamID& stream);
    void storeAssociation(int64_t arid, int64_t orid, const Arrival& arrival,
                          const std::string& vmodel = "-");
    int64_t storeMagnitude(int64_t orid, int64_t evid, const Magnitude& mag,
                           const std::string& net = "-");
    void storeStationMagnitude(int64_t magid, int64_t orid, int64_t evid,
                               const StationMagnitude& stamag);
    
    // Complete event storage (convenience method)
    bool storeCompleteEvent(const Event& event, 
                            const std::string& velocity_model = "-",
                            const std::string& network = "-");
    
    // Station/channel metadata storage
    void storeSite(const Station& station, double ondate = 2451545.0);
    void storeSitechan(const std::string& sta, const Channel& chan,
                       double ondate = 2451545.0);
    void storeAffiliation(const std::string& net, const std::string& sta);
    void storeInventory(const StationInventory& inventory, 
                        const std::string& network);
    
    // Waveform descriptor storage
    int64_t storeWfdisc(const std::string& sta, const std::string& chan,
                        double time, double endtime, int64_t nsamp,
                        double samprate, const std::string& dir,
                        const std::string& dfile);
    
    // Query methods
    std::vector<css30::Event> queryEvents(double starttime, double endtime,
                                          double minlat = -90, double maxlat = 90,
                                          double minlon = -180, double maxlon = 180,
                                          double minmag = -10, double maxmag = 10);
    
    std::vector<css30::Origin> queryOrigins(int64_t evid);
    std::vector<css30::Arrival> queryArrivals(int64_t orid);
    std::vector<css30::Assoc> queryAssociations(int64_t orid);
    std::vector<css30::Netmag> queryMagnitudes(int64_t orid);
    
    css30::Origin getPreferredOrigin(int64_t evid);
    
    // Statistics
    int64_t countEvents() const;
    int64_t countOrigins() const;
    int64_t countArrivals() const;
    std::pair<double, double> timeRange() const;
    
    // Export to CSS3.0 flat files
    bool exportToFlatFiles(const std::string& directory,
                           double starttime = 0, double endtime = 0);
    
    // Transaction support
    void beginTransaction();
    void commitTransaction();
    void rollbackTransaction();
    
    // Last error message
    std::string lastError() const { return last_error_; }

private:
    sqlite3* db_;
    std::string author_;
    std::string last_error_;
    mutable std::mutex mutex_;
    
    // ID generators
    int64_t next_evid_;
    int64_t next_orid_;
    int64_t next_arid_;
    int64_t next_magid_;
    int64_t next_wfid_;
    int64_t next_chanid_;
    int64_t next_commid_;
    
    // Prepared statements
    sqlite3_stmt* stmt_insert_event_;
    sqlite3_stmt* stmt_insert_origin_;
    sqlite3_stmt* stmt_insert_origerr_;
    sqlite3_stmt* stmt_insert_arrival_;
    sqlite3_stmt* stmt_insert_assoc_;
    sqlite3_stmt* stmt_insert_netmag_;
    sqlite3_stmt* stmt_insert_stamag_;
    sqlite3_stmt* stmt_insert_site_;
    sqlite3_stmt* stmt_insert_sitechan_;
    sqlite3_stmt* stmt_insert_affiliation_;
    sqlite3_stmt* stmt_insert_wfdisc_;
    sqlite3_stmt* stmt_update_prefor_;
    
    // Internal helpers
    bool initializeSequences();
    bool prepareStatements();
    void finalizeStatements();
    double currentLddate() const;
    int64_t epochToJdate(double epoch) const;
    bool executeSQL(const char* sql);
    void setError(const std::string& context);
};

/**
 * CSS30Writer - High-level interface for event catalog output
 */
class CSS30Writer {
public:
    CSS30Writer(CSS30Database& db) : db_(db) {}
    
    // Write a complete event with all associations
    bool writeEvent(const Event& event,
                    const std::string& velocity_model = "-",
                    const std::string& network = "-") {
        return db_.storeCompleteEvent(event, velocity_model, network);
    }
    
    // Batch write multiple events
    bool writeEvents(const std::vector<Event>& events,
                     const std::string& velocity_model = "-",
                     const std::string& network = "-") {
        db_.beginTransaction();
        try {
            for (const auto& event : events) {
                if (!writeEvent(event, velocity_model, network)) {
                    db_.rollbackTransaction();
                    return false;
                }
            }
            db_.commitTransaction();
            return true;
        } catch (...) {
            db_.rollbackTransaction();
            return false;
        }
    }
    
    // Write station inventory
    bool writeInventory(const StationInventory& inventory,
                        const std::string& network) {
        db_.storeInventory(inventory, network);
        return true;
    }

private:
    CSS30Database& db_;
};

} // namespace seisproc
