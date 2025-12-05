/**
 * CSS3.0 Database Implementation
 * 
 * SQLite-based implementation of the CSS3.0 seismic data schema.
 */

#include "realdetect/database/css30_database.hpp"
#include <sqlite3.h>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <chrono>
#include <cmath>

namespace realdetect {

CSS30Database::CSS30Database()
    : db_(nullptr)
    , author_("seisproc")
    , next_evid_(1)
    , next_orid_(1)
    , next_arid_(1)
    , next_magid_(1)
    , next_wfid_(1)
    , next_chanid_(1)
    , next_commid_(1)
    , stmt_insert_event_(nullptr)
    , stmt_insert_origin_(nullptr)
    , stmt_insert_origerr_(nullptr)
    , stmt_insert_arrival_(nullptr)
    , stmt_insert_assoc_(nullptr)
    , stmt_insert_netmag_(nullptr)
    , stmt_insert_stamag_(nullptr)
    , stmt_insert_site_(nullptr)
    , stmt_insert_sitechan_(nullptr)
    , stmt_insert_affiliation_(nullptr)
    , stmt_insert_wfdisc_(nullptr)
    , stmt_update_prefor_(nullptr)
{
}

CSS30Database::~CSS30Database() {
    close();
}

bool CSS30Database::open(const std::string& filename) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (db_) {
        close();
    }
    
    int rc = sqlite3_open(filename.c_str(), &db_);
    if (rc != SQLITE_OK) {
        setError("Failed to open database");
        db_ = nullptr;
        return false;
    }
    
    // Enable foreign keys
    executeSQL("PRAGMA foreign_keys = ON;");
    
    // Performance settings
    executeSQL("PRAGMA journal_mode = WAL;");
    executeSQL("PRAGMA synchronous = NORMAL;");
    executeSQL("PRAGMA cache_size = -64000;");  // 64MB cache
    
    // Initialize sequence counters
    if (!initializeSequences()) {
        close();
        return false;
    }
    
    // Prepare statements
    if (!prepareStatements()) {
        close();
        return false;
    }
    
    return true;
}

void CSS30Database::close() {
    if (db_) {
        finalizeStatements();
        sqlite3_close(db_);
        db_ = nullptr;
    }
}

bool CSS30Database::createSchema() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_) return false;
    
    // Create all CSS3.0 tables
    if (!executeSQL(css30::Event::CREATE_SQL)) return false;
    if (!executeSQL(css30::Origin::CREATE_SQL)) return false;
    if (!executeSQL(css30::Origerr::CREATE_SQL)) return false;
    if (!executeSQL(css30::Arrival::CREATE_SQL)) return false;
    if (!executeSQL(css30::Assoc::CREATE_SQL)) return false;
    if (!executeSQL(css30::Netmag::CREATE_SQL)) return false;
    if (!executeSQL(css30::Stamag::CREATE_SQL)) return false;
    if (!executeSQL(css30::Site::CREATE_SQL)) return false;
    if (!executeSQL(css30::Sitechan::CREATE_SQL)) return false;
    if (!executeSQL(css30::Affiliation::CREATE_SQL)) return false;
    if (!executeSQL(css30::Wfdisc::CREATE_SQL)) return false;
    if (!executeSQL(css30::Remark::CREATE_SQL)) return false;
    
    // Create indices
    if (!executeSQL(css30::CREATE_INDICES_SQL)) return false;
    
    // Create metadata table for schema version
    executeSQL(R"(
        CREATE TABLE IF NOT EXISTS css30_meta (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    )");
    
    executeSQL("INSERT OR REPLACE INTO css30_meta VALUES ('schema_version', '3.0')");
    executeSQL("INSERT OR REPLACE INTO css30_meta VALUES ('created', datetime('now'))");
    
    return true;
}

bool CSS30Database::dropSchema() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_) return false;
    
    const char* tables[] = {
        "assoc", "stamag", "netmag", "origerr", "arrival", 
        "origin", "event", "wfdisc", "sitechan", "site",
        "affiliation", "remark", "css30_meta"
    };
    
    for (const auto& table : tables) {
        std::string sql = "DROP TABLE IF EXISTS " + std::string(table);
        executeSQL(sql.c_str());
    }
    
    return true;
}

std::string CSS30Database::schemaVersion() const {
    if (!db_) return "";
    
    sqlite3_stmt* stmt;
    std::string version;
    
    if (sqlite3_prepare_v2(db_, 
            "SELECT value FROM css30_meta WHERE key = 'schema_version'",
            -1, &stmt, nullptr) == SQLITE_OK) {
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            version = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        }
        sqlite3_finalize(stmt);
    }
    
    return version;
}

int64_t CSS30Database::storeEvent(const Event& event) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_event_) return -1;
    
    int64_t evid = next_evid_++;
    double lddate = currentLddate();
    
    sqlite3_reset(stmt_insert_event_);
    sqlite3_bind_int64(stmt_insert_event_, 1, evid);
    sqlite3_bind_text(stmt_insert_event_, 2, event.id().c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt_insert_event_, 3, -1);  // prefor set later
    sqlite3_bind_text(stmt_insert_event_, 4, author_.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_event_, 5, lddate);
    
    if (sqlite3_step(stmt_insert_event_) != SQLITE_DONE) {
        setError("Failed to insert event");
        return -1;
    }
    
    return evid;
}

int64_t CSS30Database::storeOrigin(int64_t evid, const Origin& origin) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_origin_) return -1;
    
    int64_t orid = next_orid_++;
    double lddate = currentLddate();
    
    // Convert TimePoint to epoch seconds
    auto epoch = std::chrono::duration_cast<std::chrono::microseconds>(
        origin.time.time_since_epoch()).count() / 1e6;
    
    sqlite3_reset(stmt_insert_origin_);
    sqlite3_bind_double(stmt_insert_origin_, 1, origin.location.latitude);
    sqlite3_bind_double(stmt_insert_origin_, 2, origin.location.longitude);
    sqlite3_bind_double(stmt_insert_origin_, 3, origin.location.depth);
    sqlite3_bind_double(stmt_insert_origin_, 4, epoch);
    sqlite3_bind_int64(stmt_insert_origin_, 5, orid);
    sqlite3_bind_int64(stmt_insert_origin_, 6, evid);
    sqlite3_bind_int64(stmt_insert_origin_, 7, epochToJdate(epoch));
    sqlite3_bind_int64(stmt_insert_origin_, 8, origin.arrivals.size());  // nass
    sqlite3_bind_int64(stmt_insert_origin_, 9, origin.phase_count);      // ndef
    sqlite3_bind_text(stmt_insert_origin_, 10, "eq", -1, SQLITE_STATIC); // etype
    sqlite3_bind_text(stmt_insert_origin_, 11, origin.is_fixed_depth ? "g" : "f", 
                      -1, SQLITE_STATIC);  // dtype
    sqlite3_bind_text(stmt_insert_origin_, 12, origin.algorithm.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt_insert_origin_, 13, author_.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_origin_, 14, lddate);
    
    if (sqlite3_step(stmt_insert_origin_) != SQLITE_DONE) {
        setError("Failed to insert origin");
        return -1;
    }
    
    // Update event's preferred origin
    sqlite3_reset(stmt_update_prefor_);
    sqlite3_bind_int64(stmt_update_prefor_, 1, orid);
    sqlite3_bind_int64(stmt_update_prefor_, 2, evid);
    sqlite3_step(stmt_update_prefor_);
    
    return orid;
}

int64_t CSS30Database::storeOriginError(int64_t orid, const Origin& origin) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_origerr_) return -1;
    
    double lddate = currentLddate();
    
    // Convert error values to covariance matrix (simplified)
    double sxx = origin.longitude_error * origin.longitude_error;
    double syy = origin.latitude_error * origin.latitude_error;
    double szz = origin.depth_error * origin.depth_error;
    double stt = origin.time_error * origin.time_error;
    
    // Error ellipse (simplified - using lat/lon errors)
    double smajax = std::max(origin.latitude_error, origin.longitude_error);
    double sminax = std::min(origin.latitude_error, origin.longitude_error);
    double strike = origin.longitude_error > origin.latitude_error ? 90.0 : 0.0;
    
    sqlite3_reset(stmt_insert_origerr_);
    sqlite3_bind_int64(stmt_insert_origerr_, 1, orid);
    sqlite3_bind_double(stmt_insert_origerr_, 2, sxx);
    sqlite3_bind_double(stmt_insert_origerr_, 3, syy);
    sqlite3_bind_double(stmt_insert_origerr_, 4, szz);
    sqlite3_bind_double(stmt_insert_origerr_, 5, stt);
    sqlite3_bind_double(stmt_insert_origerr_, 6, origin.rms);  // sdobs
    sqlite3_bind_double(stmt_insert_origerr_, 7, smajax);
    sqlite3_bind_double(stmt_insert_origerr_, 8, sminax);
    sqlite3_bind_double(stmt_insert_origerr_, 9, strike);
    sqlite3_bind_double(stmt_insert_origerr_, 10, origin.depth_error);
    sqlite3_bind_double(stmt_insert_origerr_, 11, origin.time_error);
    sqlite3_bind_double(stmt_insert_origerr_, 12, 0.9);  // conf
    sqlite3_bind_double(stmt_insert_origerr_, 13, lddate);
    
    if (sqlite3_step(stmt_insert_origerr_) != SQLITE_DONE) {
        setError("Failed to insert origerr");
        return -1;
    }
    
    return orid;
}

int64_t CSS30Database::storeArrival(const Pick& pick, const StreamID& stream) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_arrival_) return -1;
    
    int64_t arid = next_arid_++;
    double lddate = currentLddate();
    
    // Convert TimePoint to epoch
    auto epoch = std::chrono::duration_cast<std::chrono::microseconds>(
        pick.time.time_since_epoch()).count() / 1e6;
    
    // Quality to CSS3.0 format
    const char* qual;
    switch (pick.quality) {
        case PickQuality::Impulsive: qual = "i"; break;
        case PickQuality::Emergent: qual = "e"; break;
        default: qual = "?"; break;
    }
    
    sqlite3_reset(stmt_insert_arrival_);
    sqlite3_bind_text(stmt_insert_arrival_, 1, stream.station.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_arrival_, 2, epoch);
    sqlite3_bind_int64(stmt_insert_arrival_, 3, arid);
    sqlite3_bind_int64(stmt_insert_arrival_, 4, epochToJdate(epoch));
    sqlite3_bind_text(stmt_insert_arrival_, 5, stream.channel.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt_insert_arrival_, 6, phaseTypeToString(pick.phase_type).c_str(), 
                      -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_arrival_, 7, pick.uncertainty);  // deltim
    sqlite3_bind_double(stmt_insert_arrival_, 8, pick.amplitude);
    sqlite3_bind_double(stmt_insert_arrival_, 9, pick.period);
    sqlite3_bind_double(stmt_insert_arrival_, 10, pick.snr);
    sqlite3_bind_text(stmt_insert_arrival_, 11, qual, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt_insert_arrival_, 12, author_.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_arrival_, 13, lddate);
    
    if (sqlite3_step(stmt_insert_arrival_) != SQLITE_DONE) {
        setError("Failed to insert arrival");
        return -1;
    }
    
    return arid;
}

void CSS30Database::storeAssociation(int64_t arid, int64_t orid, 
                                      const Arrival& arrival,
                                      const std::string& vmodel) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_assoc_) return;
    
    double lddate = currentLddate();
    
    // Convert distance km to degrees
    double delta = arrival.distance / 111.195;
    
    sqlite3_reset(stmt_insert_assoc_);
    sqlite3_bind_int64(stmt_insert_assoc_, 1, arid);
    sqlite3_bind_int64(stmt_insert_assoc_, 2, orid);
    sqlite3_bind_text(stmt_insert_assoc_, 3, 
                      arrival.pick->stream_id.station.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt_insert_assoc_, 4, 
                      phaseTypeToString(arrival.pick->phase_type).c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_assoc_, 5, delta);
    sqlite3_bind_double(stmt_insert_assoc_, 6, arrival.azimuth);  // seaz
    sqlite3_bind_double(stmt_insert_assoc_, 7, std::fmod(arrival.azimuth + 180.0, 360.0));  // esaz
    sqlite3_bind_double(stmt_insert_assoc_, 8, arrival.residual);
    sqlite3_bind_text(stmt_insert_assoc_, 9, arrival.used ? "d" : "n", -1, SQLITE_STATIC);
    sqlite3_bind_double(stmt_insert_assoc_, 10, arrival.weight);
    sqlite3_bind_text(stmt_insert_assoc_, 11, vmodel.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_assoc_, 12, lddate);
    
    sqlite3_step(stmt_insert_assoc_);
}

int64_t CSS30Database::storeMagnitude(int64_t orid, int64_t evid, 
                                       const Magnitude& mag,
                                       const std::string& net) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_netmag_) return -1;
    
    int64_t magid = next_magid_++;
    double lddate = currentLddate();
    
    sqlite3_reset(stmt_insert_netmag_);
    sqlite3_bind_int64(stmt_insert_netmag_, 1, magid);
    sqlite3_bind_text(stmt_insert_netmag_, 2, net.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt_insert_netmag_, 3, orid);
    sqlite3_bind_int64(stmt_insert_netmag_, 4, evid);
    sqlite3_bind_text(stmt_insert_netmag_, 5, 
                      magnitudeTypeToString(mag.type).c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt_insert_netmag_, 6, mag.station_count);
    sqlite3_bind_double(stmt_insert_netmag_, 7, mag.value);
    sqlite3_bind_double(stmt_insert_netmag_, 8, mag.uncertainty);
    sqlite3_bind_text(stmt_insert_netmag_, 9, author_.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_netmag_, 10, lddate);
    
    if (sqlite3_step(stmt_insert_netmag_) != SQLITE_DONE) {
        setError("Failed to insert netmag");
        return -1;
    }
    
    return magid;
}

void CSS30Database::storeStationMagnitude(int64_t magid, int64_t orid, int64_t evid,
                                           const StationMagnitude& stamag) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_stamag_) return;
    
    double lddate = currentLddate();
    double delta = stamag.distance / 111.195;  // km to degrees
    
    sqlite3_reset(stmt_insert_stamag_);
    sqlite3_bind_int64(stmt_insert_stamag_, 1, magid);
    sqlite3_bind_text(stmt_insert_stamag_, 2, 
                      stamag.stream_id.station.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt_insert_stamag_, 3, orid);
    sqlite3_bind_int64(stmt_insert_stamag_, 4, evid);
    sqlite3_bind_double(stmt_insert_stamag_, 5, delta);
    sqlite3_bind_text(stmt_insert_stamag_, 6, 
                      magnitudeTypeToString(stamag.type).c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_stamag_, 7, stamag.value);
    sqlite3_bind_text(stmt_insert_stamag_, 8, author_.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_stamag_, 9, lddate);
    
    sqlite3_step(stmt_insert_stamag_);
}

bool CSS30Database::storeCompleteEvent(const Event& event,
                                        const std::string& velocity_model,
                                        const std::string& network) {
    if (!db_) return false;
    
    beginTransaction();
    
    try {
        // Store event
        int64_t evid = storeEvent(event);
        if (evid < 0) {
            rollbackTransaction();
            return false;
        }
        
        // Store each origin
        for (const auto& origin : event.origins()) {
            int64_t orid = storeOrigin(evid, origin);
            if (orid < 0) {
                rollbackTransaction();
                return false;
            }
            
            // Store origin errors
            storeOriginError(orid, origin);
            
            // Store arrivals and associations
            for (const auto& arrival : origin.arrivals) {
                int64_t arid = storeArrival(*arrival.pick, arrival.pick->stream_id);
                if (arid > 0) {
                    storeAssociation(arid, orid, arrival, velocity_model);
                }
            }
        }
        
        // Store magnitudes
        int64_t last_orid = next_orid_ - 1;
        for (const auto& mag : event.magnitudes()) {
            int64_t magid = storeMagnitude(last_orid, evid, mag, network);
            
            // Store station magnitudes
            for (const auto& stamag : event.stationMagnitudes()) {
                if (stamag.type == mag.type) {
                    storeStationMagnitude(magid, last_orid, evid, stamag);
                }
            }
        }
        
        commitTransaction();
        return true;
        
    } catch (...) {
        rollbackTransaction();
        return false;
    }
}

void CSS30Database::storeSite(const Station& station, double ondate) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_site_) return;
    
    double lddate = currentLddate();
    
    sqlite3_reset(stmt_insert_site_);
    sqlite3_bind_text(stmt_insert_site_, 1, station.code().c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_site_, 2, ondate);
    sqlite3_bind_double(stmt_insert_site_, 3, -1);  // offdate
    sqlite3_bind_double(stmt_insert_site_, 4, station.latitude());
    sqlite3_bind_double(stmt_insert_site_, 5, station.longitude());
    sqlite3_bind_double(stmt_insert_site_, 6, station.elevation() / 1000.0);  // m to km
    sqlite3_bind_text(stmt_insert_site_, 7, station.code().c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt_insert_site_, 8, "ss", -1, SQLITE_STATIC);  // statype
    sqlite3_bind_double(stmt_insert_site_, 9, lddate);
    
    sqlite3_step(stmt_insert_site_);
}

void CSS30Database::storeSitechan(const std::string& sta, const Channel& chan,
                                   double ondate) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_sitechan_) return;
    
    int64_t chanid = next_chanid_++;
    double lddate = currentLddate();
    
    sqlite3_reset(stmt_insert_sitechan_);
    sqlite3_bind_text(stmt_insert_sitechan_, 1, sta.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt_insert_sitechan_, 2, chan.code.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_sitechan_, 3, ondate);
    sqlite3_bind_int64(stmt_insert_sitechan_, 4, chanid);
    sqlite3_bind_double(stmt_insert_sitechan_, 5, -1);  // offdate
    sqlite3_bind_double(stmt_insert_sitechan_, 6, chan.azimuth);  // hang
    sqlite3_bind_double(stmt_insert_sitechan_, 7, chan.dip + 90);  // vang (0=up, 90=horizontal)
    sqlite3_bind_text(stmt_insert_sitechan_, 8, chan.sensor_type.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_sitechan_, 9, lddate);
    
    sqlite3_step(stmt_insert_sitechan_);
}

void CSS30Database::storeAffiliation(const std::string& net, const std::string& sta) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_affiliation_) return;
    
    double lddate = currentLddate();
    
    sqlite3_reset(stmt_insert_affiliation_);
    sqlite3_bind_text(stmt_insert_affiliation_, 1, net.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt_insert_affiliation_, 2, sta.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_affiliation_, 3, lddate);
    
    sqlite3_step(stmt_insert_affiliation_);
}

void CSS30Database::storeInventory(const StationInventory& inventory,
                                    const std::string& network) {
    beginTransaction();
    
    for (const auto& [key, station] : inventory.stations()) {
        storeSite(*station);
        storeAffiliation(network, station->code());
        
        for (const auto& [code, chan] : station->channels()) {
            storeSitechan(station->code(), chan);
        }
    }
    
    commitTransaction();
}

int64_t CSS30Database::storeWfdisc(const std::string& sta, const std::string& chan,
                                    double time, double endtime, int64_t nsamp,
                                    double samprate, const std::string& dir,
                                    const std::string& dfile) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_ || !stmt_insert_wfdisc_) return -1;
    
    int64_t wfid = next_wfid_++;
    double lddate = currentLddate();
    
    sqlite3_reset(stmt_insert_wfdisc_);
    sqlite3_bind_text(stmt_insert_wfdisc_, 1, sta.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt_insert_wfdisc_, 2, chan.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_wfdisc_, 3, time);
    sqlite3_bind_int64(stmt_insert_wfdisc_, 4, wfid);
    sqlite3_bind_int64(stmt_insert_wfdisc_, 5, epochToJdate(time));
    sqlite3_bind_double(stmt_insert_wfdisc_, 6, endtime);
    sqlite3_bind_int64(stmt_insert_wfdisc_, 7, nsamp);
    sqlite3_bind_double(stmt_insert_wfdisc_, 8, samprate);
    sqlite3_bind_text(stmt_insert_wfdisc_, 9, dir.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt_insert_wfdisc_, 10, dfile.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(stmt_insert_wfdisc_, 11, lddate);
    
    if (sqlite3_step(stmt_insert_wfdisc_) != SQLITE_DONE) {
        setError("Failed to insert wfdisc");
        return -1;
    }
    
    return wfid;
}

std::vector<css30::Event> CSS30Database::queryEvents(double starttime, double endtime,
                                                      double minlat, double maxlat,
                                                      double minlon, double maxlon,
                                                      double minmag, double maxmag) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<css30::Event> results;
    
    if (!db_) return results;
    
    std::string sql = R"(
        SELECT e.evid, e.evname, e.prefor, e.auth, e.lddate
        FROM event e
        JOIN origin o ON e.prefor = o.orid
        WHERE o.time BETWEEN ? AND ?
        AND o.lat BETWEEN ? AND ?
        AND o.lon BETWEEN ? AND ?
    )";
    
    if (minmag > -10 || maxmag < 10) {
        sql += " AND (o.ml BETWEEN ? AND ? OR o.mb BETWEEN ? AND ? OR o.ms BETWEEN ? AND ?)";
    }
    
    sqlite3_stmt* stmt;
    if (sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
        return results;
    }
    
    sqlite3_bind_double(stmt, 1, starttime);
    sqlite3_bind_double(stmt, 2, endtime);
    sqlite3_bind_double(stmt, 3, minlat);
    sqlite3_bind_double(stmt, 4, maxlat);
    sqlite3_bind_double(stmt, 5, minlon);
    sqlite3_bind_double(stmt, 6, maxlon);
    
    if (minmag > -10 || maxmag < 10) {
        for (int i = 0; i < 3; i++) {
            sqlite3_bind_double(stmt, 7 + i*2, minmag);
            sqlite3_bind_double(stmt, 8 + i*2, maxmag);
        }
    }
    
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        css30::Event event;
        event.evid = sqlite3_column_int64(stmt, 0);
        event.evname = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        event.prefor = sqlite3_column_int64(stmt, 2);
        event.auth = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3));
        event.lddate = sqlite3_column_double(stmt, 4);
        results.push_back(event);
    }
    
    sqlite3_finalize(stmt);
    return results;
}

std::vector<css30::Origin> CSS30Database::queryOrigins(int64_t evid) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<css30::Origin> results;
    
    if (!db_) return results;
    
    sqlite3_stmt* stmt;
    const char* sql = "SELECT * FROM origin WHERE evid = ? ORDER BY orid";
    
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return results;
    }
    
    sqlite3_bind_int64(stmt, 1, evid);
    
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        css30::Origin origin;
        origin.lat = sqlite3_column_double(stmt, 0);
        origin.lon = sqlite3_column_double(stmt, 1);
        origin.depth = sqlite3_column_double(stmt, 2);
        origin.time = sqlite3_column_double(stmt, 3);
        origin.orid = sqlite3_column_int64(stmt, 4);
        origin.evid = sqlite3_column_int64(stmt, 5);
        origin.nass = sqlite3_column_int64(stmt, 7);
        origin.ndef = sqlite3_column_int64(stmt, 8);
        origin.ml = sqlite3_column_double(stmt, 20);
        results.push_back(origin);
    }
    
    sqlite3_finalize(stmt);
    return results;
}

std::vector<css30::Arrival> CSS30Database::queryArrivals(int64_t orid) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<css30::Arrival> results;
    
    if (!db_) return results;
    
    sqlite3_stmt* stmt;
    const char* sql = R"(
        SELECT a.* FROM arrival a
        JOIN assoc s ON a.arid = s.arid
        WHERE s.orid = ?
        ORDER BY a.time
    )";
    
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return results;
    }
    
    sqlite3_bind_int64(stmt, 1, orid);
    
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        css30::Arrival arr;
        arr.sta = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        arr.time = sqlite3_column_double(stmt, 1);
        arr.arid = sqlite3_column_int64(stmt, 2);
        arr.iphase = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 7));
        arr.snr = sqlite3_column_double(stmt, 21);
        results.push_back(arr);
    }
    
    sqlite3_finalize(stmt);
    return results;
}

std::vector<css30::Assoc> CSS30Database::queryAssociations(int64_t orid) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<css30::Assoc> results;
    
    if (!db_) return results;
    
    sqlite3_stmt* stmt;
    const char* sql = "SELECT * FROM assoc WHERE orid = ?";
    
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return results;
    }
    
    sqlite3_bind_int64(stmt, 1, orid);
    
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        css30::Assoc assoc;
        assoc.arid = sqlite3_column_int64(stmt, 0);
        assoc.orid = sqlite3_column_int64(stmt, 1);
        assoc.sta = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2));
        assoc.phase = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3));
        assoc.delta = sqlite3_column_double(stmt, 5);
        assoc.timeres = sqlite3_column_double(stmt, 8);
        assoc.wgt = sqlite3_column_double(stmt, 15);
        assoc.vmodel = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 16));
        results.push_back(assoc);
    }
    
    sqlite3_finalize(stmt);
    return results;
}

std::vector<css30::Netmag> CSS30Database::queryMagnitudes(int64_t orid) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<css30::Netmag> results;
    
    if (!db_) return results;
    
    sqlite3_stmt* stmt;
    const char* sql = "SELECT * FROM netmag WHERE orid = ?";
    
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return results;
    }
    
    sqlite3_bind_int64(stmt, 1, orid);
    
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        css30::Netmag mag;
        mag.magid = sqlite3_column_int64(stmt, 0);
        mag.net = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        mag.orid = sqlite3_column_int64(stmt, 2);
        mag.evid = sqlite3_column_int64(stmt, 3);
        mag.magtype = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 4));
        mag.nsta = sqlite3_column_int64(stmt, 5);
        mag.magnitude = sqlite3_column_double(stmt, 6);
        mag.uncertainty = sqlite3_column_double(stmt, 7);
        results.push_back(mag);
    }
    
    sqlite3_finalize(stmt);
    return results;
}

css30::Origin CSS30Database::getPreferredOrigin(int64_t evid) {
    std::lock_guard<std::mutex> lock(mutex_);
    css30::Origin origin;
    
    if (!db_) return origin;
    
    sqlite3_stmt* stmt;
    const char* sql = R"(
        SELECT o.* FROM origin o
        JOIN event e ON e.prefor = o.orid
        WHERE e.evid = ?
    )";
    
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return origin;
    }
    
    sqlite3_bind_int64(stmt, 1, evid);
    
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        origin.lat = sqlite3_column_double(stmt, 0);
        origin.lon = sqlite3_column_double(stmt, 1);
        origin.depth = sqlite3_column_double(stmt, 2);
        origin.time = sqlite3_column_double(stmt, 3);
        origin.orid = sqlite3_column_int64(stmt, 4);
        origin.evid = sqlite3_column_int64(stmt, 5);
    }
    
    sqlite3_finalize(stmt);
    return origin;
}

int64_t CSS30Database::countEvents() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_) return 0;
    
    sqlite3_stmt* stmt;
    int64_t count = 0;
    
    if (sqlite3_prepare_v2(db_, "SELECT COUNT(*) FROM event", -1, &stmt, nullptr) == SQLITE_OK) {
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            count = sqlite3_column_int64(stmt, 0);
        }
        sqlite3_finalize(stmt);
    }
    
    return count;
}

int64_t CSS30Database::countOrigins() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_) return 0;
    
    sqlite3_stmt* stmt;
    int64_t count = 0;
    
    if (sqlite3_prepare_v2(db_, "SELECT COUNT(*) FROM origin", -1, &stmt, nullptr) == SQLITE_OK) {
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            count = sqlite3_column_int64(stmt, 0);
        }
        sqlite3_finalize(stmt);
    }
    
    return count;
}

int64_t CSS30Database::countArrivals() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_) return 0;
    
    sqlite3_stmt* stmt;
    int64_t count = 0;
    
    if (sqlite3_prepare_v2(db_, "SELECT COUNT(*) FROM arrival", -1, &stmt, nullptr) == SQLITE_OK) {
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            count = sqlite3_column_int64(stmt, 0);
        }
        sqlite3_finalize(stmt);
    }
    
    return count;
}

std::pair<double, double> CSS30Database::timeRange() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!db_) return {0, 0};
    
    sqlite3_stmt* stmt;
    std::pair<double, double> range{0, 0};
    
    if (sqlite3_prepare_v2(db_, 
            "SELECT MIN(time), MAX(time) FROM origin", -1, &stmt, nullptr) == SQLITE_OK) {
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            range.first = sqlite3_column_double(stmt, 0);
            range.second = sqlite3_column_double(stmt, 1);
        }
        sqlite3_finalize(stmt);
    }
    
    return range;
}

bool CSS30Database::exportToFlatFiles(const std::string& directory,
                                       double starttime, double endtime) {
    if (!db_) return false;
    
    // Create output files for each table
    std::ofstream event_file(directory + "/event");
    std::ofstream origin_file(directory + "/origin");
    std::ofstream origerr_file(directory + "/origerr");
    std::ofstream arrival_file(directory + "/arrival");
    std::ofstream assoc_file(directory + "/assoc");
    std::ofstream netmag_file(directory + "/netmag");
    
    if (!event_file || !origin_file || !arrival_file) {
        return false;
    }
    
    // Export events
    sqlite3_stmt* stmt;
    std::string sql = "SELECT * FROM event";
    if (starttime > 0) {
        sql = R"(
            SELECT e.* FROM event e
            JOIN origin o ON e.prefor = o.orid
            WHERE o.time BETWEEN ? AND ?
        )";
    }
    
    if (sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        if (starttime > 0) {
            sqlite3_bind_double(stmt, 1, starttime);
            sqlite3_bind_double(stmt, 2, endtime);
        }
        
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            event_file << std::setw(8) << sqlite3_column_int64(stmt, 0) << " "
                       << std::setw(32) << sqlite3_column_text(stmt, 1) << " "
                       << std::setw(8) << sqlite3_column_int64(stmt, 2) << " "
                       << std::setw(15) << sqlite3_column_text(stmt, 3) << "\n";
        }
        sqlite3_finalize(stmt);
    }
    
    // Export origins
    sql = "SELECT * FROM origin";
    if (starttime > 0) {
        sql += " WHERE time BETWEEN ? AND ?";
    }
    
    if (sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        if (starttime > 0) {
            sqlite3_bind_double(stmt, 1, starttime);
            sqlite3_bind_double(stmt, 2, endtime);
        }
        
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            origin_file << std::fixed << std::setprecision(4)
                        << std::setw(9) << sqlite3_column_double(stmt, 0) << " "
                        << std::setw(10) << sqlite3_column_double(stmt, 1) << " "
                        << std::setw(9) << sqlite3_column_double(stmt, 2) << " "
                        << std::setprecision(5)
                        << std::setw(17) << sqlite3_column_double(stmt, 3) << " "
                        << std::setw(8) << sqlite3_column_int64(stmt, 4) << " "
                        << std::setw(8) << sqlite3_column_int64(stmt, 5) << "\n";
        }
        sqlite3_finalize(stmt);
    }
    
    return true;
}

void CSS30Database::beginTransaction() {
    if (db_) {
        executeSQL("BEGIN TRANSACTION");
    }
}

void CSS30Database::commitTransaction() {
    if (db_) {
        executeSQL("COMMIT");
    }
}

void CSS30Database::rollbackTransaction() {
    if (db_) {
        executeSQL("ROLLBACK");
    }
}

bool CSS30Database::initializeSequences() {
    // Get max values for each sequence
    sqlite3_stmt* stmt;
    
    auto getMax = [this, &stmt](const char* table, const char* column, int64_t& value) {
        std::string sql = "SELECT MAX(" + std::string(column) + ") FROM " + table;
        if (sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
            if (sqlite3_step(stmt) == SQLITE_ROW && sqlite3_column_type(stmt, 0) != SQLITE_NULL) {
                value = sqlite3_column_int64(stmt, 0) + 1;
            }
            sqlite3_finalize(stmt);
        }
    };
    
    getMax("event", "evid", next_evid_);
    getMax("origin", "orid", next_orid_);
    getMax("arrival", "arid", next_arid_);
    getMax("netmag", "magid", next_magid_);
    getMax("wfdisc", "wfid", next_wfid_);
    getMax("sitechan", "chanid", next_chanid_);
    getMax("remark", "commid", next_commid_);
    
    return true;
}

bool CSS30Database::prepareStatements() {
    int rc;
    
    // Event insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT INTO event (evid, evname, prefor, auth, lddate) VALUES (?, ?, ?, ?, ?)",
        -1, &stmt_insert_event_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Origin insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT INTO origin (lat, lon, depth, time, orid, evid, jdate, nass, ndef, "
        "etype, dtype, algorithm, auth, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_origin_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Origerr insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT INTO origerr (orid, sxx, syy, szz, stt, sdobs, smajax, sminax, "
        "strike, sdepth, stime, conf, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_origerr_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Arrival insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT INTO arrival (sta, time, arid, jdate, chan, iphase, deltim, amp, "
        "per, snr, qual, auth, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_arrival_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Assoc insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT INTO assoc (arid, orid, sta, phase, delta, seaz, esaz, timeres, "
        "timedef, wgt, vmodel, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_assoc_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Netmag insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT INTO netmag (magid, net, orid, evid, magtype, nsta, magnitude, "
        "uncertainty, auth, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_netmag_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Stamag insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT INTO stamag (magid, sta, orid, evid, delta, magtype, magnitude, "
        "auth, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_stamag_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Site insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT OR REPLACE INTO site (sta, ondate, offdate, lat, lon, elev, staname, "
        "statype, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_site_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Sitechan insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT OR REPLACE INTO sitechan (sta, chan, ondate, chanid, offdate, hang, "
        "vang, descrip, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_sitechan_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Affiliation insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT OR REPLACE INTO affiliation (net, sta, lddate) VALUES (?, ?, ?)",
        -1, &stmt_insert_affiliation_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Wfdisc insert
    rc = sqlite3_prepare_v2(db_,
        "INSERT INTO wfdisc (sta, chan, time, wfid, jdate, endtime, nsamp, samprate, "
        "dir, dfile, lddate) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        -1, &stmt_insert_wfdisc_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    // Update prefor
    rc = sqlite3_prepare_v2(db_,
        "UPDATE event SET prefor = ? WHERE evid = ?",
        -1, &stmt_update_prefor_, nullptr);
    if (rc != SQLITE_OK) return false;
    
    return true;
}

void CSS30Database::finalizeStatements() {
    if (stmt_insert_event_) sqlite3_finalize(stmt_insert_event_);
    if (stmt_insert_origin_) sqlite3_finalize(stmt_insert_origin_);
    if (stmt_insert_origerr_) sqlite3_finalize(stmt_insert_origerr_);
    if (stmt_insert_arrival_) sqlite3_finalize(stmt_insert_arrival_);
    if (stmt_insert_assoc_) sqlite3_finalize(stmt_insert_assoc_);
    if (stmt_insert_netmag_) sqlite3_finalize(stmt_insert_netmag_);
    if (stmt_insert_stamag_) sqlite3_finalize(stmt_insert_stamag_);
    if (stmt_insert_site_) sqlite3_finalize(stmt_insert_site_);
    if (stmt_insert_sitechan_) sqlite3_finalize(stmt_insert_sitechan_);
    if (stmt_insert_affiliation_) sqlite3_finalize(stmt_insert_affiliation_);
    if (stmt_insert_wfdisc_) sqlite3_finalize(stmt_insert_wfdisc_);
    if (stmt_update_prefor_) sqlite3_finalize(stmt_update_prefor_);
    
    stmt_insert_event_ = nullptr;
    stmt_insert_origin_ = nullptr;
    stmt_insert_origerr_ = nullptr;
    stmt_insert_arrival_ = nullptr;
    stmt_insert_assoc_ = nullptr;
    stmt_insert_netmag_ = nullptr;
    stmt_insert_stamag_ = nullptr;
    stmt_insert_site_ = nullptr;
    stmt_insert_sitechan_ = nullptr;
    stmt_insert_affiliation_ = nullptr;
    stmt_insert_wfdisc_ = nullptr;
    stmt_update_prefor_ = nullptr;
}

double CSS30Database::currentLddate() const {
    auto now = std::chrono::system_clock::now();
    return std::chrono::duration_cast<std::chrono::seconds>(
        now.time_since_epoch()).count();
}

int64_t CSS30Database::epochToJdate(double epoch) const {
    time_t t = static_cast<time_t>(epoch);
    struct tm* tm = gmtime(&t);
    if (!tm) return 0;
    
    // Julian date format: YYYYDDD
    return (tm->tm_year + 1900) * 1000 + tm->tm_yday + 1;
}

bool CSS30Database::executeSQL(const char* sql) {
    char* errmsg = nullptr;
    int rc = sqlite3_exec(db_, sql, nullptr, nullptr, &errmsg);
    if (rc != SQLITE_OK) {
        if (errmsg) {
            last_error_ = errmsg;
            sqlite3_free(errmsg);
        }
        return false;
    }
    return true;
}

void CSS30Database::setError(const std::string& context) {
    last_error_ = context + ": " + sqlite3_errmsg(db_);
    std::cerr << "CSS30Database error: " << last_error_ << std::endl;
}

} // namespace realdetect
