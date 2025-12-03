#pragma once

/**
 * CSS3.0 Database Schema Definitions
 * 
 * Implementation of the Center for Seismic Studies Schema 3.0 standard
 * for storing seismic event catalogs and related data.
 * 
 * Reference: Anderson et al., "Database Schema for Seismological Data"
 * Center for Seismic Studies, Science Applications International Corporation
 */

#include <string>
#include <cstdint>
#include <ctime>

namespace seisproc {
namespace css30 {

// CSS3.0 null values
constexpr int64_t NULL_INT = -1;
constexpr double NULL_DOUBLE = -999.0;
constexpr double NULL_TIME = -9999999999.999;
constexpr const char* NULL_STRING = "-";

// Maximum field lengths (CSS3.0 standard)
constexpr size_t MAX_STA = 6;      // Station code
constexpr size_t MAX_CHAN = 8;     // Channel code
constexpr size_t MAX_NET = 8;      // Network code
constexpr size_t MAX_AUTH = 15;    // Author
constexpr size_t MAX_IPHASE = 8;   // Phase type
constexpr size_t MAX_DIR = 64;     // Directory path
constexpr size_t MAX_DFILE = 32;   // Data file name
constexpr size_t MAX_ALGO = 15;    // Algorithm name
constexpr size_t MAX_DTYPE = 4;    // Data type
constexpr size_t MAX_ETYPE = 7;    // Event type

/**
 * EVENT table - Earthquake event identification
 */
struct Event {
    int64_t evid;           // Event identifier (unique)
    std::string evname;     // Event name (optional)
    int64_t prefor;         // Preferred origin ID
    std::string auth;       // Source/author
    int64_t commid;         // Comment identifier
    double lddate;          // Load date
    
    Event() : evid(NULL_INT), prefor(NULL_INT), auth("-"), 
              commid(NULL_INT), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS event (
            evid        INTEGER PRIMARY KEY,
            evname      TEXT,
            prefor      INTEGER,
            auth        TEXT DEFAULT '-',
            commid      INTEGER DEFAULT -1,
            lddate      REAL
        )
    )";
};

/**
 * ORIGIN table - Hypocenter location solutions
 */
struct Origin {
    double lat;             // Latitude (degrees)
    double lon;             // Longitude (degrees)
    double depth;           // Depth (km)
    double time;            // Origin time (epoch)
    int64_t orid;           // Origin identifier (unique)
    int64_t evid;           // Event identifier
    double jdate;           // Julian date (YYYYDDD)
    int64_t nass;           // Number associated phases
    int64_t ndef;           // Number defining phases
    int64_t ndp;            // Number depth phases
    int64_t grn;            // Geographic region number
    int64_t srn;            // Seismic region number
    std::string etype;      // Event type (eq, ex, qb, etc)
    double depdp;           // Depth from depth phases
    std::string dtype;      // Depth method
    double mb;              // Body wave magnitude
    int64_t mbid;           // Mb orid
    double ms;              // Surface wave magnitude
    int64_t msid;           // Ms orid
    double ml;              // Local magnitude
    int64_t mlid;           // ML orid
    std::string algorithm;  // Location algorithm
    std::string auth;       // Source/author
    int64_t commid;         // Comment identifier
    double lddate;          // Load date
    
    Origin() : lat(NULL_DOUBLE), lon(NULL_DOUBLE), depth(NULL_DOUBLE),
               time(NULL_TIME), orid(NULL_INT), evid(NULL_INT), jdate(NULL_INT),
               nass(NULL_INT), ndef(NULL_INT), ndp(NULL_INT), grn(NULL_INT),
               srn(NULL_INT), etype("-"), depdp(NULL_DOUBLE), dtype("-"),
               mb(NULL_DOUBLE), mbid(NULL_INT), ms(NULL_DOUBLE), msid(NULL_INT),
               ml(NULL_DOUBLE), mlid(NULL_INT), algorithm("-"), auth("-"),
               commid(NULL_INT), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS origin (
            lat         REAL,
            lon         REAL,
            depth       REAL,
            time        REAL,
            orid        INTEGER PRIMARY KEY,
            evid        INTEGER,
            jdate       INTEGER,
            nass        INTEGER DEFAULT -1,
            ndef        INTEGER DEFAULT -1,
            ndp         INTEGER DEFAULT -1,
            grn         INTEGER DEFAULT -1,
            srn         INTEGER DEFAULT -1,
            etype       TEXT DEFAULT '-',
            depdp       REAL DEFAULT -999.0,
            dtype       TEXT DEFAULT '-',
            mb          REAL DEFAULT -999.0,
            mbid        INTEGER DEFAULT -1,
            ms          REAL DEFAULT -999.0,
            msid        INTEGER DEFAULT -1,
            ml          REAL DEFAULT -999.0,
            mlid        INTEGER DEFAULT -1,
            algorithm   TEXT DEFAULT '-',
            auth        TEXT DEFAULT '-',
            commid      INTEGER DEFAULT -1,
            lddate      REAL,
            FOREIGN KEY (evid) REFERENCES event(evid)
        )
    )";
};

/**
 * ORIGERR table - Origin error information
 */
struct Origerr {
    int64_t orid;           // Origin identifier
    double sxx;             // Covariance matrix elements
    double syy;
    double szz;
    double stt;
    double sxy;
    double sxz;
    double syz;
    double stx;
    double sty;
    double stz;
    double sdobs;           // Standard error of observations
    double smajax;          // Semi-major axis of error ellipse
    double sminax;          // Semi-minor axis of error ellipse
    double strike;          // Strike of major axis
    double sdepth;          // Depth error
    double stime;           // Origin time error
    double conf;            // Confidence level
    int64_t commid;
    double lddate;
    
    Origerr() : orid(NULL_INT), sxx(NULL_DOUBLE), syy(NULL_DOUBLE),
                szz(NULL_DOUBLE), stt(NULL_DOUBLE), sxy(NULL_DOUBLE),
                sxz(NULL_DOUBLE), syz(NULL_DOUBLE), stx(NULL_DOUBLE),
                sty(NULL_DOUBLE), stz(NULL_DOUBLE), sdobs(NULL_DOUBLE),
                smajax(NULL_DOUBLE), sminax(NULL_DOUBLE), strike(NULL_DOUBLE),
                sdepth(NULL_DOUBLE), stime(NULL_DOUBLE), conf(0.9),
                commid(NULL_INT), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS origerr (
            orid        INTEGER PRIMARY KEY,
            sxx         REAL DEFAULT -999.0,
            syy         REAL DEFAULT -999.0,
            szz         REAL DEFAULT -999.0,
            stt         REAL DEFAULT -999.0,
            sxy         REAL DEFAULT -999.0,
            sxz         REAL DEFAULT -999.0,
            syz         REAL DEFAULT -999.0,
            stx         REAL DEFAULT -999.0,
            sty         REAL DEFAULT -999.0,
            stz         REAL DEFAULT -999.0,
            sdobs       REAL DEFAULT -999.0,
            smajax      REAL DEFAULT -999.0,
            sminax      REAL DEFAULT -999.0,
            strike      REAL DEFAULT -999.0,
            sdepth      REAL DEFAULT -999.0,
            stime       REAL DEFAULT -999.0,
            conf        REAL DEFAULT 0.9,
            commid      INTEGER DEFAULT -1,
            lddate      REAL,
            FOREIGN KEY (orid) REFERENCES origin(orid)
        )
    )";
};

/**
 * ARRIVAL table - Phase arrival data
 */
struct Arrival {
    std::string sta;        // Station code
    double time;            // Arrival time (epoch)
    int64_t arid;           // Arrival identifier
    double jdate;           // Julian date
    int64_t stassid;        // Station association id
    int64_t chanid;         // Channel identifier
    std::string chan;       // Channel code
    std::string iphase;     // Initial phase identification
    std::string stype;      // Signal type (l=local, r=regional, t=teleseismic)
    double deltim;          // Time uncertainty
    double azimuth;         // Observed azimuth
    double delaz;           // Azimuth uncertainty
    double slow;            // Observed slowness
    double delslo;          // Slowness uncertainty
    double ema;             // Emergence angle
    double rect;            // Rectilinearity
    double amp;             // Amplitude
    double per;             // Period
    double logat;           // Log amplitude/period
    std::string clip;       // Clipped flag
    std::string fm;         // First motion (c=compression, d=dilatation)
    double snr;             // Signal-to-noise ratio
    std::string qual;       // Quality (i=impulsive, e=emergent)
    std::string auth;       // Author
    int64_t commid;
    double lddate;
    
    Arrival() : sta("-"), time(NULL_TIME), arid(NULL_INT), jdate(NULL_INT),
                stassid(NULL_INT), chanid(NULL_INT), chan("-"), iphase("-"),
                stype("-"), deltim(NULL_DOUBLE), azimuth(NULL_DOUBLE),
                delaz(NULL_DOUBLE), slow(NULL_DOUBLE), delslo(NULL_DOUBLE),
                ema(NULL_DOUBLE), rect(NULL_DOUBLE), amp(NULL_DOUBLE),
                per(NULL_DOUBLE), logat(NULL_DOUBLE), clip("-"), fm("-"),
                snr(NULL_DOUBLE), qual("-"), auth("-"), commid(NULL_INT),
                lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS arrival (
            sta         TEXT,
            time        REAL,
            arid        INTEGER PRIMARY KEY,
            jdate       INTEGER,
            stassid     INTEGER DEFAULT -1,
            chanid      INTEGER DEFAULT -1,
            chan        TEXT DEFAULT '-',
            iphase      TEXT DEFAULT '-',
            stype       TEXT DEFAULT '-',
            deltim      REAL DEFAULT -999.0,
            azimuth     REAL DEFAULT -999.0,
            delaz       REAL DEFAULT -999.0,
            slow        REAL DEFAULT -999.0,
            delslo      REAL DEFAULT -999.0,
            ema         REAL DEFAULT -999.0,
            rect        REAL DEFAULT -999.0,
            amp         REAL DEFAULT -999.0,
            per         REAL DEFAULT -999.0,
            logat       REAL DEFAULT -999.0,
            clip        TEXT DEFAULT '-',
            fm          TEXT DEFAULT '-',
            snr         REAL DEFAULT -999.0,
            qual        TEXT DEFAULT '-',
            auth        TEXT DEFAULT '-',
            commid      INTEGER DEFAULT -1,
            lddate      REAL
        )
    )";
};

/**
 * ASSOC table - Arrival-Origin association
 */
struct Assoc {
    int64_t arid;           // Arrival identifier
    int64_t orid;           // Origin identifier
    std::string sta;        // Station code
    std::string phase;      // Associated phase type
    double belief;          // Phase identification confidence
    double delta;           // Source-receiver distance (degrees)
    double seaz;            // Station-to-event azimuth
    double esaz;            // Event-to-station azimuth
    double timeres;         // Time residual
    std::string timedef;    // Time defining (d=defining, n=non-defining)
    double azres;           // Azimuth residual
    std::string azdef;      // Azimuth defining
    double slores;          // Slowness residual
    std::string slodef;     // Slowness defining
    double emares;          // Emergence angle residual
    double wgt;             // Weight
    std::string vmodel;     // Velocity model name
    int64_t commid;
    double lddate;
    
    Assoc() : arid(NULL_INT), orid(NULL_INT), sta("-"), phase("-"),
              belief(NULL_DOUBLE), delta(NULL_DOUBLE), seaz(NULL_DOUBLE),
              esaz(NULL_DOUBLE), timeres(NULL_DOUBLE), timedef("-"),
              azres(NULL_DOUBLE), azdef("-"), slores(NULL_DOUBLE),
              slodef("-"), emares(NULL_DOUBLE), wgt(1.0), vmodel("-"),
              commid(NULL_INT), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS assoc (
            arid        INTEGER,
            orid        INTEGER,
            sta         TEXT,
            phase       TEXT,
            belief      REAL DEFAULT -999.0,
            delta       REAL DEFAULT -999.0,
            seaz        REAL DEFAULT -999.0,
            esaz        REAL DEFAULT -999.0,
            timeres     REAL DEFAULT -999.0,
            timedef     TEXT DEFAULT 'd',
            azres       REAL DEFAULT -999.0,
            azdef       TEXT DEFAULT '-',
            slores      REAL DEFAULT -999.0,
            slodef      TEXT DEFAULT '-',
            emares      REAL DEFAULT -999.0,
            wgt         REAL DEFAULT 1.0,
            vmodel      TEXT DEFAULT '-',
            commid      INTEGER DEFAULT -1,
            lddate      REAL,
            PRIMARY KEY (arid, orid),
            FOREIGN KEY (arid) REFERENCES arrival(arid),
            FOREIGN KEY (orid) REFERENCES origin(orid)
        )
    )";
};

/**
 * NETMAG table - Network magnitude
 */
struct Netmag {
    int64_t magid;          // Magnitude identifier
    std::string net;        // Network code
    int64_t orid;           // Origin identifier
    int64_t evid;           // Event identifier
    std::string magtype;    // Magnitude type (mb, ms, ml, mw)
    int64_t nsta;           // Number of stations
    double magnitude;       // Network magnitude value
    double uncertainty;     // Magnitude uncertainty
    std::string auth;
    int64_t commid;
    double lddate;
    
    Netmag() : magid(NULL_INT), net("-"), orid(NULL_INT), evid(NULL_INT),
               magtype("-"), nsta(NULL_INT), magnitude(NULL_DOUBLE),
               uncertainty(NULL_DOUBLE), auth("-"), commid(NULL_INT),
               lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS netmag (
            magid       INTEGER PRIMARY KEY,
            net         TEXT DEFAULT '-',
            orid        INTEGER,
            evid        INTEGER,
            magtype     TEXT,
            nsta        INTEGER DEFAULT -1,
            magnitude   REAL,
            uncertainty REAL DEFAULT -999.0,
            auth        TEXT DEFAULT '-',
            commid      INTEGER DEFAULT -1,
            lddate      REAL,
            FOREIGN KEY (orid) REFERENCES origin(orid),
            FOREIGN KEY (evid) REFERENCES event(evid)
        )
    )";
};

/**
 * STAMAG table - Station magnitude
 */
struct Stamag {
    int64_t magid;          // Magnitude identifier
    int64_t ampid;          // Amplitude identifier
    std::string sta;        // Station code
    int64_t arid;           // Arrival identifier
    int64_t orid;           // Origin identifier
    int64_t evid;           // Event identifier
    std::string phase;      // Associated phase
    double delta;           // Distance (degrees)
    std::string magtype;    // Magnitude type
    double magnitude;       // Station magnitude value
    double uncertainty;
    std::string auth;
    int64_t commid;
    double lddate;
    
    Stamag() : magid(NULL_INT), ampid(NULL_INT), sta("-"), arid(NULL_INT),
               orid(NULL_INT), evid(NULL_INT), phase("-"), delta(NULL_DOUBLE),
               magtype("-"), magnitude(NULL_DOUBLE), uncertainty(NULL_DOUBLE),
               auth("-"), commid(NULL_INT), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS stamag (
            magid       INTEGER,
            ampid       INTEGER DEFAULT -1,
            sta         TEXT,
            arid        INTEGER DEFAULT -1,
            orid        INTEGER,
            evid        INTEGER,
            phase       TEXT DEFAULT '-',
            delta       REAL DEFAULT -999.0,
            magtype     TEXT,
            magnitude   REAL,
            uncertainty REAL DEFAULT -999.0,
            auth        TEXT DEFAULT '-',
            commid      INTEGER DEFAULT -1,
            lddate      REAL,
            PRIMARY KEY (magid, sta),
            FOREIGN KEY (orid) REFERENCES origin(orid),
            FOREIGN KEY (evid) REFERENCES event(evid)
        )
    )";
};

/**
 * SITE table - Station location
 */
struct Site {
    std::string sta;        // Station code
    double ondate;          // Start date (Julian)
    double offdate;         // End date (Julian, -1 for open)
    double lat;             // Latitude
    double lon;             // Longitude
    double elev;            // Elevation (km)
    std::string staname;    // Station name
    std::string statype;    // Station type (ss=single, ar=array)
    std::string refsta;     // Reference station
    double dnorth;          // Offset north
    double deast;           // Offset east
    double lddate;
    
    Site() : sta("-"), ondate(NULL_INT), offdate(-1), lat(NULL_DOUBLE),
             lon(NULL_DOUBLE), elev(NULL_DOUBLE), staname("-"), statype("ss"),
             refsta("-"), dnorth(0.0), deast(0.0), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS site (
            sta         TEXT,
            ondate      REAL,
            offdate     REAL DEFAULT -1,
            lat         REAL,
            lon         REAL,
            elev        REAL,
            staname     TEXT DEFAULT '-',
            statype     TEXT DEFAULT 'ss',
            refsta      TEXT DEFAULT '-',
            dnorth      REAL DEFAULT 0.0,
            deast       REAL DEFAULT 0.0,
            lddate      REAL,
            PRIMARY KEY (sta, ondate)
        )
    )";
};

/**
 * SITECHAN table - Channel information
 */
struct Sitechan {
    std::string sta;        // Station code
    std::string chan;       // Channel code
    double ondate;          // Start date (Julian)
    int64_t chanid;         // Channel identifier
    double offdate;         // End date (Julian)
    std::string ctype;      // Channel type (n=normal, b=beam, i=inferred)
    double edepth;          // Emplacement depth (km)
    double hang;            // Horizontal angle (degrees from north)
    double vang;            // Vertical angle (degrees from horizontal)
    std::string descrip;    // Description
    double lddate;
    
    Sitechan() : sta("-"), chan("-"), ondate(NULL_INT), chanid(NULL_INT),
                 offdate(-1), ctype("n"), edepth(0.0), hang(NULL_DOUBLE),
                 vang(NULL_DOUBLE), descrip("-"), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS sitechan (
            sta         TEXT,
            chan        TEXT,
            ondate      REAL,
            chanid      INTEGER PRIMARY KEY,
            offdate     REAL DEFAULT -1,
            ctype       TEXT DEFAULT 'n',
            edepth      REAL DEFAULT 0.0,
            hang        REAL DEFAULT -999.0,
            vang        REAL DEFAULT -999.0,
            descrip     TEXT DEFAULT '-',
            lddate      REAL
        )
    )";
};

/**
 * AFFILIATION table - Network-station affiliation
 */
struct Affiliation {
    std::string net;        // Network code
    std::string sta;        // Station code
    double lddate;
    
    Affiliation() : net("-"), sta("-"), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS affiliation (
            net         TEXT,
            sta         TEXT,
            lddate      REAL,
            PRIMARY KEY (net, sta)
        )
    )";
};

/**
 * WFDISC table - Waveform descriptor
 */
struct Wfdisc {
    std::string sta;        // Station code
    std::string chan;       // Channel code
    double time;            // Start time (epoch)
    int64_t wfid;           // Waveform identifier
    int64_t chanid;         // Channel identifier
    double jdate;           // Julian date
    double endtime;         // End time (epoch)
    int64_t nsamp;          // Number of samples
    double samprate;        // Sample rate (Hz)
    double calib;           // Calibration factor
    double calper;          // Calibration period
    std::string instype;    // Instrument type
    std::string segtype;    // Segment type (o=original, s=segmented)
    std::string datatype;   // Data type (s4, i4, f4, etc)
    std::string clip;       // Clipping flag
    std::string dir;        // Directory path
    std::string dfile;      // Data file name
    int64_t foff;           // File offset (bytes)
    int64_t commid;
    double lddate;
    
    Wfdisc() : sta("-"), chan("-"), time(NULL_TIME), wfid(NULL_INT),
               chanid(NULL_INT), jdate(NULL_INT), endtime(NULL_TIME),
               nsamp(NULL_INT), samprate(NULL_DOUBLE), calib(1.0),
               calper(1.0), instype("-"), segtype("o"), datatype("f4"),
               clip("-"), dir("-"), dfile("-"), foff(0), commid(NULL_INT),
               lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS wfdisc (
            sta         TEXT,
            chan        TEXT,
            time        REAL,
            wfid        INTEGER PRIMARY KEY,
            chanid      INTEGER DEFAULT -1,
            jdate       INTEGER,
            endtime     REAL,
            nsamp       INTEGER,
            samprate    REAL,
            calib       REAL DEFAULT 1.0,
            calper      REAL DEFAULT 1.0,
            instype     TEXT DEFAULT '-',
            segtype     TEXT DEFAULT 'o',
            datatype    TEXT DEFAULT 'f4',
            clip        TEXT DEFAULT '-',
            dir         TEXT,
            dfile       TEXT,
            foff        INTEGER DEFAULT 0,
            commid      INTEGER DEFAULT -1,
            lddate      REAL
        )
    )";
};

/**
 * REMARK table - Comment storage
 */
struct Remark {
    int64_t commid;         // Comment identifier
    int64_t lineno;         // Line number
    std::string remark;     // Comment text
    double lddate;
    
    Remark() : commid(NULL_INT), lineno(1), remark("-"), lddate(NULL_TIME) {}
    
    static constexpr const char* CREATE_SQL = R"(
        CREATE TABLE IF NOT EXISTS remark (
            commid      INTEGER,
            lineno      INTEGER,
            remark      TEXT,
            lddate      REAL,
            PRIMARY KEY (commid, lineno)
        )
    )";
};

// Index creation SQL for performance
constexpr const char* CREATE_INDICES_SQL = R"(
    CREATE INDEX IF NOT EXISTS idx_origin_evid ON origin(evid);
    CREATE INDEX IF NOT EXISTS idx_origin_time ON origin(time);
    CREATE INDEX IF NOT EXISTS idx_origin_latlon ON origin(lat, lon);
    CREATE INDEX IF NOT EXISTS idx_arrival_sta ON arrival(sta);
    CREATE INDEX IF NOT EXISTS idx_arrival_time ON arrival(time);
    CREATE INDEX IF NOT EXISTS idx_assoc_orid ON assoc(orid);
    CREATE INDEX IF NOT EXISTS idx_netmag_orid ON netmag(orid);
    CREATE INDEX IF NOT EXISTS idx_netmag_evid ON netmag(evid);
    CREATE INDEX IF NOT EXISTS idx_stamag_orid ON stamag(orid);
    CREATE INDEX IF NOT EXISTS idx_site_sta ON site(sta);
    CREATE INDEX IF NOT EXISTS idx_wfdisc_sta_time ON wfdisc(sta, time);
)";

} // namespace css30
} // namespace seisproc
