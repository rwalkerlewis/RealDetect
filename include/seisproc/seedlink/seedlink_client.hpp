#pragma once

#include "../core/types.hpp"
#include "../core/waveform.hpp"
#include "../core/miniseed.hpp"
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <condition_variable>
#include <memory>

namespace seisproc {

/**
 * SeedLink packet from server
 */
struct SeedLinkPacket {
    std::string sequence;     // 6-character sequence number
    StreamID stream_id;
    TimePoint start_time;
    double sample_rate;
    SampleVector samples;
    
    SeedLinkPacket() : sample_rate(0) {}
};

/**
 * SeedLink stream selector
 */
struct StreamSelector {
    std::string network;      // Network code or wildcard (*)
    std::string station;      // Station code or wildcard
    std::string selector;     // SEED selector (e.g., "BH?" for all BH channels)
    
    StreamSelector() = default;
    StreamSelector(const std::string& net, const std::string& sta, 
                   const std::string& sel = "???")
        : network(net), station(sta), selector(sel) {}
    
    std::string toString() const {
        return network + "_" + station + ":" + selector;
    }
};

/**
 * Callback type for received packets
 */
using PacketCallback = std::function<void(const SeedLinkPacket&)>;
using WaveformCallback = std::function<void(WaveformPtr)>;
using ErrorCallback = std::function<void(const std::string&)>;

/**
 * SeedLinkClient - Connect to SeedLink servers for real-time data
 */
class SeedLinkClient {
public:
    enum class State {
        Disconnected,
        Connecting,
        Connected,
        Streaming,
        Error
    };
    
    SeedLinkClient();
    ~SeedLinkClient();
    
    // Connection
    bool connect(const std::string& host, int port = 18000);
    void disconnect();
    bool isConnected() const { return state_ == State::Connected || 
                                       state_ == State::Streaming; }
    State state() const { return state_; }
    
    // Server info
    std::string serverInfo() const { return server_info_; }
    std::vector<std::string> availableStations();
    
    // Stream selection
    void addStream(const StreamSelector& selector);
    void addStream(const std::string& network, const std::string& station,
                   const std::string& selector = "???");
    void clearStreams();
    
    // Start/stop streaming
    bool startStreaming();
    void stopStreaming();
    bool isStreaming() const { return state_ == State::Streaming; }
    
    // Callbacks
    void setPacketCallback(PacketCallback cb) { packet_callback_ = std::move(cb); }
    void setWaveformCallback(WaveformCallback cb) { waveform_callback_ = std::move(cb); }
    void setErrorCallback(ErrorCallback cb) { error_callback_ = std::move(cb); }
    
    // Configuration
    void setNetworkTimeout(int seconds) { network_timeout_ = seconds; }
    void setReconnectInterval(int seconds) { reconnect_interval_ = seconds; }
    void setAutoReconnect(bool enable) { auto_reconnect_ = enable; }
    
    // Statistics
    size_t packetsReceived() const { return packets_received_; }
    size_t bytesReceived() const { return bytes_received_; }

private:
    int socket_fd_;
    std::string host_;
    int port_;
    State state_;
    std::string server_info_;
    
    std::vector<StreamSelector> selectors_;
    
    std::thread receive_thread_;
    std::atomic<bool> running_;
    std::mutex mutex_;
    
    PacketCallback packet_callback_;
    WaveformCallback waveform_callback_;
    ErrorCallback error_callback_;
    
    int network_timeout_;
    int reconnect_interval_;
    bool auto_reconnect_;
    
    std::atomic<size_t> packets_received_;
    std::atomic<size_t> bytes_received_;
    
    // Protocol methods
    bool sendCommand(const std::string& cmd);
    std::string readLine();
    bool readPacket(SeedLinkPacket& packet);
    void receiveLoop();
    
    // Reconnection handling
    void handleDisconnect();
    bool attemptReconnect();
};

/**
 * StreamBuffer - Buffer incoming data per stream
 */
class StreamBuffer {
public:
    StreamBuffer(size_t max_seconds = 300, double sample_rate = 100.0);
    
    void append(const SeedLinkPacket& packet);
    void append(WaveformPtr waveform);
    
    // Get continuous waveform for time range
    WaveformPtr getWaveform(TimePoint start, TimePoint end) const;
    
    // Get latest N seconds
    WaveformPtr getLatest(double seconds) const;
    
    // Stream info
    const StreamID& streamId() const { return stream_id_; }
    TimePoint startTime() const { return start_time_; }
    TimePoint endTime() const { return end_time_; }
    double sampleRate() const { return sample_rate_; }
    
    // Clear buffer
    void clear();
    void trimTo(TimePoint earliest);

private:
    StreamID stream_id_;
    double sample_rate_;
    size_t max_samples_;
    TimePoint start_time_;
    TimePoint end_time_;
    SampleVector buffer_;
    mutable std::mutex mutex_;
};

/**
 * MultiStreamBuffer - Manage buffers for multiple streams
 */
class MultiStreamBuffer {
public:
    MultiStreamBuffer(size_t max_seconds = 300);
    
    void addPacket(const SeedLinkPacket& packet);
    void addWaveform(WaveformPtr waveform);
    
    StreamBuffer* getBuffer(const StreamID& id);
    const StreamBuffer* getBuffer(const StreamID& id) const;
    
    std::vector<StreamID> streams() const;
    
    void clear();
    void trimAll(TimePoint earliest);

private:
    size_t max_seconds_;
    std::map<std::string, std::unique_ptr<StreamBuffer>> buffers_;
    mutable std::mutex mutex_;
};

} // namespace seisproc
