#include "realdetect/seedlink/seedlink_client.hpp"
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>
#include <fcntl.h>
#include <poll.h>
#include <cstring>
#include <iostream>
#include <sstream>

namespace realdetect {

SeedLinkClient::SeedLinkClient()
    : socket_fd_(-1)
    , port_(18000)
    , state_(State::Disconnected)
    , running_(false)
    , network_timeout_(30)
    , reconnect_interval_(10)
    , auto_reconnect_(true)
    , packets_received_(0)
    , bytes_received_(0)
{
}

SeedLinkClient::~SeedLinkClient() {
    disconnect();
}

bool SeedLinkClient::connect(const std::string& host, int port) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (state_ != State::Disconnected) {
        disconnect();
    }
    
    host_ = host;
    port_ = port;
    state_ = State::Connecting;
    
    // Resolve hostname
    struct addrinfo hints, *result;
    std::memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;
    
    std::string port_str = std::to_string(port);
    int err = getaddrinfo(host.c_str(), port_str.c_str(), &hints, &result);
    if (err != 0) {
        std::cerr << "Failed to resolve host: " << gai_strerror(err) << std::endl;
        state_ = State::Error;
        return false;
    }
    
    // Create socket
    socket_fd_ = socket(result->ai_family, result->ai_socktype, result->ai_protocol);
    if (socket_fd_ < 0) {
        std::cerr << "Failed to create socket: " << strerror(errno) << std::endl;
        freeaddrinfo(result);
        state_ = State::Error;
        return false;
    }
    
    // Set socket options
    int flags = fcntl(socket_fd_, F_GETFL, 0);
    fcntl(socket_fd_, F_SETFL, flags | O_NONBLOCK);
    
    // Connect
    err = ::connect(socket_fd_, result->ai_addr, result->ai_addrlen);
    freeaddrinfo(result);
    
    if (err < 0 && errno != EINPROGRESS) {
        std::cerr << "Failed to connect: " << strerror(errno) << std::endl;
        close(socket_fd_);
        socket_fd_ = -1;
        state_ = State::Error;
        return false;
    }
    
    // Wait for connection with timeout
    struct pollfd pfd;
    pfd.fd = socket_fd_;
    pfd.events = POLLOUT;
    
    int poll_result = poll(&pfd, 1, network_timeout_ * 1000);
    if (poll_result <= 0) {
        std::cerr << "Connection timeout" << std::endl;
        close(socket_fd_);
        socket_fd_ = -1;
        state_ = State::Error;
        return false;
    }
    
    // Check connection status
    int so_error;
    socklen_t len = sizeof(so_error);
    getsockopt(socket_fd_, SOL_SOCKET, SO_ERROR, &so_error, &len);
    if (so_error != 0) {
        std::cerr << "Connection failed: " << strerror(so_error) << std::endl;
        close(socket_fd_);
        socket_fd_ = -1;
        state_ = State::Error;
        return false;
    }
    
    // Set back to blocking with timeout
    fcntl(socket_fd_, F_SETFL, flags);
    
    struct timeval tv;
    tv.tv_sec = network_timeout_;
    tv.tv_usec = 0;
    setsockopt(socket_fd_, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
    setsockopt(socket_fd_, SOL_SOCKET, SO_SNDTIMEO, &tv, sizeof(tv));
    
    // Read server greeting
    server_info_ = readLine();
    if (server_info_.empty()) {
        std::cerr << "Failed to read server greeting" << std::endl;
        close(socket_fd_);
        socket_fd_ = -1;
        state_ = State::Error;
        return false;
    }
    
    std::cout << "Connected to SeedLink: " << server_info_ << std::endl;
    state_ = State::Connected;
    return true;
}

void SeedLinkClient::disconnect() {
    running_ = false;
    
    if (receive_thread_.joinable()) {
        receive_thread_.join();
    }
    
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (socket_fd_ >= 0) {
        sendCommand("BYE");
        close(socket_fd_);
        socket_fd_ = -1;
    }
    
    state_ = State::Disconnected;
}

std::vector<std::string> SeedLinkClient::availableStations() {
    std::vector<std::string> stations;
    
    if (state_ != State::Connected) return stations;
    
    if (!sendCommand("CAT")) return stations;
    
    std::string line;
    while (!(line = readLine()).empty()) {
        if (line == "END") break;
        stations.push_back(line);
    }
    
    return stations;
}

void SeedLinkClient::addStream(const StreamSelector& selector) {
    selectors_.push_back(selector);
}

void SeedLinkClient::addStream(const std::string& network, const std::string& station,
                                const std::string& selector) {
    addStream(StreamSelector(network, station, selector));
}

void SeedLinkClient::clearStreams() {
    selectors_.clear();
}

bool SeedLinkClient::startStreaming() {
    if (state_ != State::Connected) {
        std::cerr << "Not connected" << std::endl;
        return false;
    }
    
    // Configure streams
    for (const auto& sel : selectors_) {
        std::ostringstream cmd;
        cmd << "STATION " << sel.station << " " << sel.network;
        if (!sendCommand(cmd.str())) return false;
        
        std::string response = readLine();
        if (response.find("OK") == std::string::npos) {
            std::cerr << "Station command failed: " << response << std::endl;
            continue;
        }
        
        cmd.str("");
        cmd << "SELECT " << sel.selector;
        if (!sendCommand(cmd.str())) return false;
        
        response = readLine();
        // SELECT doesn't always return OK
    }
    
    // Start data transfer
    if (!sendCommand("DATA")) return false;
    
    std::string response = readLine();
    if (response.find("OK") == std::string::npos && 
        response.find("ERROR") != std::string::npos) {
        std::cerr << "DATA command failed: " << response << std::endl;
        return false;
    }
    
    state_ = State::Streaming;
    running_ = true;
    
    receive_thread_ = std::thread(&SeedLinkClient::receiveLoop, this);
    
    return true;
}

void SeedLinkClient::stopStreaming() {
    running_ = false;
    
    if (receive_thread_.joinable()) {
        receive_thread_.join();
    }
    
    state_ = State::Connected;
}

bool SeedLinkClient::sendCommand(const std::string& cmd) {
    if (socket_fd_ < 0) return false;
    
    std::string full_cmd = cmd + "\r\n";
    ssize_t sent = send(socket_fd_, full_cmd.c_str(), full_cmd.size(), 0);
    
    return sent == static_cast<ssize_t>(full_cmd.size());
}

std::string SeedLinkClient::readLine() {
    std::string line;
    char c;
    
    while (recv(socket_fd_, &c, 1, 0) == 1) {
        if (c == '\r') continue;
        if (c == '\n') break;
        line += c;
        
        if (line.size() > 1024) break;  // Prevent overflow
    }
    
    return line;
}

bool SeedLinkClient::readPacket(SeedLinkPacket& packet) {
    // SeedLink packet format:
    // "SL" (2 bytes) + sequence (6 bytes) + miniSEED record (512 bytes)
    
    uint8_t header[8];
    ssize_t received = 0;
    
    while (received < 8) {
        ssize_t r = recv(socket_fd_, header + received, 8 - received, 0);
        if (r <= 0) {
            if (errno == EAGAIN || errno == EWOULDBLOCK) continue;
            return false;
        }
        received += r;
    }
    
    // Check for "SL" header
    if (header[0] != 'S' || header[1] != 'L') {
        // Might be an info packet
        if (header[0] == 'I' && header[1] == 'N' && 
            header[2] == 'F' && header[3] == 'O') {
            // Read and discard info packet
            std::string info = readLine();
            return false;
        }
        return false;
    }
    
    packet.sequence = std::string(reinterpret_cast<char*>(header + 2), 6);
    
    // Read miniSEED record (512 bytes default)
    std::vector<uint8_t> mseed(512);
    received = 0;
    
    while (received < 512) {
        ssize_t r = recv(socket_fd_, mseed.data() + received, 512 - received, 0);
        if (r <= 0) {
            if (errno == EAGAIN || errno == EWOULDBLOCK) continue;
            return false;
        }
        received += r;
    }
    
    bytes_received_ += 520;  // Header + record
    
    // Parse miniSEED
    MiniSeedRecord record;
    if (!record.parse(mseed.data(), mseed.size())) {
        return false;
    }
    
    packet.stream_id = record.streamId();
    packet.start_time = record.startTime();
    packet.sample_rate = record.sampleRate();
    packet.samples = record.samples();
    
    packets_received_++;
    
    return true;
}

void SeedLinkClient::receiveLoop() {
    while (running_) {
        SeedLinkPacket packet;
        
        if (readPacket(packet)) {
            if (packet_callback_) {
                packet_callback_(packet);
            }
            
            if (waveform_callback_) {
                auto wf = std::make_shared<Waveform>(
                    packet.stream_id, packet.sample_rate, packet.start_time);
                wf->data() = packet.samples;
                waveform_callback_(wf);
            }
        } else {
            // Check if still connected
            if (!running_) break;
            
            // Possible disconnect
            int error = 0;
            socklen_t len = sizeof(error);
            int result = getsockopt(socket_fd_, SOL_SOCKET, SO_ERROR, &error, &len);
            
            if (result != 0 || error != 0) {
                handleDisconnect();
            }
        }
    }
}

void SeedLinkClient::handleDisconnect() {
    state_ = State::Disconnected;
    
    if (error_callback_) {
        error_callback_("Disconnected from server");
    }
    
    if (auto_reconnect_) {
        attemptReconnect();
    }
}

bool SeedLinkClient::attemptReconnect() {
    std::cerr << "Attempting to reconnect..." << std::endl;
    
    for (int attempt = 0; attempt < 10 && running_; attempt++) {
        std::this_thread::sleep_for(std::chrono::seconds(reconnect_interval_));
        
        if (!running_) break;
        
        if (connect(host_, port_)) {
            // Reconfigure streams
            for (const auto& sel : selectors_) {
                std::ostringstream cmd;
                cmd << "STATION " << sel.station << " " << sel.network;
                sendCommand(cmd.str());
                readLine();
                
                cmd.str("");
                cmd << "SELECT " << sel.selector;
                sendCommand(cmd.str());
                readLine();
            }
            
            sendCommand("DATA");
            readLine();
            
            state_ = State::Streaming;
            return true;
        }
    }
    
    return false;
}

} // namespace realdetect
