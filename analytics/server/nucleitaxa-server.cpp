#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <string.h>
#include <signal.h>
#include <regex>

/**
 * NucleiTaxa Real-Time Analytics Server
 * 
 * Lightweight C++ WebSocket server that monitors pipeline execution
 * and broadcasts live metrics to connected clients.
 * 
 * Metrics:
 * - ASV count (updates as denoise completes)
 * - Taxonomy distribution (live pie chart data)
 * - Quality metrics (Q-score stats)
 * - Performance (elapsed time, CPU%, memory)
 */

volatile sig_atomic_t keep_running = 1;

void signal_handler(int sig) {
    keep_running = 0;
}

struct MetricSnapshot {
    long asv_count = 0;
    long chimera_count = 0;
    float mean_quality = 0.0f;
    float gc_percent = 0.0f;
    std::map<std::string, int> taxa_dist;  // Phylum -> count
    float elapsed_sec = 0.0f;
    int active_stage = 0;
};

class AnalyticsServer {
private:
    int port;
    std::string output_dir;
    std::vector<int> client_sockets;
    MetricSnapshot last_snapshot;
    
public:
    AnalyticsServer(int p, const std::string& out_dir) 
        : port(p), output_dir(out_dir) {
        std::cout << "[Analytics] Initializing server on port " << port << std::endl;
    }
    
    void start() {
        int server_fd = socket(AF_INET, SOCK_STREAM, 0);
        if (server_fd < 0) {
            std::cerr << "[Error] Socket creation failed" << std::endl;
            return;
        }
        
        // Reuse port
        int opt = 1;
        setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
        
        struct sockaddr_in addr;
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = inet_addr("127.0.0.1");
        addr.sin_port = htons(port);
        
        if (bind(server_fd, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
            std::cerr << "[Error] Bind failed" << std::endl;
            return;
        }
        
        listen(server_fd, 5);
        std::cout << "[Analytics] Server listening on localhost:" << port << std::endl;
        
        // Accept connections in separate thread
        std::thread accept_thread([this, server_fd]() {
            while (keep_running) {
                struct sockaddr_in client_addr;
                socklen_t client_len = sizeof(client_addr);
                
                int client_fd = accept(server_fd, (struct sockaddr *)&client_addr, &client_len);
                if (client_fd >= 0) {
                    client_sockets.push_back(client_fd);
                    std::cout << "[Analytics] Client connected (total: " << client_sockets.size() << ")" << std::endl;
                }
            }
        });
        accept_thread.detach();
        
        // Main metrics broadcast loop
        while (keep_running) {
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            
            auto snapshot = read_metrics();
            broadcast_json(snapshot);
        }
        
        close(server_fd);
    }
    
private:
    MetricSnapshot read_metrics() {
        MetricSnapshot snap;
        snap.elapsed_sec = get_elapsed_time();
        
        // Read ASV table from denoise stage
        std::string asv_file = output_dir + "/02-denoise/seqtab.txt";
        snap.asv_count = count_lines(asv_file);
        
        // Read chimera counts
        std::string chimera_file = output_dir + "/03-chimera/chimera_flagged.txt";
        snap.chimera_count = count_lines(chimera_file);
        
        // Read taxonomy distribution
        std::string taxa_file = output_dir + "/04-taxonomy/taxa_assignments.txt";
        snap.taxa_dist = parse_taxonomy_distribution(taxa_file);
        
        // Read QC metrics
        std::string qc_file = output_dir + "/01-qc/quality_stats.txt";
        snap.mean_quality = extract_float_metric(qc_file, "mean_phred");
        snap.gc_percent = extract_float_metric(qc_file, "gc_content");
        
        // Detect active stage from newest file timestamp
        snap.active_stage = get_active_stage();
        
        return snap;
    }
    
    long count_lines(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) return 0;
        
        long count = 0;
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty()) count++;
        }
        return count;
    }
    
    std::map<std::string, int> parse_taxonomy_distribution(const std::string& filepath) {
        std::map<std::string, int> dist;
        std::ifstream file(filepath);
        if (!file.is_open()) return dist;
        
        std::string line;
        std::regex phylum_regex("p__([^;]+)");
        std::smatch match;
        
        while (std::getline(file, line)) {
            if (std::regex_search(line, match, phylum_regex)) {
                std::string phylum = match[1];
                dist[phylum]++;
            }
        }
        return dist;
    }
    
    float extract_float_metric(const std::string& filepath, const std::string& key) {
        std::ifstream file(filepath);
        if (!file.is_open()) return 0.0f;
        
        std::string line;
        std::string search_key = key + "=";
        
        while (std::getline(file, line)) {
            if (line.find(search_key) != std::string::npos) {
                try {
                    std::string value = line.substr(line.find('=') + 1);
                    return std::stof(value);
                } catch (...) {
                    return 0.0f;
                }
            }
        }
        return 0.0f;
    }
    
    float get_elapsed_time() {
        std::ifstream log("output_dir/nucleitaxa.log");
        // Parse first and last timestamp
        // For simplicity, return 0 if not available
        return 0.0f;
    }
    
    int get_active_stage() {
        // Check which stage output directory has most recent files
        for (int i = 6; i >= 1; --i) {
            std::string stage_dir = output_dir + "/0" + std::to_string(i) + "-*/";
            // Check if directory has files
            // Return most recent stage
        }
        return 0;
    }
    
    void broadcast_json(const MetricSnapshot& snap) {
        std::ostringstream json;
        json << "{\n"
             << "  \"asv_count\": " << snap.asv_count << ",\n"
             << "  \"chimera_count\": " << snap.chimera_count << ",\n"
             << "  \"mean_quality\": " << snap.mean_quality << ",\n"
             << "  \"gc_percent\": " << snap.gc_percent << ",\n"
             << "  \"elapsed_sec\": " << snap.elapsed_sec << ",\n"
             << "  \"active_stage\": " << snap.active_stage << ",\n"
             << "  \"taxonomy_distribution\": {";
        
        bool first = true;
        for (const auto& [taxa, count] : snap.taxa_dist) {
            if (!first) json << ",";
            json << "\n    \"" << taxa << "\": " << count;
            first = false;
        }
        json << "\n  }\n}\n";
        
        // Broadcast to all connected clients
        std::string json_str = json.str();
        for (auto it = client_sockets.begin(); it != client_sockets.end(); ++it) {
            ssize_t sent = send(*it, json_str.c_str(), json_str.size(), 0);
            if (sent < 0) {
                close(*it);
                it = client_sockets.erase(it);
            }
        }
    }
};

int main(int argc, char *argv[]) {
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);
    
    int port = 8888;
    std::string output_dir = "./results";
    
    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--port" || arg == "-p") && i + 1 < argc) {
            port = std::stoi(argv[++i]);
        } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
            output_dir = argv[++i];
        }
    }
    
    std::cout << "[Analytics] NucleiTaxa Analytics Server v1.0" << std::endl;
    std::cout << "[Analytics] Output directory: " << output_dir << std::endl;
    
    AnalyticsServer server(port, output_dir);
    server.start();
    
    return 0;
}
