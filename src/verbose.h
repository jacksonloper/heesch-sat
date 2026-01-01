#pragma once

#include <iostream>
#include <chrono>
#include <string>
#include <iomanip>

// Global verbose flag - set from command line
inline bool g_verbose = false;

// Verbose log macro - only logs when g_verbose is true
#define VLOG(msg) do { if (g_verbose) { std::cerr << "[VERBOSE] " << msg << std::endl; } } while(0)
#define VLOG_NOENDL(msg) do { if (g_verbose) { std::cerr << "[VERBOSE] " << msg; } } while(0)

// Timer class for profiling
class VerboseTimer {
public:
    using clock = std::chrono::steady_clock;
    using time_point = clock::time_point;
    using duration = std::chrono::duration<double>;

    VerboseTimer(const std::string& name)
        : name_(name)
        , start_(clock::now())
    {}

    ~VerboseTimer() {
        if (g_verbose) {
            auto end = clock::now();
            double elapsed = duration(end - start_).count();
            std::cerr << "[TIMER] " << name_ << ": "
                      << std::fixed << std::setprecision(4) << elapsed << "s" << std::endl;
        }
    }

    double elapsed() const {
        return duration(clock::now() - start_).count();
    }

    void checkpoint(const std::string& label) {
        if (g_verbose) {
            double elapsed = duration(clock::now() - start_).count();
            std::cerr << "[TIMER] " << name_ << " @ " << label << ": "
                      << std::fixed << std::setprecision(4) << elapsed << "s" << std::endl;
        }
    }

private:
    std::string name_;
    time_point start_;
};

// RAII timer that logs on destruction
#define VTIMER(name) VerboseTimer _vtimer##__LINE__(name)

// Simple timer for manual timing
class ManualTimer {
public:
    using clock = std::chrono::steady_clock;
    using time_point = clock::time_point;
    using duration = std::chrono::duration<double>;

    ManualTimer() : start_(clock::now()) {}

    void reset() { start_ = clock::now(); }

    double elapsed() const {
        return duration(clock::now() - start_).count();
    }

private:
    time_point start_;
};
