#pragma once

#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

namespace sweepmap {

using hash_t     = uint64_t;
using rpos_t     = int32_t;  // reference position; klib anyway doesn't support 64-bit
using qpos_t     = int32_t;  // query position
using segm_t     = int32_t;

class Timer {
public:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_point_, end_time_point_;
    double accumulated_time_;
    double min_, max_;
    bool running_;

    Timer() : accumulated_time_(0.0), min_(1e9), max_(-1.0), running_(false) {}

    void start() {
		if (running_)
            throw std::runtime_error("Timer cannot be started since it is already running.");
        start_time_point_ = std::chrono::high_resolution_clock::now();
        running_ = true;
    }

    void stop() {
		if (!running_)
            throw std::runtime_error("Timer cannot be stopped since it is not running.");
        end_time_point_ = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration<double>(end_time_point_ - start_time_point_).count();
        accumulated_time_ += diff;
        running_ = false;
        if (diff < min_) min_ = diff;
        if (diff > max_) max_ = diff;
    }

    double secs() const {
		if (running_)
            throw std::runtime_error("Timer is still running.");
		return accumulated_time_;
    }

    double range_ratio() const {
		if (running_)
            throw std::runtime_error("Timer is still running.");
		if (max_ < 0.0)
            throw std::runtime_error("There is no range ratio since the timer has not been run.");
        return max_ / min_;
    }
};

class Timers {
public:
    std::unordered_map<std::string, Timer> timers_;

    void start(const std::string& name) {
        //std::cerr << "Starting timer " << name << std::endl;
        timers_[name].start();
    }

    void stop(const std::string& name) {
        auto it = timers_.find(name);
        if (it != timers_.end()) {
            it->second.stop();
        }
    }

    double secs(const std::string& name) const {
        auto it = timers_.find(name);
		assert(it != timers_.end());
        if (it != timers_.end()) {
            return it->second.secs();
        }
        assert(false);
        return 0.0;
    }

    double range_ratio(const std::string& name) const {
        auto it = timers_.find(name);
		if (it == timers_.end())
            throw std::runtime_error("Timer " + name + " not found.");
        return it->second.range_ratio();
    }

    double perc(const std::string& name, const std::string& total) const {
        auto it = timers_.find(name);
		assert(it != timers_.end());
        auto it_total = timers_.find(total);
		assert(it_total != timers_.end());
        if (it != timers_.end()) {
			if (it_total != timers_.end()) {
				return it->second.secs() / it_total->second.secs() * 100.0;
			}
        }
        return 0.0;
    }
};

class Counter {
private:
    int64_t count_;

public:
    Counter() : count_(0) {}
    void inc(int64_t value = 1) { count_ += value; }
    int64_t count() const { return count_; }
    
    // Non-const operator[] that returns reference to count
    int64_t& operator[](int) { return count_; }

    Counter& operator=(int64_t value) {
        count_ = value;
        return *this;
    }

    Counter& operator+=(int64_t value) {
        count_ += value;
        return *this;
    }

    Counter& operator+=(const Counter& other) {
        count_ += other.count_;
        return *this;
    }

    // Implicit conversion operator to int64_t
    operator int64_t() const { return count_; }
};

class Counters {
private:
    std::unordered_map<std::string, Counter> counters_;

public:
    void inc(const std::string& name, int64_t value = 1) {
		auto it = counters_.find(name);
        if (it == counters_.end()) {
            counters_[name] = Counter();
			it = counters_.find(name);
        }
        counters_[name].inc(value);
    }

    int64_t count(const std::string& name) const {
        assert(counters_.find(name) != counters_.end());
        return counters_.at(name).count();
    }

    double frac(const std::string& name, const std::string& total) const {
        assert(counters_.find(name) != counters_.end());
        assert(counters_.find(total) != counters_.end());
        return double(counters_.at(name).count()) / double(counters_.at(total).count());
	}

    double perc(const std::string& name, const std::string& total) const {
		return 100.0 * frac(name, total);
	}

    // Non-const operator[] that returns reference to Counter
    Counter& operator[](const std::string& name) {
        auto it = counters_.find(name);
        if (it == counters_.end()) {
            counters_[name] = Counter();
        }
        return counters_[name];
    }

    //Const operator[] that returns const int64_t
    int64_t operator[](const std::string& name) const {
        auto it = counters_.find(name);
        if (it == counters_.end()) {
            return 0;
        }
        return it->second.count();
    }
};

} // namespace sweepmap
