#pragma once

#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

#include "types.h"

//            std::terminate(); \

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            throw std::runtime_error("Assertion failure"); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

namespace sweepmap {

/*
class TimerEmpty {
public:
    void start() {}
    void stop() {}
    double secs() const { return 0.0; }
    double range_ratio() const { return 0.0; }
};

class TimersEmpty {
public:
    void start(const std::string& name) {}
    void stop(const std::string& name) {}
    double secs(const std::string& name) const { return 0.0; }
    double range_ratio(const std::string& name) const { return 0.0; }
    double perc(const std::string& name, const std::string& total) const { return 0.0; }
    void clear() {}
    TimersEmpty& operator+=(const TimersEmpty& other) { return *this; }
};

typedef TimerEmpty Timer;
typedef TimersEmpty Timers;
*/

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
    
    void update_range(double diff) {
        if (diff < min_) min_ = diff;
        if (diff > max_) max_ = diff;
    }

    void stop() {
		if (!running_)
            throw std::runtime_error("Timer cannot be stopped since it is not running.");
        end_time_point_ = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration<double>(end_time_point_ - start_time_point_).count();
        accumulated_time_ += diff;
        running_ = false;
        update_range(diff);
    }

    double secs() const {
		if (running_) {
        //    throw std::runtime_error("Timer is still running.");
            return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time_point_).count();
        }
		return accumulated_time_;
    }

    double range_ratio() const {
		if (running_)
            throw std::runtime_error("Timer is still running.");

		if (max_ < 0.0)
            throw std::runtime_error("There is no range ratio since the timer has not been run.");
        return max_ / min_;
    }

    Timer& operator+=(const Timer& other) {
        accumulated_time_ += other.accumulated_time_;
        update_range(other.accumulated_time_);
        return *this;
    }
};

class Timers {
public:
    std::unordered_map<std::string, Timer> timers_;

    // TODO: change to use std::convertible_to<std::string> when C++20 is available
    template <typename... Args>
    void init(Args&&... s) {
        ((timers_[std::string(std::forward<Args>(s))] = Timer()), ...);
    }

    void start(const std::string& name) {
        //if (name == "mapping")
        //    std::cerr << "Starting timer " << name << std::endl;
        if (!timers_.contains(name))
            timers_[name] = Timer();
        timers_[name].start();
    }

    void stop(const std::string& name) {
        auto it = timers_.find(name);
        if (it != timers_.end()) {
            it->second.stop();
        }
    }

    double secs(const std::string& name) const {
		//ASSERT(timers_.contains(name), "Timer \"" << name << "\" not found.");
        try {
            return timers_.at(name).secs();
        } catch (const std::runtime_error& e) {
            std::cerr << "Error with timer \"" << name << "\": " << e.what() << std::endl;
            return 0.0;
        }
    }

    double range_ratio(const std::string& name) const {
        auto it = timers_.find(name);
		if (it == timers_.end())
            throw std::runtime_error("Timer " + name + " not found.");
        try {
            return it->second.range_ratio();
        } catch (const std::runtime_error& e) {
            std::cerr << "Error with timer \"" << name << "\": " << e.what() << std::endl;
            return 0.0;
        }
    }

    double perc(const std::string& name, const std::string& total) const {
        auto it = timers_.find(name);
		ASSERT(it != timers_.end(), "Timer \"" << name << "\" not found.");
        auto it_total = timers_.find(total);
		ASSERT(it_total != timers_.end(), "Timer \"" << total << "\" not found.");
        if (it != timers_.end()) {
			if (it_total != timers_.end()) {
				return it->second.secs() / it_total->second.secs() * 100.0;
			}
        }
        return 0.0;
    }

    void clear() {
        timers_.clear();
    }

    Timers& operator+=(const Timers& other) {
        for (const auto& [name, timer] : other.timers_) {
            timers_[name] += timer;
        }
        return *this;
    }
};

class Counter {
private:
    int64_t count_;

public:
    Counter() : count_(0) {}
    Counter(int64_t value) : count_(value) {}
    void inc(int64_t value = 1) { count_ += value; }
    int64_t count() const { return count_; }

    Counter& operator+=(const Counter& other) {
        count_ += other.count_;
        return *this;
    }
    
    Counter& operator=(const Counter& other) {
        count_ = other.count_;
        return *this;
    }

    bool operator<(const Counter& other) const {
        return count_ < other.count_;
    }
};

class Counters {
private:
    std::unordered_map<std::string, Counter> counters_;

public:
    template <typename... Args>
    void init(Args&&... s) {
        ((counters_[std::string(std::forward<Args>(s))] = Counter()), ...);
    }

    void inc(const std::string& name, int64_t value = 1) {
        //ASSERT(counters_.contains(name), "Counter \"" << name << "\" not found.");
        if (!counters_.contains(name))
            counters_[name] = Counter();
        counters_[name].inc(value);
    }

    int64_t count(const std::string& name) const {
        ASSERT(counters_.contains(name), "Counter \"" << name << "\" not found.");
        return counters_.at(name).count();
    }

    double frac(const std::string& name, const std::string& total) const {
        ASSERT(counters_.contains(name), "Counter \"" << name << "\" not found.");
        ASSERT(counters_.contains(total), "Counter \"" << total << "\" not found.");
        return double(counters_.at(name).count()) / double(counters_.at(total).count());
	}

    double perc(const std::string& name, const std::string& total) const {
		return 100.0 * frac(name, total);
	}

    void clear() {
        counters_.clear();
    }

    Counters& operator+=(const Counters& other) {
        for (const auto& [name, counter] : other.counters_) {
            counters_[name] += counter;
        }
        return *this;
    } 
    
    Counter& operator[](const std::string& name) {
        return counters_[name];
    }
};

class ProgressBar {
private:
    static constexpr const char* PBSTR = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
    static constexpr int PBWIDTH = 60;
    std::string message_;
    std::ostream& out_;

public:
    explicit ProgressBar(std::string message = "Progress", std::ostream& out = std::cerr) 
        : message_(std::move(message)), out_(out) {}

    void update(double progress) const {
        int val = static_cast<int>(progress * 100);
        int lpad = static_cast<int>(progress * PBWIDTH);
        int rpad = PBWIDTH - lpad;
        out_ << "\r" << message_ << " " << std::setw(3) << val << "% ["
             << std::string(lpad, '|') << std::string(rpad, ' ') << "]" 
             << std::flush;
    }
};

struct ParsedQueryId {
	bool valid;
	std::string segm_id;
	rpos_t start_pos;
	rpos_t end_pos;
	char strand;
	
	static ParsedQueryId parse(const std::string& query_id) {
		// query_id is in the form "S1_21!NC_060948.1!57693539!57715501!+"
		ParsedQueryId result{false, "", 0, 0, '?'};
		
		// Skip if query_id doesn't contain correct mapping info
		if (query_id.find('!') == std::string::npos) return result;

		// Split query_id on '!' character
		std::vector<std::string> parts;
		size_t start = 0;
		size_t end = query_id.find('!');
		while (end != std::string::npos) {
			parts.push_back(query_id.substr(start, end - start));
			start = end + 1;
			end = query_id.find('!', start);
		}
		parts.push_back(query_id.substr(start));

		// Need 5 parts: read_id, segment, start, end, strand
		if (parts.size() != 5) return result;

		try {
			result.segm_id = parts[1];
			result.start_pos = stoll(parts[2]);
			result.end_pos = stoll(parts[3]);
			result.strand = parts[4][0];
			result.valid = true;
		} catch (const std::exception& e) {
			std::cerr << "Error parsing query_id with start_pos " << parts[2] << " and end_pos " << parts[3] << ": " << e.what() << std::endl;
			result.valid = false;
			throw;
		}
		
		return result;
	}
};

} // namespace sweepmap