#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <chrono>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <variant>
#include <vector>

#include <iostream>
#include <ostream>
#include <sstream>



// --- --- --- --- --- ---
// --- Series constants
// --- --- --- --- --- ---

/// Constant to be use when no window is required
constexpr size_t NO_WINDOW{std::numeric_limits<size_t>::max()};

/// Constant representing the maximum length allowed for a series.
/// Account for extra columns/lines that may need to be allocated and representing the window.
constexpr size_t MAX_SERIES_LENGTH{NO_WINDOW - 2};

// --- --- --- --- --- ---
// --- Floating point constants
// --- --- --- --- --- ---

/// Positive infinity
constexpr double POSITIVE_INFINITY{std::numeric_limits<double>::infinity()};


// --- --- --- --- --- ---
// --- Tooling
// --- --- --- --- --- ---

/// Minimum of 3 values using std::min<T>
template<typename T>
[[nodiscard]] inline T min(T a, T b, T c) { return std::min<T>(a, std::min<T>(b, c)); }

/// Maximum of 3 values using std::min<T>
template<typename T>
[[nodiscard]] inline T max(T a, T b, T c) { return std::max<T>(a, std::max<T>(b, c)); }


// --- --- --- --- --- ---
// --- Should not happen
// --- --- --- --- --- ---

/// Throw an exception "should not happen". Used as default case in switches.
[[noreturn]] void inline should_not_happen(){ throw std::logic_error("Should not happen"); }

// --- --- --- --- --- ---
// --- Initialisation tool
// --- --- --- --- --- ---

namespace initBlock_detail {
    struct tag { };

    template <class F>
    decltype(auto) operator + (tag, F &&f) {
        return std::forward<F>(f)();
    }
}

#define initBlock initBlock_detail::tag{} + [&]() -> decltype(auto)

// --- --- ---
// --- --- --- Basic set and map "contains" function
// --- --- ---

/** sset "contains" function */
template<typename Key, class Compare = std::less<Key>, class Allocator = std::allocator<Key>>
[[nodiscard]] inline bool contains(const std::set<Key, Compare, Allocator> &s, const Key &key) {
    return s.find(key) != s.end();
}

/** smap "contains" function */
template<typename Key, typename T, typename Compare = std::less<Key>, typename Allocator = std::allocator<std::pair<const Key, T> > >
[[nodiscard]] inline bool contains(const std::map<Key, T, Compare, Allocator> &m, const Key &key) {
    return m.find(key) != m.end();
}


// --- --- ---
// --- --- --- Uncopyable
// --- --- ---

/** Private inherit from this class to create an uncopyable (but still movable) class */
class Uncopyable {
protected:
    // Protect the constructor and the destructor (class not usable by itself)
    Uncopyable() = default;

    ~Uncopyable() = default;

    // Still movable
    Uncopyable(Uncopyable &&) = default;

    Uncopyable &operator=(Uncopyable &&) = default;

    // Delete copy and copy-assignment operator
    Uncopyable(const Uncopyable &other) = delete;

    Uncopyable &operator=(const Uncopyable &other) = delete;
};


// --- --- ---
// --- --- --- Clock related types aliases
// --- --- ---

using myclock_t = std::chrono::steady_clock;
using duration_t = myclock_t::duration;
using time_point_t = myclock_t::time_point;

// --- --- ---
// --- --- --- Functions
// --- --- ---

/** Create a time point for "now" */
inline time_point_t now() { return myclock_t::now(); }

/** Print a duration in a human readable form (from nanoseconds to hours) in an output stream. */
inline void printDuration(std::ostream &out, const duration_t &elapsed) {
    namespace c = std::chrono;
    auto execution_time_ns = c::duration_cast<c::nanoseconds>(elapsed).count();
    auto execution_time_us = c::duration_cast<c::microseconds>(elapsed).count();
    auto execution_time_ms = c::duration_cast<c::milliseconds>(elapsed).count();
    auto execution_time_sec = c::duration_cast<c::seconds>(elapsed).count();
    auto execution_time_min = c::duration_cast<c::minutes>(elapsed).count();
    auto execution_time_hour = c::duration_cast<c::hours>(elapsed).count();

    bool first = true;

    if (execution_time_hour > 0) {
        first = false; // no need to test, if above condition is true, this is the first
        out << execution_time_hour << "h";
    }
    if (execution_time_min > 0) {
        if (first) { first = false; } else { out << " "; }
        out << execution_time_min % 60 << "m";
    }
    if (execution_time_sec > 0) {
        if (first) { first = false; } else { out << " "; }
        out << "" << execution_time_sec % 60 << "s";
    }
    if (execution_time_ms > 0) {
        if (first) { first = false; } else { out << " "; }
        out << "" << execution_time_ms % long(1E+3) << "ms";
    }
    if (execution_time_us > 0) {
        if (first) { first = false; } else { out << " "; }
        out << "" << execution_time_us % long(1E+3) << "us";
    }
    if (execution_time_ns >= 0) {
        if (first) { first = false; } else { out << " "; }
        out << "" << execution_time_ns % long(1E+3) << "ns";
    }
}


/** Shortcut for the above function, converting two time points into a duration. */
inline void printExecutionTime(std::ostream &out, time_point_t start_time, time_point_t end_time) {
    const auto elapsed = end_time - start_time;
    printDuration(out, elapsed);
}

/** Shortcut to print in a string */
[[nodiscard]] inline std::string as_string(time_point_t start_time, time_point_t end_time) {
    std::stringstream ss;
    printExecutionTime(ss, start_time, end_time);
    return ss.str();
}


// --- --- ---
// --- --- --- istream
// --- --- ---

/** Test if a char is a whitespace (excluding a newline). */
static inline bool is_white(int c) {
    static auto w = std::string(" \f\r\t\v");
    return w.find(c) != std::string::npos;
}

/** Skip while whitespaces (excluding newline) are read. */
static inline void skip_white(std::istream &input) {
    int c{EOF};
    while ((c = input.peek()) != EOF && is_white(c)) { input.ignore(1); }
}

/** Skip the current line (read until newline) */
static inline void skip_line(std::istream &input) {
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

/** Read a word in the buffer (read until a whitespace, including newline, is found).
 *  Return the char that ended the read, keeping it in the stream. */
static inline int read_word(std::istream &input, std::string &buffer) {
    int c = EOF;
    while ((c = input.peek()) != EOF && !std::isspace(c)) { buffer.push_back(input.get()); }
    return c;
}

/** Read a word in the buffer (read until a whitespace, including newline, is found), converting to lower case.
 *  Return the char that ended the read, keeping it in the stream. */
static inline int read_word_to_lower(std::istream &input, std::string &buffer) {
    int c{EOF};
    while ((c = input.peek()) != EOF && !std::isspace(c)) { buffer.push_back(std::tolower(input.get())); }
    return c;
}

/** Test if a char is a delimiter ',' or ':' or '\n' */
static inline bool is_delim(int c) {
    return c == ',' || c == ':' || c == '\n';
}

/** Read into buffer until a delimiter is found; returns that delimiter (taken out opf the stream). */
static inline int read_until_delim(std::istream &input, std::string &buffer) {
    int c{EOF};
    while ((c = input.peek()) != EOF && !is_delim(c)) { buffer.push_back(input.get()); }
    if (c != EOF) { input.ignore(1); }
    return c;
}

// --- --- ---
// --- --- --- string/istringstream
// --- --- ---

/** String splitter on a delimiter. Accept a istringstream */
static inline std::vector<std::string> split(std::istringstream &&input, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    while (std::getline(input, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

/** trim from start (in place) */
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

/** trim from end (in place) */
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

/** trim from both ends (in place) */
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

/** Attempt to convert a string into a bool */
static inline std::optional<bool> as_bool(const std::string &str) {
    if (str == "true") {
        return {true};
    } else if (str == "false") {
        return {false};
    } else {
        return {};
    }
}

/** Attempt to convert a string into an integer */
static inline std::optional<int> as_int(const std::string &str) {
    try {
        int i = std::stoi(str);
        return {i};
    } catch (...) {
        return {};
    }
}

/** Attempt to convert a string into an size_t */
static inline std::optional<size_t> as_size_t(const std::string &str) {
    try {
        size_t i = std::stoul(str);
        return {i};
    } catch (...) {
        return {};
    }
}

/** Attempt to convert a string into an double */
static inline std::optional<double> as_double(const std::string &str) {
    try {
        double d = std::stod(str);
        return {d};
    } catch (...) {
        return {};
    }
}

