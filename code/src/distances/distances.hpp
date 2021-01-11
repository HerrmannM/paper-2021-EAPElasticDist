#pragma once

#include <functional>
#include <tuple>
#include <variant>

#include <cmath>
#include "../utils.hpp"

/** Unsigned arithmetic:
 * Given an 'index' and a 'window', get the start index corresponding to std::max(0, index-window) */
[[nodiscard]] inline size_t cap_start_index_to_window(size_t index, size_t window) {
    if (index > window) { return index - window; } else { return 0; }
}

/** Unsigned arithmetic:
 * Given an 'index', a 'window' and an 'end', get the stop index corresponding to std::min(end, index+window+1).
 * The expression index+window+1 is illegal for any index>0 as window could be MAX-1
 * */
[[nodiscard]] inline size_t
cap_stop_index_to_window_or_end(size_t index, size_t window, size_t end) {
    // end-window is valid when window<end
    if (window < end && index + 1 < end - window) { return index + window + 1; } else { return end; }
}

/** Absolute value for any comparable and subtractive type, without overflowing risk for unsigned types.
 *  Also work for signed type. */
template<typename T>
[[nodiscard]] inline T absdiff(T a, T b) { return (a > b) ? a - b : b - a; }

/// Type alias for a tuple representing (lines, nblines, cols, nbcosl)
using lico_t = std::tuple<const double*, size_t, const double*, size_t>;

/// Helper function checking and ordering the length of the series
[[nodiscard]] inline std::variant<double, lico_t> check_order_series(
        const double *series1, size_t length1,
        const double *series2, size_t length2
){
    // Pre-conditions. Accept nullptr if length is 0
    assert((series1 != nullptr || length1 == 0) && length1 < MAX_SERIES_LENGTH);
    assert((series2 != nullptr || length2 == 0) && length2 < MAX_SERIES_LENGTH);
    // Check sizes. If both series are empty, return 0, else if one is empty and not the other, maximal error.
    if (length1 == 0 && length2 == 0) { return {double(0.0)}; }
    else if ((length1 == 0) != (length2 == 0)) { return POSITIVE_INFINITY; }
    // Use the smallest size as the columns (which will be the allocation size)
    return (length1 > length2) ?
           std::tuple(series1, length1, series2, length2) :
           std::tuple(series2, length2, series1, length1);
}

/// Square square_distances between two numeric values
template<typename T>
[[nodiscard]] inline T square_dist(T a, T b){
    const auto d = a-b;
    return d*d;
}

