#pragma once

#include "../distances.hpp"

// --- --- --- --- --- ---
// Element Wise
// --- --- --- --- --- ---

/** Element wise square_distance. Default to squared euclidean square_distance.
 * Only defined for same length series (return +INF if different length).
 * @param series1   First series
 * @param length1   Length of the first series
 * @param series2   Second series
 * @param length2   Length of the second series
 * @return Sum of element wise square_distances or +INF if different lengths
 */
[[nodiscard]] inline double elementwise(
        const double *series1, size_t length1,
        const double *series2, size_t length2
) {
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Pre-conditions. Accept nullptr if length is 0
    assert((series1 != nullptr || length1 == 0) && length1 < MAX_SERIES_LENGTH);
    assert((series2 != nullptr || length2 == 0) && length2 < MAX_SERIES_LENGTH);
    // Check sizes. If both series are empty, return 0, else if one is empty and not the other, maximal error.
    if (length1 != length2) { return POSITIVE_INFINITY; }
    // Compute the Euclidean-like square_distance
    double cost = 0.0;
    for (size_t i{0}; i < length1; ++i) { cost += square_dist(series1[i], series2[i]); }
    return cost;
}


// --- --- --- --- --- ---
// Element Wise with cut-off
// --- --- --- --- --- ---

/** Element wise square_distance with cut-off point for early abandoning.
 * Default to squared euclidean square_distance.
 * Only defined for same length series (return +INF if different length).
 * @param series1   First series
 * @param length1   Length of the first series
 * @param series2   Second series
 * @param length2   Length of the second series
 * @return Sum of element wise square_distances or +INF if different lengths or early abandoned
 */
template<bool doLA=false>
[[nodiscard]] inline double elementwise(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        double cutoff
) {
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Pre-conditions. Accept nullptr if length is 0
    assert((series1 != nullptr || length1 == 0) && length1 < MAX_SERIES_LENGTH);
    assert((series2 != nullptr || length2 == 0) && length2 < MAX_SERIES_LENGTH);
    // Check sizes. If both series are empty, return 0, else if one is empty and not the other, maximal error.
    if (length1 != length2) { return POSITIVE_INFINITY; }

    if constexpr (doLA)
    {

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Create a new tighter upper bounds (most commonly used in the code).
        // First, take the "next float" after "cutoff" to deal with numerical instability.
        // Then, subtract the cost of the last alignment.
        // Adjust the lower bound, taking the last alignment into account
        const double lastA = square_dist(series1[length1 - 1], series2[length1 - 1]);
        const double ub = std::nextafter(cutoff, POSITIVE_INFINITY) - lastA;

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Compute the Euclidean-like square_distance up to, excluding, the last alignment
        double cost = 0;
        for (size_t i{0}; i < length1 - 1; ++i) { // Stop before the last: counted in the bound!
            cost += square_dist(series1[i], series2[i]);
            if (cost > ub) { return POSITIVE_INFINITY; }
        }
        // Add the last alignment and check the result
        cost += lastA;
        if (cost > cutoff) { return POSITIVE_INFINITY; } else { return cost; }

    } else {

        double cost = 0;
        for (size_t i{0}; i < length1; ++i) {
            cost += square_dist(series1[i], series2[i]);
            if (cost > cutoff) { return POSITIVE_INFINITY; }
        }
        return cost;
    }
}
