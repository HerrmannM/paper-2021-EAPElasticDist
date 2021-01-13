#pragma once

#include "../distances.hpp"

namespace internal {

    /** Weighted Dynamic Time Warping with cutoff point for early abandoning and pruning.
     *  Double buffered implementation using O(n) space.
     *  Worst case scenario has a O(n²) time complexity (no pruning nor early abandoning).
     *  A tight cutoff can allow a lot of pruning, speeding up the process considerably.
     *  Actual implementation assuming that some pre-conditions are fulfilled.
     * @param lines     Pointer to the "line series". Must be the longest series. Cannot be null.
     * @param nblines   Length of the line series. Must be 0 < nbcols <= nblines < tempo::MAX_SERIES_LENGTH.
     * @param cols      Pointer to the "column series". Must be the shortest series. Cannot be null.
     * @param nbcols    Length of the column series. Must be 0 < nbcols <= nblines < tempo::MAX_SERIES_LENGTH.
     * @param weights   Pointer to the weights. Must be at least as long as nblines.
     * @param cutoff.   Attempt to prune computation of alignments with cost > cutoff.
     *                  May lead to early abandoning.
     * @return WTW between the two series or +INF if early abandoned.
     */
    template<bool doLA=false>
    [[nodiscard]] inline double wdtw(
            const double *lines, size_t nblines,
            const double *cols, size_t nbcols,
            const double *weights,
            double cutoff
    ) {
        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // In debug mode, check preconditions
        assert(lines != nullptr && nblines != 0 && nblines < MAX_SERIES_LENGTH);
        assert(cols != nullptr && nbcols != 0 && nbcols < MAX_SERIES_LENGTH);
        assert(nbcols <= nblines);

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Create a new tighter upper bounds (most commonly used in the code).
        // First, take the "next float" after "cutoff" to deal with numerical instability.
        // Then, subtract the cost of the last alignment.
        double ub = initBlock{
            const auto ll = nblines-1;
            const auto lc = nbcols-1;  // Precondition: ll>=lc, so ll-lc>=0, well defined for unsigned size_t.
            return nextafter(cutoff, POSITIVE_INFINITY) -square_dist(lines[ll], cols[lc])*weights[ll-lc];
        };

        if constexpr (!doLA){ ub = cutoff; }

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Double buffer allocation, init to +INF.
        // Base indices for the 'c'urrent row and the 'p'revious row. Account for the extra cell (+1 and +2)
        std::vector<double> buffers_v((1 + nbcols) * 2, POSITIVE_INFINITY);
        auto* buffers = buffers_v.data();
        size_t c{0+1}, p{nbcols+2};

        // Line & column counters
        size_t i{0}, j{0};

        // Cost accumulator. Also used as the "left neighbour".
        double cost{0};

        // EAP variables: track where to start the next line, and the position of the previous pruning point.
        // Must be init to 0: index 0 is the next starting index and also the "previous pruning point"
        size_t next_start{0}, prev_pp{0};

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Initialisation of the top border: already initialized to +INF. Initialise the left corner to 0.
        buffers[c-1] = 0;

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Main loop
        for (; i < nblines; ++i) {
            // --- --- --- Swap and variables init
            std::swap(c, p);
            const double li = lines[i];
            size_t curr_pp = next_start; // Next pruning point init at the start of the line
            j = next_start;
            // --- --- --- Stage 0: Initialise the left border
            {
                cost = POSITIVE_INFINITY;
                buffers[c+next_start-1] = POSITIVE_INFINITY;
            }
            // --- --- --- Stage 1: Up to the previous pruning point while advancing next_start: diag and top
            for (; j == next_start && j < prev_pp; ++j) {
                const auto d = square_dist(li, cols[j])*weights[absdiff(i,j)];
                cost = std::min(buffers[p + j - 1], buffers[p + j]) + d;
                buffers[c + j] = cost;
                if (cost <= ub) { curr_pp = j + 1; } else { ++next_start; }
            }
            // --- --- --- Stage 2: Up to the previous pruning point without advancing next_start: left, diag and top
            for (; j < prev_pp; ++j) {
                const auto d = square_dist(li, cols[j])*weights[absdiff(i,j)];
                cost = min(cost, buffers[p + j - 1], buffers[p + j]) + d;
                buffers[c + j] = cost;
                if (cost <= ub) { curr_pp = j + 1; }
            }
            // --- --- --- Stage 3: At the previous pruning point. Check if we are within bounds.
            if (j < nbcols) { // If so, two cases.
                const auto d = square_dist(li, cols[j])*weights[absdiff(i,j)];
                if (j == next_start) { // Case 1: Advancing next start: only diag.
                    cost = buffers[p + j - 1] + d;
                    buffers[c + j] = cost;
                    if (cost <= ub) { curr_pp = j + 1; }
                    else {
                        // Special case if we are on the last alignment: return the actual cost if we are <= cutoff
                        if (i == nblines - 1 && j == nbcols - 1 && cost <= cutoff) { return cost; }
                        else { return POSITIVE_INFINITY; }
                    }
                } else { // Case 2: Not advancing next start: possible path in previous cells: left and diag.
                    cost = std::min(cost, buffers[p + j - 1]) + d;
                    buffers[c + j] = cost;
                    if (cost <= ub) { curr_pp = j + 1; }
                }
                ++j;
            } else { // Previous pruning point is out of bound: exit if we extended next start up to here.
                if (j == next_start) {
                    // But only if we are above the original UB
                    // Else set the next starting point to the last valid column
                    if (cost > cutoff) { return POSITIVE_INFINITY; }
                    else { next_start = nbcols - 1; }
                }
            }
            // --- --- --- Stage 4: After the previous pruning point: only prev.
            // Go on while we advance the curr_pp; if it did not advance, the rest of the line is guaranteed to be > ub.
            for (; j == curr_pp && j < nbcols; ++j) {
                const auto d = square_dist(li, cols[j])*weights[absdiff(i,j)];
                cost = cost + d;
                buffers[c + j] = cost;
                if (cost <= ub) { ++curr_pp; }
            }
            // --- --- ---
            prev_pp = curr_pp;
        } // End of main loop for(;i<nblines;++i)

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Finalisation
        // Check for last alignment (i==nblines implied, Stage 4 implies j<=nbcols). Cost must be <= original bound.
        if (j == nbcols && cost <= cutoff) { return cost; }
        else { return POSITIVE_INFINITY; }
    }

} // End of namespace internal

// --- --- --- --- ---
// Weights generation
// --- --- --- --- ---

/// From the paper, changing this values does not change the results (scaling), so keep to 1
constexpr double WDTW_MAX_WEIGHT = 1;

/** Compute a weight at index i in a sequence 1..m
 * @param g "Controls the level of penalization for the points with larger phase difference".
 *        range [0, +inf), usually in [0.01, 0.6].
 *        Some examples:
 *        * 0: constant weight
 *        * 0.05: nearly linear weights
 *        * 0.25: sigmoid weights
 *        * 3: two square_distinct weights between half sequences
 *
 * @param half_max_length Mid point of the sequence (m/2)
 * @param i Index of the point in [1..m] (m=length of the sequence)
 * @return the weight for index i
 */
[[nodiscard]] inline double compute_weight(double g, double half_max_length, double i) {
    return WDTW_MAX_WEIGHT / (1 + exp(-g * (i - half_max_length)));
}

/// Populate the weights_array of size length with weights derive from the g factor
inline void populate_weights(double g, double *weights_array, size_t length) {
    double half_max_length = double(length) / 2;
    for (size_t i{0}; i < length; ++i) {
        weights_array[i] = compute_weight(g, half_max_length, double(i));
    }
}

/// Create a vector of weights
inline std::vector<double> generate_weights(double g, size_t length){
    std::vector<double> weights(length, 0);
    populate_weights(g, weights.data(), length);
    return weights;
}


template<bool doLA=false>
[[nodiscard]] double wdtw(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        const double *weights
) {
    const auto check_result = check_order_series(series1, length1, series2, length2);
    switch (check_result.index()) {
        case 0: { return std::get<0>(check_result);}
        case 1: {
            const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Compute a cutoff point using the diagonal
            double cutoff{0};
            // We have less columns than lines: cover all the columns first.
            for (size_t i{0}; i < nbcols; ++i) { cutoff += square_dist(lines[i], cols[i])*weights[0]; }
            // Then go down in the last column
            if(nbcols<nblines) {
                const auto lc = cols[nbcols - 1];
                for (size_t i{nbcols}; i < nblines; ++i) { cutoff += square_dist(lines[i], lc)*weights[i-nbcols+1]; }
            }

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            return internal::wdtw<doLA>(lines, nblines, cols, nbcols, weights, cutoff);
        }
        default: should_not_happen();
    }
}

template<bool doLA=false>
[[nodiscard]] double wdtw(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        const double *weights,
        double cutoff
) {
    const auto check_result = check_order_series(series1, length1, series2, length2);
    switch (check_result.index()) {
        case 0: { return std::get<0>(check_result);}
        case 1: {
            const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);
            return internal::wdtw<doLA>(lines, nblines, cols, nbcols, weights, cutoff);
        }
        default: should_not_happen();
    }
}






/** Weigthed DTW, implemented on a double buffer
 * @param series1 Pointer to the first series' values
 * @param length1 Length of the first series
 * @param series2 Pointer to the second series' values
 * @param length2 Length of the second series
 * @param weights Pointer to the weights. Underlying array must be at least as long as the longest series.
 * @return WDTW value
 */
[[nodiscard]] double wdtw_base(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        const double *weights
) {
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Pre-conditions. Accept nullptr if length is 0
    assert((series1 != nullptr || length1==0) && length1 < MAX_SERIES_LENGTH);
    assert((series2 != nullptr || length2==0) && length2 < MAX_SERIES_LENGTH);
    // Check sizes. If both series are empty, return 0, else if one is empty and not the other, maximal error.
    if (length1 == 0 && length2 == 0) { return 0; }
    else if ((length1 == 0) != (length2 == 0)) { return POSITIVE_INFINITY; }
    // Use the smallest size as the columns (which will be the allocation size)
    const double *cols = (length1 < length2) ? series1 : series2;
    const double *lines = (length1 < length2) ? series2 : series1;
    const size_t nbcols = std::min(length1, length2);
    const size_t nblines = std::max(length1, length2);
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

    // Double buffer allocation, no initialisation required (border condition manage in the code).
    // Base indices for the 'c'urrent row and the 'p'revious row.
    auto buffers = std::unique_ptr<double[]>(new double[nbcols*2]);
    std::vector<double> buffers_v(nbcols * 2, 0);
    size_t c{0}, p{nbcols};

    // Line counter
    size_t i{0};
    // Variable storing the cost (also act as the cost of the left neighbour)
    double cost;

    // --- Initialization of the first line, take the border condition into account
    {
        const double li = lines[i];
        // First cell is a special case
        {
            cost = square_dist(li, cols[0])*weights[0];
            buffers[c + 0] = cost;
        }
        // Rest of the line, a cell only depends on the previous cell
        for (size_t j{1}; j < nbcols; ++j) {
            cost = cost + square_dist(li, cols[j])*weights[j];
            buffers[c + j] = cost;
        }
        // Finalisation
        ++i;
    }

    // --- Main loop
    for (; i < nblines; ++i) {
        // --- --- --- Swap and variables init
        std::swap(c, p);
        const double li = lines[i];
        // --- --- --- Handle border condition: compute the first column which can only depends on the "top" item
        {
            cost = buffers[p + 0] + square_dist(li, cols[0])*weights[i];
            buffers[c + 0] = cost;
        }
        // --- --- --- Iterate through the columns
        for (size_t j{1}; j < nbcols; ++j) {
            cost = min(cost, buffers[p + j - 1], buffers[p + j]) + square_dist(li, cols[j]) * weights[absdiff(i, j)];
            buffers[c + j] = cost;
        }
    }

    // --- Finalisation
    return buffers[c + nbcols - 1];
}

[[nodiscard]] double wdtw_base_ea(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        const double *weights,
        double cutoff
) {
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Pre-conditions. Accept nullptr if length is 0
    assert((series1 != nullptr || length1 == 0) && length1 < MAX_SERIES_LENGTH);
    assert((series2 != nullptr || length2 == 0) && length2 < MAX_SERIES_LENGTH);
    // Check sizes. If both series are empty, return 0, else if one is empty and not the other, maximal error.
    if (length1 == 0 && length2 == 0) { return 0; }
    else if ((length1 == 0) != (length2 == 0)) { return POSITIVE_INFINITY; }
    // Use the smallest size as the columns (which will be the allocation size)
    const double *cols = (length1 < length2) ? series1 : series2;
    const double *lines = (length1 < length2) ? series2 : series1;
    const size_t nbcols = std::min(length1, length2);
    const size_t nblines = std::max(length1, length2);
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

    // Double buffer allocation, no initialisation required (border condition manage in the code).
    // Base indices for the 'c'urrent row and the 'p'revious row.
    auto buffers = std::unique_ptr<double[]>(new double[nbcols * 2]);
    std::vector<double> buffers_v(nbcols * 2, 0);
    size_t c{0}, p{nbcols};

    // Line counter
    size_t i{0};
    // Variable storing the cost (also act as the cost of the left neighbour)
    double cost;
    double minv = POSITIVE_INFINITY;

    // --- Initialization of the first line, take the border condition into account
    {
        const double li = lines[i];
        // First cell is a special case
        {
            cost = square_dist(li, cols[0]) * weights[0];
            buffers[c + 0] = cost;
        }
        // Rest of the line, a cell only depends on the previous cell
        for (size_t j{1}; j < nbcols; ++j) {
            cost = cost + square_dist(li, cols[j]) * weights[j];
            buffers[c + j] = cost;
        }
        // Finalisation
        ++i;
    }

    // --- Main loop
    for (; i < nblines; ++i) {
        // --- --- --- Swap and variables init
        std::swap(c, p);
        const double li = lines[i];
        // --- --- --- Handle border condition: compute the first column which can only depends on the "top" item
        {
            cost = buffers[p + 0] + square_dist(li, cols[0]) * weights[i];
            buffers[c + 0] = cost;
            minv = cost;
        }
        // --- --- --- Iterate through the columns
        for (size_t j{1}; j < nbcols; ++j) {
            cost = min(cost, buffers[p + j - 1], buffers[p + j]) + square_dist(li, cols[j]) * weights[absdiff(i, j)];
            buffers[c + j] = cost;
            minv = std::min(minv, cost);
        }
        if (minv > cutoff) { return POSITIVE_INFINITY; }
    }

    // --- Finalisation
    return buffers[c + nbcols - 1];
}
