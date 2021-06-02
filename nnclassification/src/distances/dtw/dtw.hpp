#pragma once

#include "../distances.hpp"

namespace internal {

    /** Dynamic Time Warping with cutoff point for early abandoning and pruning.
     *  Double buffered implementation using O(n) space.
     *  Worst case scenario has a O(n²) time complexity (no pruning nor early abandoning).
     *  A tight cutoff can allow a lot of pruning, speeding up the process considerably.
     *  Actual implementation assuming that some pre-conditions are fulfilled.
     * @param lines     Pointer to the "line series". Must be the longest series. Cannot be null.
     * @param nblines   Length of the line series. Must be 0 < nbcols <= nblines < tempo::MAX_SERIES_LENGTH.
     * @param cols      Pointer to the "column series". Must be the shortest series. Cannot be null.
     * @param nbcols    Length of the column series. Must be 0 < nbcols <= nblines < tempo::MAX_SERIES_LENGTH.
     * @param cutoff.   Attempt to prune computation of alignments with cost > cutoff.
     *                  May lead to early abandoning.
     * @return DTW between the two series or +INF if early abandoned.
     */
    template<bool doLA=false>
    [[nodiscard]] inline double dtw(
            const double *lines, size_t nblines,
            const double *cols, size_t nbcols,
            const double cutoff
    ) {
        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // In debug mode, check preconditions
        assert(lines != nullptr && nblines != 0 && nblines < MAX_SERIES_LENGTH);
        assert(cols != nullptr && nbcols != 0 && nbcols < MAX_SERIES_LENGTH);
        assert(nbcols <= nblines);

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Create a new tighter upper bounds (most commonly used in the code).
        // First, take the "next float" after "cutoff" to deal with numerical instability.
        // Then, subtract the cost of the last alignment.
        double ub = nextafter(cutoff, POSITIVE_INFINITY) - square_dist(lines[nblines - 1], cols[nbcols - 1]);
        if constexpr (!doLA){ ub = cutoff; }

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Double buffer allocation, no initialisation required (border condition manage in the code).
        // Base indices for the 'c'urrent row and the 'p'revious row.
        auto buffers = std::unique_ptr<double[]>(new double[nbcols * 2]);
        size_t c{0}, p{nbcols};

        // Line & column counters
        size_t i{0}, j{0};

        // Cost accumulator. Also used as the "left neighbour".
        double cost;

        // EAP variables: track where to start the next line, and the position of the previous pruning point.
        // Must be init to 0: index 0 is the next starting index and also the "previous pruning point"
        size_t next_start{0}, prev_pp{0};

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Initialisation of the first line.
        {
            const double l0 = lines[0];
            // Fist cell is a special case.
            // Check against the original upper bound dealing with the case where we have both series of length 1.
            cost = square_dist(l0, cols[0]);
            if (cost > cutoff) { return POSITIVE_INFINITY; }
            buffers[c + 0] = cost;
            // All other cells. Checking against "ub" is OK as the only case where the last cell of this line is the
            // last alignment is taken are just above (1==nblines==nbcols, and we have nblines >= nbcols).
            size_t curr_pp = 1;
            for (j = 1; j == curr_pp && j < nbcols; ++j) {
                cost = cost + square_dist(l0, cols[j]);
                buffers[c + j] = cost;
                if (cost <= ub) { ++curr_pp; }
            }
            ++i;
            prev_pp = curr_pp;
        }

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Main loop
        for (; i < nblines; ++i) {
            // --- --- --- Swap and variables init
            std::swap(c, p);
            const double li = lines[i];
            size_t curr_pp = next_start; // Next pruning point init at the start of the line
            j = next_start;
            // --- --- --- Stage 0: Special case for the first column. Can only look up (border on the left)
            {
                cost = buffers[p + j] + square_dist(li, cols[j]);
                buffers[c + j] = cost;
                if (cost <= ub) { curr_pp = j + 1; } else { ++next_start; }
                ++j;
            }
            // --- --- --- Stage 1: Up to the previous pruning point while advancing next_start: diag and top
            for (; j == next_start && j < prev_pp; ++j) {
                cost = std::min(buffers[p + j - 1], buffers[p + j]) + square_dist(li, cols[j]);
                buffers[c + j] = cost;
                if (cost <= ub) { curr_pp = j + 1; } else { ++next_start; }
            }
            // --- --- --- Stage 2: Up to the previous pruning point without advancing next_start: left, diag and top
            for (; j < prev_pp; ++j) {
                cost = min(cost, buffers[p + j - 1], buffers[p + j]) + square_dist(li, cols[j]);
                buffers[c + j] = cost;
                if (cost <= ub) { curr_pp = j + 1; }
            }
            // --- --- --- Stage 3: At the previous pruning point. Check if we are within bounds.
            if (j < nbcols) { // If so, two cases.
                if (j == next_start) { // Case 1: Advancing next start: only diag.
                    cost = buffers[p + j - 1] + square_dist(li, cols[j]);
                    buffers[c + j] = cost;
                    if (cost <= ub) { curr_pp = j + 1; }
                    else {
                        // Special case if we are on the last alignment: return the actual cost if we are <= cutoff
                        if (i == nblines - 1 && j == nbcols - 1 && cost <= cutoff) { return cost; }
                        else { return POSITIVE_INFINITY; }
                    }
                } else { // Case 2: Not advancing next start: possible path in previous cells: left and diag.
                    cost = std::min(cost, buffers[p + j - 1]) + square_dist(li, cols[j]);
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
                cost = cost + square_dist(li, cols[j]);
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

template<bool doLA=false>
[[nodiscard]] double dtw(
        const double *series1, size_t length1,
        const double *series2, size_t length2
) {
    const auto check_result = check_order_series(series1, length1, series2, length2);
    switch (check_result.index()) {
        case 0: { return std::get<0>(check_result); }
        case 1: {
            const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Compute a cutoff point using the diagonal
            double cutoff{0};
            // We have less columns than lines: cover all the columns first.
            for (size_t i{0}; i < nbcols; ++i) { cutoff += square_dist(lines[i], cols[i]); }
            // Then go down in the last column
            if(nbcols<nblines) {
                const auto lc = cols[nbcols - 1];
                for (size_t i {nbcols}; i < nblines; ++i) { cutoff += square_dist(lines[i], lc); }
            }

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            return internal::dtw<doLA>(lines, nblines, cols, nbcols, cutoff);
        }
        default: should_not_happen();
    }
}

template<bool doLA=false>
[[nodiscard]] double dtw(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        double cutoff
) {
    const auto check_result = check_order_series(series1, length1, series2, length2);
    switch (check_result.index()) {
        case 0: { return std::get<0>(check_result);}
        case 1: {
            const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);
            return internal::dtw<doLA>(lines, nblines, cols, nbcols, cutoff);
        }
        default: should_not_happen();
    }
}



/** Classic DTW implemented on a double buffer.
 * @param series1 Pointer to the first series' values
 * @param length1 Length of the first series
 * @param series2 Pointer to the second series' values
 * @param length2 Length of the second series
 * @return DTW between the two series
 */
[[nodiscard]] double dtw_base(
        const double* series1, size_t length1,
        const double* series2, size_t length2
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
            cost = square_dist(li, cols[0]);
            buffers[c + 0] = cost;
        }
        // Rest of the line, a cell only depends on the previous cell
        for (size_t j{1}; j < nbcols; ++j) {
            cost = cost + square_dist(li, cols[j]);
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
            cost = buffers[p + 0] + square_dist(li, cols[0]);
            buffers[c + 0] = cost;
        }
        // --- --- --- Iterate through the columns
        for (size_t j{1}; j < nbcols; ++j) {
            cost = square_dist(li, cols[j]) + min(cost, buffers[p + j - 1], buffers[p + j]);
            buffers[c + j] = cost;
        }
    }

    // --- Finalisation
    return buffers[c + nbcols - 1];
}

[[nodiscard]] double dtw_base_ea(
        const double* series1, size_t length1,
        const double* series2, size_t length2,
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
            cost = square_dist(li, cols[0]);
            buffers[c + 0] = cost;
        }
        // Rest of the line, a cell only depends on the previous cell
        for (size_t j{1}; j < nbcols; ++j) {
            cost = cost + square_dist(li, cols[j]);
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
            cost = buffers[p + 0] + square_dist(li, cols[0]);
            buffers[c + 0] = cost;
            minv = cost;
        }
        // --- --- --- Iterate through the columns
        for (size_t j{1}; j < nbcols; ++j) {
            cost = square_dist(li, cols[j]) + min(cost, buffers[p + j - 1], buffers[p + j]);
            buffers[c + j] = cost;
            minv = std::min(minv, cost);
        }
        if (minv > cutoff) { return POSITIVE_INFINITY; }
    }

    // --- Finalisation
    return buffers[c + nbcols - 1];

}