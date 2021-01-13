#pragma once

#include "../distances.hpp"

namespace internal {

    /** Cost function used when transforming X=(x1, x2, ... xi) into Y = (y1, ..., yj) by Split or Merge (symmetric)
     * @param new_point in either X or Y
     * @param xi Last point of X
     * @param yj Last point of Y
     * @param c cost of split and merge operation
     * @return msm cost of the xi-yj alignment (without "recursive" part)
     */
    inline double split_merge_cost(double new_point, double xi, double yj, double c) {
        if (((xi <= new_point) && (new_point <= yj)) || ((yj <= new_point) && (new_point <= xi))) {
            return c;
        } else {
            return c + std::min(std::abs(new_point - xi), std::abs(new_point - yj));
        }
    }

    /** Move Split Merge metric with cutoff point for early abandoning and pruning.
     *  Double buffered implementation using O(n) space.
     *  Worst case scenario has a O(n²) time complexity (no pruning nor early abandoning).
     *  A tight cutoff can allow a lot of pruning, speeding up the process considerably.
     *  Actual implementation assuming that some pre-conditions are fulfilled.
     * @param lines     Pointer to the "line series". Must be the longest series. Cannot be null.
     * @param nblines   Length of the line series. Must be 0 < nbcols <= nblines < tempo::MAX_SERIES_LENGTH.
     * @param cols      Pointer to the "column series". Must be the shortest series. Cannot be null.
     * @param nbcols    Length of the column series. Must be 0 < nbcols <= nblines < tempo::MAX_SERIES_LENGTH.
     * @param co        Cost of the Split and Merge operations, also added to Move operations. Must have c>=0.
     * @param cutoff.   Attempt to prune computation of alignments with cost > cutoff.
     *                  May lead to early abandoning.
     * @return MSM between the two series or +INF if early abandoned.
     */
    template<bool doLA=false>
    [[nodiscard]] double msm(
            const double *lines, size_t nblines,
            const double *cols, size_t nbcols,
            const double co,
            double cutoff
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
        double ub = initBlock{
            // The last alignment can only computed if we have nbcols >= 2
            if(nbcols>=2) {
                const auto li = lines[nblines - 1];
                const auto li1 = lines[nblines - 2];
                const auto cj = cols[nbcols - 1];
                const auto cj1 = cols[nbcols - 2];
                const auto la = min(
                        std::abs(li - cj),                  // Diag: Move
                        split_merge_cost(cj, li, cj1, co),  // Previous: Split/Merge
                        split_merge_cost(li, li1, cj, co)   // Above: Split/Merge
                );
                return double(nextafter(cutoff, POSITIVE_INFINITY) - la);
            } else {
                return double(cutoff); // Force type to prevent auto-deduction failure
            }
        };
        if constexpr (!doLA){ ub = cutoff; }

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Double buffer allocation, no initialisation required (border condition manage in the code).
        // Base indices for the 'c'urrent row and the 'p'revious row.
        auto buffers = std::unique_ptr<double[]>(new double[nbcols * 2]);
        size_t c{0}, p{nbcols};

        // Line & column counters
        size_t i{0}, j{0};

        // Cost accumulator. Also used as the "left neighbour".
        double cost{0};

        // EAP variables: track where to start the next line, and the position of the previous pruning point.
        // Must be init to 0: index 0 is the next starting index and also the "previous pruning point"
        size_t next_start{0}, prev_pp{0};

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Initialisation: compute the first line. Required as the main loop starts at line=1, not 0.
        {
            const auto l0 = lines[0];
            // First cell (0,0) is a special case. Early abandon if above the cut-off point.
            {
                cost = std::abs(l0 - cols[0]); // Very first cell
                buffers[c + 0] = cost;
                if (cost <= ub) { prev_pp = 1; } else { return POSITIVE_INFINITY; }
            }
            // Rest of the line, a cell only depends on the previous cell. Stop when > ub, update prev_pp.
            for (j=1; j < nbcols; ++j) {
                cost = cost + split_merge_cost(cols[j], l0, cols[j - 1], co);
                if (cost <= ub) { buffers[c + j] = cost; prev_pp = j + 1; } else { break; }
            }
            // Next line.
            ++i;
        }

        // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        // Main loop
        for (; i < nblines; ++i) {
            // --- --- --- Swap and variables init
            std::swap(c, p);
            const double li = lines[i];
            const double li1 = lines[i - 1];
            size_t curr_pp = next_start; // Next pruning point init at the start of the line
            j = next_start;
            // --- --- --- Stage 0: Special case for the first column. Can only look up (border on the left)
            {
                cost = buffers[p + j] + split_merge_cost(li, li1, cols[j], co);
                buffers[c + j] = cost;
                if (cost <= ub) { curr_pp = j + 1; } else { ++next_start; }
                ++j;
            }
            // --- --- --- Stage 1: Up to the previous pruning point while advancing next_start: diag and top
            for (; j == next_start && j < prev_pp; ++j) {
                const double cj = cols[j];
                cost = std::min(
                        buffers[p + j - 1] + std::abs(li - cj),             // Diag: Move
                        buffers[p + j] + split_merge_cost(li, li1, cj, co)  // Above: Split/Merge
                );
                buffers[c + j] = cost;
                if (cost <= ub) { curr_pp = j + 1; } else { ++next_start; }
            }
            // --- --- --- Stage 2: Up to the previous pruning point without advancing next_start: left, diag and top
            for (; j < prev_pp; ++j) {
                const double cj = cols[j];
                cost = min(
                        buffers[p + j - 1] + std::abs(li - cj),               // Diag: Move
                        cost + split_merge_cost(cj, li, cols[j - 1], co),     // Previous: Split/Merge
                        buffers[p + j] + split_merge_cost(li, li1, cj, co)    // Above: Split/Merge
                );
                buffers[c + j] = cost;
                if (cost <= ub) { curr_pp = j + 1; }
            }
            // --- --- --- Stage 3: At the previous pruning point. Check if we are within bounds.
            if (j < nbcols) { // If so, two cases.
                const double cj = cols[j];
                if (j == next_start) { // Case 1: Advancing next start: only diag.
                    cost = buffers[p + j - 1] + std::abs(li - cj);            // Diag: Move
                    buffers[c + j] = cost;
                    if (cost <= ub) { curr_pp = j + 1; }
                    else {
                        // Special case if we are on the last alignment: return the actual cost if we are <= cutoff
                        if (i == nblines - 1 && j == nbcols - 1 && cost <= cutoff) { return cost; }
                        else { return POSITIVE_INFINITY; }
                    }
                } else { // Case 2: Not advancing next start: possible path in previous cells: left and diag.
                    cost = std::min(
                            buffers[p + j - 1] + std::abs(li - cj),               // Diag: Move
                            cost + split_merge_cost(cj, li, cols[j - 1], co)      // Previous: Split/Merge
                    );
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
                cost = cost + split_merge_cost(cols[j], li, cols[j - 1], co);      // Previous: Split/Merge
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
[[nodiscard]] double msm(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        const double co
) {
    const auto check_result = check_order_series(series1, length1, series2, length2);
    switch (check_result.index()) {
        case 0: { return std::get<0>(check_result);}
        case 1: {
            const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Compute a cutoff point using the diagonal.
            double cutoff{0};
            // We have less columns than lines: cover all the columns first.
            for (size_t i{0}; i < nbcols; ++i) { cutoff += std::abs(lines[i] - cols[i]); } // Diag: Move
            // Then go down in the last column
            if(nbcols<nblines) {
                const auto lc = cols[nbcols - 1];
                for (size_t i{nbcols}; i < nblines; ++i) {
                    cutoff += internal::split_merge_cost(lines[i], lines[i - 1], lc, co);   // Above: Split/Merge
                }
            }

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            return internal::msm<doLA>(lines, nblines, cols, nbcols, co, cutoff);
        }
        default: should_not_happen();
    }
}


template<bool doLA=false>
[[nodiscard]] double msm(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        const double co,
        double cutoff
) {
    const auto check_result = check_order_series(series1, length1, series2, length2);
    switch (check_result.index()) {
        case 0: { return std::get<0>(check_result);}
        case 1: {
            const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);
            return internal::msm<doLA>(lines, nblines, cols, nbcols, co, cutoff);
        }
        default: should_not_happen();
    }
}



/** Move-Split-Merge metric, implemented on a double buffer.
 * @param series1 Pointer to the first series' values
 * @param length1 Length of the first series
 * @param series2 Pointer to the second series' values
 * @param length2 Length of the second series
 * @param co Cost of the Split and Merge operations, also added to Move operations. Must have c>=0.
 * @return MSM cost
 */
[[nodiscard]] static double msm_base (
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        double co
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

    // Double buffer allocation initialized to POSITIVE_INFINITY.
    // Base indices for the 'c'urrent row and the 'p'revious row.
    std::vector<double> buffers_v(nbcols * 2, POSITIVE_INFINITY);
    double *buffers = buffers_v.data();
    size_t c{0}, p{nbcols};

    // --- Accesses
    const double c0 = cols[0]; // Use to compute the column border

    // --- Initialisation: compute the first line, implicitly dealing with the line border condition
    {
        const auto l0 = lines[0];
        // Column border is handle in the main loop
        buffers[c + 0] = std::abs(l0 - c0); // Very first cell
        for (size_t j{1}; j < nbcols; ++j) {
            buffers[c + j] = buffers[c + j - 1] + internal::split_merge_cost(cols[j], l0, cols[j - 1], co);
        }
    }

    // --- Main loop
    for (size_t i{1}; i < nblines; ++i) {
        // --- --- --- Swap and variables init
        std::swap(c, p);
        const double li = lines[i];
        const double li1 = lines[i - 1];
        // --- --- --- Init the border (very first column)
        buffers[c + 0] = buffers[p + 0] + internal::split_merge_cost(li, li1, c0, co);
        // --- --- --- Iterate through the columns (start at the 2nd column)
        for (size_t j{1}; j < nbcols; ++j) {
            const double cj = cols[j];
            buffers[c + j] = min(
                    buffers[p + j - 1] + std::abs(li - cj),                                     // Diag: Move
                    buffers[c + j - 1] + internal::split_merge_cost(cj, li, cols[j - 1], co),   // Previous: Split/Merge
                    buffers[p + j] + internal::split_merge_cost(li, li1, cj, co)                // Above: Split/Merge
            );
        } // End for over columns
    }

    // --- Finalisation
    return buffers[c + nbcols - 1];
}

[[nodiscard]] static double msm_base_ea (
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        double co,
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

    // Double buffer allocation initialized to POSITIVE_INFINITY.
    // Base indices for the 'c'urrent row and the 'p'revious row.
    std::vector<double> buffers_v(nbcols * 2, POSITIVE_INFINITY);
    double *buffers = buffers_v.data();
    size_t c{0}, p{nbcols};
    double minv = POSITIVE_INFINITY;

    // --- Accesses
    const double c0 = cols[0]; // Use to compute the column border

    // --- Initialisation: compute the first line, implicitly dealing with the line border condition
    {
        const auto l0 = lines[0];
        // Column border is handle in the main loop
        buffers[c + 0] = std::abs(l0 - c0); // Very first cell
        for (size_t j{1}; j < nbcols; ++j) {
            buffers[c + j] = buffers[c + j - 1] + internal::split_merge_cost(cols[j], l0, cols[j - 1], co);
        }
    }

    // --- Main loop
    for (size_t i{1}; i < nblines; ++i) {
        // --- --- --- Swap and variables init
        std::swap(c, p);
        const double li = lines[i];
        const double li1 = lines[i - 1];
        // --- --- --- Init the border (very first column)
        buffers[c + 0] = buffers[p + 0] + internal::split_merge_cost(li, li1, c0, co);
        minv = buffers_v[c+0];
        // --- --- --- Iterate through the columns (start at the 2nd column)
        for (size_t j{1}; j < nbcols; ++j) {
            const double cj = cols[j];
            buffers[c + j] = min(
                    buffers[p + j - 1] + std::abs(li - cj),                                     // Diag: Move
                    buffers[c + j - 1] + internal::split_merge_cost(cj, li, cols[j - 1], co),   // Previous: Split/Merge
                    buffers[p + j] + internal::split_merge_cost(li, li1, cj, co)                // Above: Split/Merge
            );
            minv = std::min(minv, buffers[c+j]);
        } // End for over columns
        if (minv > cutoff) { return POSITIVE_INFINITY; }
    }

    // --- Finalisation
    return buffers[c + nbcols - 1];
}
