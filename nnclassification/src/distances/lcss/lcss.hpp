#pragma once

#include "../distances.hpp"

namespace internal {
    /// Check if two double numbers are within EPSILON (1 = similar, 0 = not similar)
    [[nodiscard]] bool sim(double a, double b, double e) { return std::fabs(a - b) < e; }
} // End of namespace internal


[[nodiscard]] double lcss(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        double epsilon,
        size_t w
) {
    const auto check_result = check_order_series(series1, length1, series2, length2);
    switch (check_result.index()) {
        case 0: { return std::get<0>(check_result); }
        case 1: {
            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);
            // Cap the windows and check that, given the constraint, an alignment is possible
            if (w > nblines) { w = nblines; }
            if (nblines - nbcols > w) { return POSITIVE_INFINITY; }

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Double buffer allocation, no initialisation required (border condition manage in the code).
            // Base indices for the 'c'urrent row and the 'p'revious row. Account for the extra cell (+1 and +2)
            std::vector<double> buffers_v((1 + nbcols) * 2, 0);
            double *buffers = buffers_v.data();
            size_t c{0 + 1}, p{nbcols + 2};

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Initialisation: OK, border line and "first diag" init to 0

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Main loop
            for (size_t i{0}; i < nblines; ++i) {
                // --- --- --- Swap and variables init
                std::swap(c, p);
                const double li = lines[i];
                const size_t jStart = cap_start_index_to_window(i, w);
                const size_t jStop = cap_stop_index_to_window_or_end(i, w, nbcols);
                // --- --- --- Init the border (very first column)
                buffers[c + jStart - 1] = 0;
                // --- --- --- Iterate through the columns
                for (size_t j{jStart}; j < jStop; ++j) {
                    if (internal::sim(li, cols[j], epsilon)) {
                        buffers[c + j] = buffers[p + j - 1] + 1; // Diag + 1
                    } else { // Note: Diagonal lookup required, e.g. when w=0
                        buffers[c + j] = max(buffers[c + j - 1], buffers[p + j - 1], buffers[p + j]);
                    }
                }
            }

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Finalisation: put the result on a [0 - 1] range
            return 1.0 - (double(buffers[c + nbcols - 1]) / nbcols);
        }
        default: should_not_happen();
    }
}


[[nodiscard]] double lcss(
        const double *series1, size_t length1,
        const double *series2, size_t length2,
        double epsilon,
        size_t w,
        double cutoff
) {
    const auto check_result = check_order_series(series1, length1, series2, length2);
    switch (check_result.index()) {
        case 0: { return std::get<0>(check_result); }
        case 1: {
            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);
            // Cap the windows and check that, given the constraint, an alignment is possible
            if (w > nblines) { w = nblines; }
            if (nblines - nbcols > w) { return POSITIVE_INFINITY; }

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Double buffer allocation, init to 0.
            // Base indices for the 'c'urrent row and the 'p'revious row. Account for the extra cell (+1 and +2)
            std::vector<double> buffers_v((1 + nbcols) * 2, 0);
            double *buffers = buffers_v.data();
            size_t c{0 + 1}, p{nbcols + 2};

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Score to reach to equal ub, to beat to do better
            if (cutoff > 1) { cutoff = 1; }
            const size_t to_reach = std::floor((1 - cutoff) * nbcols);
            size_t current_max = 0;

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Initialisation: OK, border line and "first diag" init to 0

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Main loop
            for (size_t i{0}; i < nblines; ++i) {
                // --- --- --- Stop if not enough remaining lines to reach the target (by taking the diagonal)
                const size_t lines_left = nblines - i;
                if (current_max + lines_left < to_reach) { return POSITIVE_INFINITY; }
                // --- --- --- Swap and variables init
                std::swap(c, p);
                const double li = lines[i];
                const size_t jStart = cap_start_index_to_window(i, w);
                const size_t jStop = cap_stop_index_to_window_or_end(i, w, nbcols);
                // --- --- --- Init the border (very first column)
                buffers[c + jStart - 1] = 0;
                // --- --- --- Iterate through the columns
                for (size_t j{jStart}; j < jStop; ++j) {
                    if (internal::sim(li, cols[j], epsilon)) {
                        const size_t cost = buffers[p + j - 1] + 1; // Diag + 1
                        current_max = std::max(current_max, cost);
                        buffers[c + j] = cost;
                    } else { // Note: Diagonal lookup required, e.g. when w=0
                        buffers[c + j] = max(buffers[c + j - 1], buffers[p + j - 1], buffers[p + j]);
                    }
                }
            }

            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // Finalisation: put the result on a [0 - 1] range
            return 1.0 - (double(buffers[c + nbcols - 1]) / nbcols);
        }
        default: should_not_happen();
    }
}
