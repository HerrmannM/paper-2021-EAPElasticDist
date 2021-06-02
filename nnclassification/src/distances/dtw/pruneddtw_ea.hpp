#pragma once

/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/**                                                                   **/
/** This code hase been extracted from the UCR-USP suite by           **/
/** Matthieu Herrmann. Original disclaimer reproduce as is.           **/
/**                                                                   **/
/***********************************************************************/
/**                                                                   **/
/** This suite is a modification of the UCR Suite and, consequently,  **/
/** follows the same copyright (Transcribed below, without any change)**/
/**                                                                   **/
/** This modified version is responsability of Diego Furtado Silva,   **/
/** Rafael Giusti, Eamonn Keogh and Gustavo E. A. P. A. Batista       **/
/**                                                                   **/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected © 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/

namespace internal {

#define min(x, y) ((x)<(y)?(x):(y))
#define max(x, y) ((x)>(y)?(x):(y))
#define dist(x, y) ((x-y)*(x-y))

#include <cmath>

  /// Calculate Dynamic Time Wrapping distance
  /// A,B: data and query, respectively
  /// cb : cummulative bound used for early abandoning
  /// m: length of the series
  /// r  : size of Sakoe-Chiba warpping band
  /// Code from UCR-USP with following modifications:
  /// - Renamed to pruneddtw2018
  /// - Formatting
  /// - Redefinition of INF using numeric_limits
  /// - Check & correct the window size
  double pruneddtw2018(const double* A, const double* B, const double* cb, int m, int r, double bsf) {
    using namespace std;
    constexpr double INF = numeric_limits<double>::infinity();
    if(r>=m){r=m-1;} else if(r<0){r=0;}

    double* cost;
    double* cost_prev;
    double* cost_tmp;
    int i, j;
    double x, y, z, min_cost;

    // Variables to implement the pruning - PrunedDTW
    int sc = 0, ec = 0, next_ec, lp; //lp stands for last pruning
    double UB = bsf-cb[r+1];
    bool foundSC, prunedEC = false;
    int iniJ;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(m).
    cost = (double*) malloc(sizeof(double)*(m));
    cost_prev = (double*) malloc(sizeof(double)*(m));

    for (j = 0; j<m; j++) {
      cost[j] = INF;
      cost_prev[j] = INF;
    }

    for (i = 0; i<m; i++) {

      min_cost = INF;

      foundSC = false;
      prunedEC = false;
      next_ec = i+r+1;

      iniJ = max(0, max(i-r, sc));

      for (j = iniJ; j<=min(m-1, i+r); j++) {
        /// Initialize all row and column
        if ((i==0) && (j==0)) {
          cost[j] = dist(A[0], B[0]);
          min_cost = cost[j];
          foundSC = true;
          continue;
        }

        if (j==iniJ) { y = INF; }
        else { y = cost[j-1]; }
        if ((i==0) || (j==i+r) || (j>=lp)) { x = INF; }
        else { x = cost_prev[j]; }
        if ((i==0) || (j==0) || (j>lp)) { z = INF; }
        else { z = cost_prev[j-1]; }

        /// Classic DTW calculation
        cost[j] = min(min(x, y), z)+dist(A[i], B[j]);

        /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
        if (cost[j]<min_cost) {
          min_cost = cost[j];
        }

        /// Pruning criteria
        if (!foundSC && cost[j]<=UB) {
          sc = j;
          foundSC = true;
        }

        if (cost[j]>UB) {
          if (j>ec) {
            lp = j;
            prunedEC = true;
            break;
          }
        } else {
          next_ec = j+1;
        }

      }
      if (i+r<m-1) {
        UB = bsf-cb[i+r+1];
        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (min_cost+cb[i+r+1]>=bsf) {
          free(cost);
          free(cost_prev);
          return INF;
        }
      }

      /// Move current array to previous array.
      cost_tmp = cost;
      cost = cost_prev;
      cost_prev = cost_tmp;

      if (sc>0) { cost_prev[sc-1] = INF; }
      if (!prunedEC) { lp = i+r+1; }
      ec = next_ec;
    }

    // If pruned in the last row
    if (prunedEC) { cost_prev[m-1] = INF; }

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[m-1];
    free(cost);
    free(cost_prev);
    return final_dtw;
  }

#undef min
#undef max
#undef dist
} // End of namespace internal

[[nodiscard]] double pruneddtw2018(
  const double* series1, size_t length1,
  const double* series2, size_t length2,
  const double* cb,
  size_t w,
  double bsf
) {
  const auto check_result = check_order_series(series1, length1, series2, length2);
  switch (check_result.index()) {
    case 0: { return std::get<0>(check_result); }
    case 1: {
      const auto[lines, nblines, cols, nbcols] = std::get<1>(check_result);
      return internal::pruneddtw2018(lines, cols, cb, (int) nblines, (int) w, bsf);
    }
    default: should_not_happen();
  }
}

double lb_Keogh_prunedDTW2018(
  const double* query, size_t length_query,
  const double* upper, const double* lower,
  double ub,
  double* dist_to_env
) {
  double lb{0};
  double d;

  for (size_t i = 0; i<length_query && lb<ub; i++) {
    d = 0;
    double qi{query[i]};
    if (const auto ui{upper[i]}; qi>ui) { d = square_dist(qi, ui); }
    else if (const auto li{lower[i]}; qi<li) { d = square_dist(qi, li); }
    lb += d;
    dist_to_env[i] = d;
  }
  return lb;
}

/// Cumulate "from the back", write result in 'cumulative'
void cumulate_prunedDTW2018(const double* dist_to_env, double* cumulative, int length) {
  cumulative[length-1] = dist_to_env[length-1];
  for (int i = length-2; i>=0; i--) {
    cumulative[i] = cumulative[i+1]+dist_to_env[i];
  }
}
