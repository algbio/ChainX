#ifndef CHAINX_ALGO_H
#define CHAINX_ALGO_H

#include <iostream>
#include <cassert>
#include <vector>
#include <tuple>
#include <algorithm>
#include <zlib.h>  
#include <string>
#include <chrono>

#undef VERBOSE
#define VERBOSE 0

namespace chainx
{
  /**
   * @brief   compute anchor-restricted edit distance using strong precedence criteria
   * 			    optimized to run faster using engineering trick(s), comparison mode: global
   **/
  int compute_global(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    int n = anchors.size();
    std::vector<int> costs(n, 0);

    int bound_redit = 100; //distance assumed to be <= 100
    int revisions = 0;
    //with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

    while (true) 
    {
      int inner_loop_start = 0;

      for(int j=1; j<n; j++)
      {
        //compute cost[i] here
        int find_min_cost = std::numeric_limits<int>::max();

        int j_a = std::get<0>(anchors[j]);
        int j_b = std::get<0>(anchors[j]) + std::get<2>(anchors[j]) - 1;
        int j_c = std::get<1>(anchors[j]);
        int j_d = std::get<1>(anchors[j]) + std::get<2>(anchors[j]) - 1;

        // anchor i < anchor j 

        while (j_a - std::get<0>(anchors[inner_loop_start]) - 1 > bound_redit)
          inner_loop_start++;

        for(int i=j-1; i>=inner_loop_start; i--)
        {
          int i_a = std::get<0>(anchors[i]);
          int i_b = std::get<0>(anchors[i]) + std::get<2>(anchors[i]) - 1;
          int i_c = std::get<1>(anchors[i]);
          int i_d = std::get<1>(anchors[i]) + std::get<2>(anchors[i]) - 1;

          if (costs[i] < std::numeric_limits<int>::max() && i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
          {
            int gap1 = std::max(0, j_a - i_b - 1);
            int gap2 = std::max(0, j_c - i_d - 1);
            int g = std::max(gap1,gap2);

            int overlap1 = std::max(0, i_b - j_a + 1);
            int overlap2 = std::max(0, i_d - j_c + 1);
            int o = std::abs(overlap1 - overlap2);

            find_min_cost = std::min(find_min_cost, costs[i] + g + o);
          }
        }
        //save optimal cost at offset j
        costs[j] = find_min_cost;
      }

      if (costs[n-1] > bound_redit)
      {
        bound_redit = bound_redit * 4;
        revisions++;
      }
      else
        break;
    }

    if (VERBOSE)
      std::cerr << "Cost array = " << costs << "\n";

    if (VERBOSE)
      std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
    return costs[n-1];
  }

  /**
   * @brief   version of compute_global that finds the optimal chain using diagonal distance
   **/
  int compute_global_optimal(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    // anchors are sorted by starting position in first sequence (the reference? TODO check)
    const int n = anchors.size();
    std::vector<int> costs(n, 0);
    std::vector<int> start_endpoints; // sorted (index) list of anchors startpoints and endpoints together
    start_endpoints.reserve(2*n);
    std::vector<int> diagonal(n, 0); // diagonal (index) of each anchor
    std::vector<int> diagonal_bucket_value;
    int d = -1; // number of distinct diagonals

    int bound_redit = 100; //distance assumed to be <= 100
    int revisions = 0;
    //with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

    // sort anchors by start and endpoint
    for (int j=0; j<n; j++)
    {
        start_endpoints.push_back(j);
        start_endpoints.push_back(-j); // negative index means endpoint
    }
    std::sort (start_endpoints.begin(), start_endpoints.end(),
        [&](const int i,
          const int j) -> bool
        {
        return (((i >= 0) ? std::get<0>(anchors[i]) : std::get<0>(anchors[-i]) + std::get<2>(anchors[-i]) - 1) <
                ((j >= 0) ? std::get<0>(anchors[j]) : std::get<0>(anchors[-j]) + std::get<2>(anchors[-j]) - 1));
        });

    // sort anchors by diagonal to figure out each anchor's diagonal index
    {
        // TODO test if unordered_set is faster
        std::vector<int> diagonal_order(n);
        for (int j=0; j<n; j++)
            diagonal_order[j]=j;
        std::sort (diagonal_order.begin(), diagonal_order.end(),
        [&](const int i,
          const int j) -> bool
        {
        return (std::get<0>(anchors[i]) - std::get<1>(anchors[i])) <
               (std::get<0>(anchors[j]) - std::get<1>(anchors[j]));
        });
        int dd = -1;
        int prec_diagonal = std::numeric_limits<int>::max();
        for (int k=0; k<n; k++)
        {
            const int curr_diagonal = std::get<0>(anchors[diagonal_order[k]]) - std::get<1>(anchors[diagonal_order[k]]);
            if (curr_diagonal != prec_diagonal)
            {
                dd += 1;
                diagonal_bucket_value.push_back(curr_diagonal);
            }
            diagonal[diagonal_order[k]] = dd;
            prec_diagonal = curr_diagonal;
        }
        d = dd + 1;
    }

    std::vector<int> active_anchor(d, -1); // active anchor per diagonal
    int stats_overlap_comparisons = 0;
    while (true)
    {
      int inner_loop_start = 0;
      stats_overlap_comparisons = 0;

      // main loop
      for(int j=2; j<2*n; j++)
      {
        if (start_endpoints[j] >= 0) // if startpoint
        {
          const int anchorj = start_endpoints[j];
          const int curr_diagonal = std::get<0>(anchors[anchorj]) - std::get<1>(anchors[anchorj]);
          int find_min_cost = std::numeric_limits<int>::max();
          int find_min_cost_gap = std::numeric_limits<int>::max();
          int find_min_cost_overlap = std::numeric_limits<int>::max();

          // handle start of diagonal
          if (active_anchor[diagonal[anchorj]] >= 0) {
            find_min_cost = costs[active_anchor[diagonal[anchorj]]];
          }
          active_anchor[diagonal[anchorj]] = anchorj;

          int j_a = std::get<0>(anchors[anchorj]);
          int j_b = std::get<0>(anchors[anchorj]) + std::get<2>(anchors[anchorj]) - 1;
          int j_c = std::get<1>(anchors[anchorj]);
          int j_d = std::get<1>(anchors[anchorj]) + std::get<2>(anchors[anchorj]) - 1;

          // anchor i < anchor j

          while (inner_loop_start < j &&
                  (start_endpoints[inner_loop_start] > 0 ||
                    j_a - (std::get<0>(anchors[-start_endpoints[inner_loop_start]]) + std::get<2>(anchors[-start_endpoints[inner_loop_start]]) - 1) - 1 > bound_redit))
            inner_loop_start++;

          // process anchors with gap in first sequence
          for(int i=j-1; i>=inner_loop_start; i--)
          {
            if (start_endpoints[i] > 0)
              continue;
            const int anchori = -start_endpoints[i];
            int i_a = std::get<0>(anchors[anchori]);
            int i_b = std::get<0>(anchors[anchori]) + std::get<2>(anchors[anchori]) - 1;
            int i_c = std::get<1>(anchors[anchori]);
            int i_d = std::get<1>(anchors[anchori]) + std::get<2>(anchors[anchori]) - 1;

            if (costs[anchori] < std::numeric_limits<int>::max() && i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
            {
              int gap1 = std::max(0, j_a - i_b - 1);
              int gap2 = std::max(0, j_c - i_d - 1);
              int g = std::max(gap1,gap2);

              int overlap1 = std::max(0, i_b - j_a + 1);
              int overlap2 = std::max(0, i_d - j_c + 1);
              int o = std::abs(overlap1 - overlap2);

              find_min_cost_gap = std::min(find_min_cost_gap, costs[anchori] + g + o);
            }
          }

          // process anchors overlapping in first sequence
          for (int dd=diagonal[anchorj]+1; dd<d; dd++)
          {
            const int diagonal_distance = std::abs(curr_diagonal - diagonal_bucket_value[dd]);
            if (diagonal_distance > bound_redit)
              break;
            if (active_anchor[dd] != -1)
            {
              if (costs[active_anchor[dd]] < std::numeric_limits<int>::max())
              {
                stats_overlap_comparisons+=1;
                find_min_cost_overlap = std::min(find_min_cost_overlap, costs[active_anchor[dd]] + diagonal_distance);
              }
            }
          }
          for (int dd=diagonal[anchorj]-1; dd>0; dd--)
          {
            const int diagonal_distance = std::abs(curr_diagonal - diagonal_bucket_value[dd]);
            if (diagonal_distance > bound_redit)
              break;
            if (active_anchor[dd] != -1)
            {
              if (costs[active_anchor[dd]] < std::numeric_limits<int>::max())
              {
                stats_overlap_comparisons+=1;
                find_min_cost_overlap = std::min(find_min_cost_overlap, costs[active_anchor[dd]] + diagonal_distance);
              }
            }
          }

          //save optimal cost at offset j
          costs[anchorj] = std::min(find_min_cost_gap, find_min_cost_overlap);
          costs[anchorj] = std::min(costs[anchorj], find_min_cost);
        } else { // if start_endpoints[j] < 0
          const int anchorj = -start_endpoints[j];
          active_anchor[diagonal[anchorj]] = -1;
        }
      }

      if (costs[n-1] > bound_redit)
      {
        bound_redit = bound_redit * 4;
        revisions++;
      }
      else
        break;
    }

    if (VERBOSE)
      std::cerr << "Cost array = " << costs << "\n";

    if (VERBOSE)
      std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
    return costs[n-1];
  }

  /**
   * @brief   compute anchor-restricted (semi-global) edit distance using strong precedence criteria
   * 			    optimized to run faster using engineering trick(s)
   **/
  int compute_semiglobal(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    int n = anchors.size();
    std::vector<int> costs(n, 0);

    int bound_redit = 100; //distance assumed to be <= 100
    int revisions = 0;
    //with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

    while (true) 
    {
      int inner_loop_start = 0;

      for(int j=1; j<n; j++)
      {
        //compute cost[i] here
        int find_min_cost = std::numeric_limits<int>::max();

        int j_a = std::get<0>(anchors[j]);
        int j_b = std::get<0>(anchors[j]) + std::get<2>(anchors[j]) - 1;
        int j_c = std::get<1>(anchors[j]);
        int j_d = std::get<1>(anchors[j]) + std::get<2>(anchors[j]) - 1;

        // anchor i < anchor j 

        while (j_a - std::get<0>(anchors[inner_loop_start]) - 1 > bound_redit)
          inner_loop_start++;

        {
          //always consider the first dummy anchor 
          //connection to first dummy anchor is done with modified cost to allow free gaps
          int i_d = std::get<1>(anchors[0]) + std::get<2>(anchors[0]) - 1;
          int qry_gap = j_c - i_d - 1;
          find_min_cost = std::min(find_min_cost, costs[0] + qry_gap);
        }

        //process all anchors in array for the final last dummy anchor
        if (j == n-1) inner_loop_start=0;

        for(int i=j-1; i>=inner_loop_start; i--)
        {
          int i_a = std::get<0>(anchors[i]);
          int i_b = std::get<0>(anchors[i]) + std::get<2>(anchors[i]) - 1;
          int i_c = std::get<1>(anchors[i]);
          int i_d = std::get<1>(anchors[i]) + std::get<2>(anchors[i]) - 1;

          if (costs[i] < std::numeric_limits<int>::max() && i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
          {
            int gap1 = std::max(0, j_a - i_b - 1);
            int gap2 = std::max(0, j_c - i_d - 1);

            if (j == n-1) gap1=0; //modified cost for the last dummy anchor to allow free gaps
            int g = std::max(gap1,gap2);

            int overlap1 = std::max(0, i_b - j_a + 1);
            int overlap2 = std::max(0, i_d - j_c + 1);
            int o = std::abs(overlap1 - overlap2);

            find_min_cost = std::min(find_min_cost, costs[i] + g + o);
          }
        }

        //save optimal cost at offset j
        costs[j] = find_min_cost;
      }

      if (costs[n-1] > bound_redit)
      {
        bound_redit = bound_redit * 4;
        revisions++;
      }
      else
        break;
    }

    if (VERBOSE)
      std::cerr << "Cost array = " << costs << "\n";

    if (VERBOSE)
      std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
    return costs[n-1];
  }

  /**
   * @brief   version of compute_semiglobal that finds the optimal chain using diagonal
   * 		distance
   **/
  int compute_semiglobal_optimal(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    // anchors are sorted by starting position in first sequence (the reference? TODO check)
    const int n = anchors.size();
    std::vector<int> costs(n, 0);
    std::vector<int> start_endpoints; // sorted (index) list of anchors startpoints and endpoints together
    start_endpoints.reserve(2*n);
    std::vector<int> diagonal(n, 0); // diagonal (index) of each anchor
    std::vector<int> diagonal_bucket_value;
    int d = -1; // number of distinct diagonals

    int bound_redit = 100; //distance assumed to be <= 100
    int revisions = 0;
    //with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

    // sort anchors by start and endpoint
    for (int j=0; j<n; j++)
    {
        start_endpoints.push_back(j);
        start_endpoints.push_back(-j); // negative index means endpoint
    }
    std::sort (start_endpoints.begin(), start_endpoints.end(),
        [&](const int i,
          const int j) -> bool
        {
        return (((i >= 0) ? std::get<0>(anchors[i]) : std::get<0>(anchors[-i]) + std::get<2>(anchors[-i]) - 1) <
                ((j >= 0) ? std::get<0>(anchors[j]) : std::get<0>(anchors[-j]) + std::get<2>(anchors[-j]) - 1));
        });

    // sort anchors by diagonal to figure out each anchor's diagonal index
    {
        // TODO test if unordered_set is faster for bucketing
        std::vector<int> diagonal_order(n);
        for (int j=0; j<n; j++)
            diagonal_order[j]=j;
        std::sort (diagonal_order.begin(), diagonal_order.end(),
        [&](const int i,
          const int j) -> bool
        {
        return (std::get<0>(anchors[i]) - std::get<1>(anchors[i])) <
               (std::get<0>(anchors[j]) - std::get<1>(anchors[j]));
        });
        int dd = -1;
        int prec_diagonal = std::numeric_limits<int>::max();
        for (int k=0; k<n; k++)
        {
            const int curr_diagonal = std::get<0>(anchors[diagonal_order[k]]) - std::get<1>(anchors[diagonal_order[k]]);
            if (curr_diagonal != prec_diagonal)
            {
                dd += 1;
                diagonal_bucket_value.push_back(curr_diagonal);
            }
            diagonal[diagonal_order[k]] = dd;
            prec_diagonal = curr_diagonal;
        }
        d = dd + 1;
    }

    std::vector<int> active_anchor(d, -1); // active anchor per diagonal
    while (true)
    {
      int inner_loop_start = 0;

      // main loop
      for(int j=2; j<2*n; j++)
      {
        if (start_endpoints[j] >= 0) // if startpoint
        {
          const int anchorj = start_endpoints[j];
          const int curr_diagonal = std::get<0>(anchors[anchorj]) - std::get<1>(anchors[anchorj]);
          int find_min_cost = std::numeric_limits<int>::max();
          int find_min_cost_gap = std::numeric_limits<int>::max();
          int find_min_cost_overlap = std::numeric_limits<int>::max();

          // handle start of diagonal
          if (active_anchor[diagonal[anchorj]] >= 0) {
            find_min_cost = costs[active_anchor[diagonal[anchorj]]];
          }
          active_anchor[diagonal[anchorj]] = anchorj;

          int j_a = std::get<0>(anchors[anchorj]);
          int j_b = std::get<0>(anchors[anchorj]) + std::get<2>(anchors[anchorj]) - 1;
          int j_c = std::get<1>(anchors[anchorj]);
          int j_d = std::get<1>(anchors[anchorj]) + std::get<2>(anchors[anchorj]) - 1;

          // anchor i < anchor j

          while (inner_loop_start < j &&
                  (start_endpoints[inner_loop_start] > 0 ||
                    j_a - (std::get<0>(anchors[-start_endpoints[inner_loop_start]]) + std::get<2>(anchors[-start_endpoints[inner_loop_start]]) - 1) - 1 > bound_redit))
            inner_loop_start++;

          {
            //always consider the first dummy anchor 
            //connection to first dummy anchor is done with modified cost to allow free gaps
            int i_d = std::get<1>(anchors[0]) + std::get<2>(anchors[0]) - 1;
            int qry_gap = j_c - i_d - 1;
            find_min_cost = std::min(find_min_cost, costs[0] + qry_gap);
          }

          //process all anchors in array for the final last dummy anchor
          if (anchorj == n-1) inner_loop_start=0;

          // process anchors with gap in first sequence
          for(int i=j-1; i>=inner_loop_start; i--)
          {
            if (start_endpoints[i] > 0)
              continue;
            const int anchori = -start_endpoints[i];
            int i_a = std::get<0>(anchors[anchori]);
            int i_b = std::get<0>(anchors[anchori]) + std::get<2>(anchors[anchori]) - 1;
            int i_c = std::get<1>(anchors[anchori]);
            int i_d = std::get<1>(anchors[anchori]) + std::get<2>(anchors[anchori]) - 1;

            if (costs[anchori] < std::numeric_limits<int>::max() && i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
            {
              int gap1 = std::max(0, j_a - i_b - 1);
              int gap2 = std::max(0, j_c - i_d - 1);
	      if (anchorj == n-1) gap1=0; //modified cost for the last dummy anchor to allow free gaps

              int g = std::max(gap1,gap2);

              int overlap1 = std::max(0, i_b - j_a + 1);
              int overlap2 = std::max(0, i_d - j_c + 1);
              int o = std::abs(overlap1 - overlap2);

              find_min_cost_gap = std::min(find_min_cost_gap, costs[anchori] + g + o);
            }
          }

          // process anchors overlapping in first sequence
          for (int dd=diagonal[anchorj]+1; dd<d; dd++)
          {
            const int diagonal_distance = std::abs(curr_diagonal - diagonal_bucket_value[dd]);
            if (diagonal_distance > bound_redit)
              break;
            if (active_anchor[dd] != -1 && costs[active_anchor[dd]] < std::numeric_limits<int>::max())
            {
              find_min_cost_overlap = std::min(find_min_cost_overlap, costs[active_anchor[dd]] + diagonal_distance);
            }
          }
          for (int dd=diagonal[anchorj]-1; dd>0; dd--)
          {
            const int diagonal_distance = std::abs(curr_diagonal - diagonal_bucket_value[dd]);
            if (diagonal_distance > bound_redit)
              break;
            if (active_anchor[dd] != -1 && costs[active_anchor[dd]] < std::numeric_limits<int>::max())
            {
              find_min_cost_overlap = std::min(find_min_cost_overlap, costs[active_anchor[dd]] + diagonal_distance);
            }
          }

          //save optimal cost at offset j
          costs[anchorj] = std::min(find_min_cost_gap, find_min_cost_overlap);
          costs[anchorj] = std::min(costs[anchorj], find_min_cost); // edge cases
        } else { // if endpoint start_endpoints[j] < 0
          const int anchorj = -start_endpoints[j];
	  active_anchor[diagonal[anchorj]] = -1;
        }
      }

      if (costs[n-1] > bound_redit)
      {
        bound_redit = bound_redit * 4;
        revisions++;
      }
      else
        break;
    }

    if (VERBOSE)
      std::cerr << "Cost array = " << costs << "\n";

    if (VERBOSE)
      std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
    return costs[n-1];
  }

  /**
   * @brief   compute anchor-restricted edit distance using standard edit-distance like dynamic programming 
   **/
  int DP_global(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    int n = anchors.size();

    //get sequence lengths from end dummy anchor
    //assuming anchors are already sorted
    int len_ref = std::get<0>(anchors[n-1]);
    int len_qry = std::get<1>(anchors[n-1]);

    //initialize a boolean matrix (len_ref+1 x len_qry+1)
    //we will offset by 1 to be consistent with DP matrix
    std::vector<std::vector<bool> > matchAllowed(len_ref+1);
    for(int i=0; i<len_ref+1; i++) matchAllowed[i] = std::vector<bool>(len_qry+1, false);
    for(int i = 0; i<n-1; i++) //use all (except end dummy) anchors
    {
      int e_a = std::get<0>(anchors[i]);
      int e_c = std::get<1>(anchors[i]);
      int e_len = std::get<2>(anchors[i]);

      for(int j=0; j<e_len; j++)
      {
        matchAllowed[e_a+1  + j][e_c+1  + j] = true;
      }
    }


    //initialize dp_matrix (len_ref+1 x len_qry+1) 
    std::vector<std::vector<int> > dp_matrix(len_ref+1);
    for(int i=0; i<=len_ref; i++) dp_matrix[i] = std::vector<int>(len_qry+1);
    for(int i=0; i<=len_ref; i++) dp_matrix[i][0] = i;
    for(int j=0; j<=len_qry; j++) dp_matrix[0][j] = j;

    for(int i=1; i<=len_ref; i++)
    {
      for(int j=1; j<=len_qry; j++)
      {
        if (matchAllowed[i][j])
          dp_matrix[i][j] = std::min({dp_matrix[i - 1][j - 1], dp_matrix[i - 1][j] + 1, dp_matrix[i][j - 1] + 1});
        else
          dp_matrix[i][j] = std::min({dp_matrix[i - 1][j - 1] + 1, dp_matrix[i - 1][j] + 1, dp_matrix[i][j - 1] + 1});
      }
    }

    return dp_matrix[len_ref][len_qry];
  }

  /**
   * @brief   compute anchor-restricted (semi-global) edit distance using standard edit-distance like dynamic programming 
   **/
  int DP_semiglobal(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    int n = anchors.size();

    //get sequence lengths from end dummy anchor
    //assuming anchors are already sorted
    int len_ref = std::get<0>(anchors[n-1]);
    int len_qry = std::get<1>(anchors[n-1]);

    //initialize a boolean matrix (len_ref+1 x len_qry+1)
    //we will offset by 1 to be consistent with DP matrix
    std::vector<std::vector<bool> > matchAllowed(len_ref+1);
    for(int i=0; i<len_ref+1; i++) matchAllowed[i] = std::vector<bool>(len_qry+1, false);
    for(int i = 0; i<n-1; i++) //use all (except end dummy) anchors
    {
      int e_a = std::get<0>(anchors[i]);
      int e_c = std::get<1>(anchors[i]);
      int e_len = std::get<2>(anchors[i]);

      for(int j=0; j<e_len; j++)
      {
        matchAllowed[e_a+1  + j][e_c+1  + j] = true;
      }
    }


    //initialize dp_matrix (len_ref+1 x len_qry+1) 
    std::vector<std::vector<int> > dp_matrix(len_ref+1);
    for(int i=0; i<=len_ref; i++) dp_matrix[i] = std::vector<int>(len_qry+1);
    for(int i=0; i<=len_ref; i++) dp_matrix[i][0] = 0; //changed from i to 0 for free gaps
    for(int j=0; j<=len_qry; j++) dp_matrix[0][j] = j;

    for(int i=1; i<=len_ref; i++)
    {
      for(int j=1; j<=len_qry; j++)
      {
        if (matchAllowed[i][j])
          dp_matrix[i][j] = std::min({dp_matrix[i - 1][j - 1], dp_matrix[i - 1][j] + 1, dp_matrix[i][j - 1] + 1});
        else
          dp_matrix[i][j] = std::min({dp_matrix[i - 1][j - 1] + 1, dp_matrix[i - 1][j] + 1, dp_matrix[i][j - 1] + 1});
      }
    }

    int final_distance = std::numeric_limits<int>::max();
    for(int i=0; i<=len_ref; i++) final_distance = std::min (final_distance, dp_matrix[i][len_qry]);
    return final_distance;
  }
}

#endif
