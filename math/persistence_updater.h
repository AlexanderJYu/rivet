/**********************************************************************
Copyright 2014-2016 The RIVET Developers. See the COPYRIGHT file at
the top-level directory of this distribution.

This file is part of RIVET.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
/**
 * \class	PersistenceUpdater
 * \brief	Computes barcode templates (using the mathematics of "vineyard updates" to store in the Arrangement)
 * \author	Matthew L. Wright
 * \date	March 2015
 */

#ifndef __PERSISTENCE_UPDATER_H__
#define __PERSISTENCE_UPDATER_H__

//forward declarations
class Face;
class Halfedge;
class IndexMatrix;
class MapMatrix_Perm;
class MapMatrix_RowPriority_Perm;
class Arrangement;
class MultiBetti;
class FIRep;
class TemplatePoint;
struct TemplatePointsMatrixEntry;
//struct Subset_map;

#include "template_points_matrix.h"

#include <interface/progress.h>
#include <map>
#include <vector>
#include <boost/unordered_map.hpp>
#include "subset.h"
#include "time_root.h"
#include "computing_s.h"
#include "debug.h"
#include "math/bifiltration_data.h"


class PersistenceUpdater {
public:
    PersistenceUpdater(Arrangement& m, FIRep& b, BifiltrationData& b_data, std::vector<TemplatePoint>& xi_pts, unsigned verbosity); //constructor for when we must compute all of the barcode templates
    //PersistenceUpdater(Arrangement& m, std::vector<TemplatePoint>& xi_pts); //constructor for when we load the pre-computed barcode templates from a RIVET data file

    //functions to compute and store barcode templates in each 2-cell of the arrangement
    void store_barcodes_with_reset(std::vector<std::shared_ptr<Halfedge>>& path, Progress& progress); //hybrid approach -- for expensive crossings, resets the matrices and does a standard persistence calculation
    void store_barcodes_quicksort(std::vector<std::shared_ptr<Halfedge>>& path); ///TODO -- for expensive crossings, rearranges columns via quicksort and fixes the RU-decomposition globally

    typedef std::vector<std::shared_ptr<TemplatePointsMatrixEntry>> zeta_seq;
    typedef std::pair<std::pair<unsigned,unsigned>,int> bigrade_root;
    typedef std::pair<unsigned,unsigned> bigrade;

    void print_zeta_seq(zeta_seq& zs);

    void print_iterator(std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it);

    void print_oracle_image(std::map<int, Subset_map>& UF);

    void print_true_bigrade(std::pair<int,int> big);

    void print_T_next(boost::unordered::unordered_map<bigrade,
                      boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                      bool with_comp_labels = false);

    void recursively_print_T_cur(Time_root& T_cur);

    static bool zsComparator(const std::shared_ptr<TemplatePointsMatrixEntry> e,
                             const std::shared_ptr<TemplatePointsMatrixEntry> f);

    /*static bool bigradeComparator(const std::pair<double,double> b1,
                                  const std::pair<double,double> b2);*/

    bool bigrade_geq(const std::pair<double,double> b1,
                     const std::pair<double,double> b2)
    {
        return (b1.first >= b2.first && b1.second >= b2.second);
    }

    typedef std::pair<double,double> double_bigrade;

    //template <class ForwardIterator>
    //ForwardIterator bigrade_lower_bound (ForwardIterator first, ForwardIterator last, const double_bigrade& val)
    std::vector<double_bigrade>::iterator bigrade_lower_bound(std::vector<double_bigrade>::iterator first,
                                                              std::vector<double_bigrade>::iterator last,
                                                              const double_bigrade& val)

    {
        std::vector<double_bigrade>::iterator it;
        std::vector<double_bigrade>::iterator last_copy = last;
        std::iterator_traits<std::vector<double_bigrade>::iterator>::difference_type count, step;
        count = distance(first,last);
        while (count>0)
        {
            it = last; step=count/2; std::advance (it,-1*step);
            if (bigrade_geq(*it,val))
            {
                last=--it;
                count-=step+1;
                //debug() << "last = (" << last->first << "," << last->second << ")";
            }
            else count=step;
        }
        if (bigrade_geq(*last,val))
            return last;
        else // found element is < val
        {
            if (last == last_copy) // all elements in list are < val
                return last;
            else                     // there exists an element >= val
                return last+1;
        }

    }

    void draw_dendrogram(boost::unordered::unordered_map<std::pair<double,int>, Time_root>& time_root_to_tr,
                                             Time_root& last_upper_tr);

    // functions for testing templates vs. naive algorithm
    std::string int_set_to_string(std::set<int> comp_label);
    std::set<int> populate_ordered_component_label(Time_root& tr);
    bool is_dendrogram_isomorphic(Time_root& tr1, Time_root& tr2);
    Time_root compute_dendrogram_template(zeta_seq& zs);
    void relabel_dendrogram_with_oracle(Time_root& tr);






    //functions to compute and store dendrogram templates in each 2-cell of the arrangement
    void store_dendrogram_templates(std::vector<std::shared_ptr<Halfedge>>& path, Progress& progress);

    //stores a dendrogram template in a 2-cell of the arrangement
    void store_initial_dendrogram_template(std::shared_ptr<Face> cell, zeta_seq& zs);

    // instantiate a level (at the given bigrade) of T_next
    void instantiate_level(std::map<int, Subset_map>& oracle,
                           boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                           std::pair<unsigned,unsigned> big);

    // populate component labels of a level (at the given bigrade) in T_next
    void populate_component_labels(std::map<int, Subset_map>& oracle,
                                   boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                                   bigrade big);

    // populate num_leaves of all nodes in the dendrogram rooted at tr
    int populate_num_leaves_and_birth_label_str(Time_root& tr);

    // connect different levels (big_1 and big_2) of T_next
    void connect_levels(std::map<int, Subset_map>& oracle,
                        boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                        bigrade big_1,
                        bigrade big_2);

    // copy upper part of T_cur
    void copy_upper_dendrogram(Time_root& T_cur,
                               boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                               bigrade big,
                               std::vector<bigrade_root>& parents_of_iplus1,
                               std::map<int, Subset_map>& oracle,
                               std::shared_ptr<Time_root> parent_ptr);

    // convert component labels to birth_labels
    void convert_component_to_birth_labels(
            boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
            std::vector<bigrade> functor_bigrades);

    // removes nodes with in and out degree 1 and null birth label
    void remove_redundant_nodes(
            std::shared_ptr<Time_root> tr_ptr,
            bigrade lower_cutoff,
            boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next);

    //classify nodes of dendrogram during updates
    void classify_nodes(std::vector<Time_root>& L,
                        std::vector<Time_root>& i_minus_1,
                        std::vector<Time_root>& i_plus_1,
                        std::vector<Time_root>& U,
                        Time_root& tr,
                        std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_prime,
                        bool it_i_minus_1_exists,
                        bool it_i_plus_1_exists);

    //get bigrade from TemplatePointsMatrixEntry
    std::pair<unsigned,unsigned> get_bigrade(std::shared_ptr<TemplatePointsMatrixEntry> u);

    //set oracle
    void set_oracle(boost::unordered::unordered_map<std::pair<int,int>,std::map<int, Subset_map>>& oracle);
    void set_cs(computing_s* cs);

    //function to set all edge weights to uniform value for dendrogram testing
    void set_uniform_anchor_weights(std::vector<std::shared_ptr<Halfedge>>& path);

    //function to set the "edge weights" for each anchor line
    void set_anchor_weights(std::vector<std::shared_ptr<Halfedge>>& path);

    //function to clear the levelset lists -- e.g., following the edge-weight calculation
    void clear_levelsets();

private:
    //data structures

    Arrangement& arrangement; //pointer to the DCEL arrangement in which the barcodes will be stored
    FIRep& bifiltration; //pointer to the bifiltration
    //SimplexTree& bifiltration; //pointer to the bifiltration
    BifiltrationData& bifiltration_data;
    boost::unordered::unordered_map<std::pair<int,int>,std::map<int, Subset_map>> oracle;
    computing_s* cs;
    int dim; //dimension of homology to be computed

    unsigned verbosity;

    TemplatePointsMatrix template_points_matrix; //sparse matrix to hold xi support points -- used for finding anchors (to build the arrangement) and tracking simplices during the vineyard updates (when computing barcodes to store in the arrangement)

    std::map<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>> lift_low; //map from "low" columns to xiMatrixEntrys
    std::map<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>> lift_high; //map from "high" columns to xiMatrixEntrys

    MapMatrix_Perm* R_low; //boundary matrix for "low" simplices
    MapMatrix_Perm* R_high; //boundary matrix for "high" simplices
    MapMatrix_RowPriority_Perm* U_low; //upper-trianglular matrix that records the reductions for R_low
    MapMatrix_RowPriority_Perm* U_high; //upper-trianglular matrix that records the reductions for R_high

    ///TODO: is there a way to avoid maintaining the following permutation vectors?
    std::vector<unsigned> perm_low; //map from column index at initial cell to column index at current cell
    std::vector<unsigned> inv_perm_low; //inverse of the previous map
    std::vector<unsigned> perm_high; //map from column index at initial cell to column index at current cell
    std::vector<unsigned> inv_perm_high; //inverse of the previous map

    ///TESTING ONLY
    //bool testing;
    //MapMatrix_Perm* D_low;
    //MapMatrix_Perm* D_high;

    //functions

    typedef std::vector<unsigned> Perm; //for storing permutations

    //stores multigrade info for the persistence computations (data structures prepared with respect to a near-vertical line positioned to the right of all \xi support points)
    //  low is true for simplices of dimension hom_dim, false for simplices of dimension hom_dim+1
    void store_multigrades(IndexMatrix* ind, bool low);

    //finds the proper order of simplexes for the persistence calculation (with respect to a near-vertical line positioned to the right of all \xi support points)
    //  low is true for simplices of dimension hom_dim, false for simplices of dimension hom_dim+1
    //  simplex_order will be filled with a map : dim_index --> order_index for simplices of the given dimension
    //  NOTE: If a simplex with dim_index i does not appear in the order (i.e. its grade is not less than the LUB of all xi support points), then simplex_order[i] = -1.
    //  returns the number of simplices in the order
    unsigned build_simplex_order(IndexMatrix* ind, bool low, std::vector<int>& simplex_order);

    //counts the number of transpositions that will happen if we cross an anchor and do vineyeard-updates
    unsigned long count_transpositions(std::shared_ptr<TemplatePointsMatrixEntry> at_anchor, bool from_below);

    //counts the number of transpositions that result from separations; used in the above function
    void count_transpositions_from_separations(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser, bool horiz, bool low, unsigned long& count_trans, unsigned& count_lesser);

    //moves grades associated with TemplatePointsMatrixEntry greater, that come before TemplatePointsMatrixEntry lesser in R^2, so that they become associated with lesser
    //   horiz is true iff greater and lesser are on the same horizontal line (i.e., they share the same y-coordinate)
    //   returns a count of the number transpositions performed
    unsigned long split_grade_lists(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser, bool horiz);

    //splits grade lists and updates the permutation vectors, but does NOT do vineyard updates
    void split_grade_lists_no_vineyards(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser, bool horiz);

    //moves all grades associated with TemplatePointsMatrixEntry lesser so that they become associated with TemplatePointsMatrixEntry greater
    void merge_grade_lists(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser);

    //moves columns from an equivalence class given by std::shared_ptr<TemplatePointsMatrixEntry> first to their new positions after or among the columns in the equivalence class given by std::shared_ptr<TemplatePointsMatrixEntry> second
    //  the boolean argument indicates whether an anchor is being crossed from below (or from above)
    //  returns a count of the number of transpositions performed
    unsigned long move_columns(std::shared_ptr<TemplatePointsMatrixEntry> first, std::shared_ptr<TemplatePointsMatrixEntry> second, bool from_below);

    //moves a block of n columns, the rightmost of which is column s, to a new position following column t (NOTE: assumes s <= t)
    //  returns a count of the number of transpositions performed
    unsigned long move_low_columns(int s, unsigned n, int t);
    unsigned long move_high_columns(int s, unsigned n, int t);

    //performs a vineyard update corresponding to the transposition of columns a and (a + 1)
    void vineyard_update_low(unsigned a);
    void vineyard_update_high(unsigned a);

    //swaps two blocks of columns by updating the total order on columns, then rebuilding the matrices and computing a new RU-decomposition
    void update_order_and_reset_matrices(std::shared_ptr<TemplatePointsMatrixEntry> first, std::shared_ptr<TemplatePointsMatrixEntry> second, bool from_below, MapMatrix_Perm* RL_initial, MapMatrix_Perm* RH_initial);

    //updates the total order on columns, rebuilds the matrices, and computing a new RU-decomposition for a NON-STRICT anchor
    void update_order_and_reset_matrices(MapMatrix_Perm* RL_initial, MapMatrix_Perm* RH_initial);

    //swaps two blocks of simplices in the total order, and counts switches and separations
    void count_switches_and_separations(std::shared_ptr<TemplatePointsMatrixEntry> at_anchor, bool from_below, unsigned long& switches, unsigned long& seps);

    //used by the previous function to split grade lists at each anchor crossing
    void do_separations(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser, bool horiz);

    //removes entries corresponding to an TemplatePointsMatrixEntry from lift_low and lift_high
    void remove_lift_entries(std::shared_ptr<TemplatePointsMatrixEntry> entry);

    //creates the appropriate entries in lift_low and lift_high for an TemplatePointsMatrixEntry with nonempty sets of "low" or "high" simplices
    void add_lift_entries(std::shared_ptr<TemplatePointsMatrixEntry> entry);

    //stores a barcode template in a 2-cell of the arrangement
    ///TODO: IMPROVE THIS -- track most recent barcode at the simplicial level and re-examine only the necessary columns!!!
    void store_barcode_template(std::shared_ptr<Face> cell);

    //chooses an initial threshold by timing vineyard updates corresponding to random transpositions
    void choose_initial_threshold(unsigned decomp_time, unsigned long & num_trans, unsigned & trans_time, unsigned long & threshold);

    ///TESTING ONLY
    //void check_low_matrix(MapMatrix_Perm* RL, MapMatrix_RowPriority_Perm* UL);
    //void check_high_matrix(MapMatrix_Perm* RH, MapMatrix_RowPriority_Perm* UH);
    void print_perms(Perm& per, Perm& inv);
    void print_high_partition();
};

#endif // __PERSISTENCE_UPDATER_H__
