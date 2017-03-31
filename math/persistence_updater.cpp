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

#include "persistence_updater.h"

#include "../dcel/anchor.h"
#include "../dcel/barcode_template.h"
#include "../dcel/dcel.h"
#include "dcel/arrangement.h"
#include "debug.h"
#include "index_matrix.h"
#include "map_matrix.h"
#include "multi_betti.h"
#include "firep.h"
//#include "graph.h"
#include "dendrogram_data.h"

#include <chrono>
#include <stdexcept> //for error-checking and debugging
#include <stdlib.h> //for rand()
#include <timer.h>

#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
//#include <boost/filesystem/path.hpp>
//#include <boost/filesystem/operations.hpp>
//#include <boost/filesystem/exception.hpp>
#include <boost/graph/graphviz.hpp>
#include "dendrogram_viz.h"
#include <queue>
/*
//constructor for when we must compute all of the barcode templates
PersistenceUpdater::PersistenceUpdater(Arrangement& m, FIRep& b, std::vector<TemplatePoint>& xi_pts, unsigned verbosity)
    : arrangement(m)
    , bifiltration(b)
    , dim(b.hom_dim)
    , verbosity(verbosity)
    , template_points_matrix(m.x_exact.size(), m.y_exact.size())

//    , testing(false)
{
    //fill the xiSupportMatrix with the xi support points and anchors
    //  also stores the anchors in xi_pts
    std::cout << "========================= COMPUTE template_points_matrix: " << std::endl;
    for (auto& matrix_entry : template_points_matrix.fill_and_find_anchors(xi_pts))
    {
        //Add anchors to arrangement also
        //std::cout << "(" << matrix_entry->x << "," << matrix_entry->y << ")";
        m.add_anchor(Anchor(matrix_entry));
    }
    std::cout << std::endl;
}*/

//constructor for multi-critical dendrogram
PersistenceUpdater::PersistenceUpdater(Arrangement& m, FIRep& b, BifiltrationData& b_data, std::vector<TemplatePoint>& xi_pts, unsigned verbosity)
    : arrangement(m)
    , bifiltration(b)
    , bifiltration_data(b_data)
    , dim(b.hom_dim)
    , verbosity(verbosity)
    , template_points_matrix(m.x_exact.size(), m.y_exact.size())

//    , testing(false)
{
    //fill the xiSupportMatrix with the xi support points and anchors
    //  also stores the anchors in xi_pts
    std::cout << "========================= COMPUTE template_points_matrix: " << std::endl;
    for (auto& matrix_entry : template_points_matrix.fill_and_find_anchors(xi_pts))
    {
        //Add anchors to arrangement also
        //std::cout << "(" << matrix_entry->x << "," << matrix_entry->y << ")";
        m.add_anchor(Anchor(matrix_entry));
    }
    std::cout << std::endl;
}

////constructor for when we load the pre-computed barcode templates from a RIVET data file
//PersistenceUpdater::PersistenceUpdater(Arrangement& m, std::vector<TemplatePoint>& xi_pts) :
//    arrangement(m),
//    template_points_matrix(m.x_grades.size(), m.y_grades.size()),
//    testing(false)
//{
//    //fill the xiSupportMatrix with the xi support points
//    template_points_matrix.fill_and_find_anchors(xi_pts, m);
//}

//typedef boost::GraphvizGraph Graph;
//typedef boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;

/*
void PersistenceUpdater::draw_dendrogram(boost::unordered::unordered_map<std::pair<double,int>, Time_root>& time_root_to_tr,
                                         Time_root& last_upper_tr)
{

}*/

enum {is_forest_connector = -2}; // MAGIC NUMBER

void PersistenceUpdater::print_zeta_seq(zeta_seq& zs)
{
    for (auto entry : zs)
    {
        std::cout << "(" << entry->x << "," << entry->y << ")";
    }
    std::cout << std::endl;
}

void PersistenceUpdater::print_iterator(std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it)
{
    std::cout << "(x,y) = (" << (*it)->x << "," << (*it)->y << ")" << std::endl;
}

void PersistenceUpdater::print_oracle_image(std::map<int, Subset_map>& UF)
{
    for (auto v_subset : UF)
    {
        std::cout << "v = " << v_subset.first << std::endl;
        std::cout << "v.parent = " << v_subset.second.parent << std::endl;
        std::cout << "v.rank = " << v_subset.second.rank << std::endl;
    }
}

void PersistenceUpdater::print_true_bigrade(std::pair<int,int> big)
{
    std::cout << "(" << arrangement.x_grades[big.first] << ","
              << arrangement.y_grades[big.second] << ")" << std::endl;
}

void PersistenceUpdater::print_T_next(boost::unordered::unordered_map<bigrade,
                                      boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                                      bool with_comp_labels)
{
    for(auto big_r_tr : T_next)
    {
        //std::pair<unsigned,unsigned> big = big_r_tr.first;
        for (auto r_tr : big_r_tr.second)
        {
            std::cout << "===============" << std::endl;
            Time_root::print_tr(*r_tr.second, with_comp_labels);
            /*if (with_comp_labels)
            {
                std::stringstream ss_cl;
                ss_cl << "(";
                for (auto v : r_tr.second.component_label)
                    ss_cl << v << ",";
                ss_cl << ")";
                std::cout << "cl = " << ss_cl.str() << std::endl;
            }*/
            std::cout << "children: " << std::endl;
            Time_root::print_children((*r_tr.second).children, with_comp_labels);
        }
    }
}

/*
bool PersistenceUpdater::bigradeComparator(const std::pair<double,double> b1,
                                           const std::pair<double,double> b2)
{
    return (b1.first < b2.first) || ((b1.first == b2.first) && (b1.second < b2.second));
}*/


bool PersistenceUpdater::zsComparator(const std::shared_ptr<TemplatePointsMatrixEntry> e,
                                      const std::shared_ptr<TemplatePointsMatrixEntry> f)
{
    return (e->x < f->x) || (e->x == f->x && e->y < f->y);
}

bool PersistenceUpdater::is_dendrogram_isomorphic(Time_root& tr1, Time_root& tr2)
{
    //std::string comp_label_1 = int_set_to_string(tr1.ordered_component_label);
    //std::string comp_label_2 = int_set_to_string(tr2.ordered_component_label);
    //bool comp_labels_match = (comp_label_1 == comp_label_2);
    //std::cout << "comp_label_1 = " << comp_label_1 << std::endl;
    //std::cout << "comp_label_2 = " << comp_label_2 << std::endl;
    std::string cl1 = int_set_to_string(tr1.ordered_component_label);
    std::string cl2 = int_set_to_string(tr2.ordered_component_label);

    bool is_isomorphic = (tr1.bigrade == tr2.bigrade) &&
                         (cl1 == cl2) &&
                         (tr1.num_leaves == tr2.num_leaves);

    std::cout << "is_isomorphic = " << is_isomorphic << std::endl;
    std::cout << "tr1.bigrade = (" << tr1.bigrade.first << "," << tr1.bigrade.second << ")" << std::endl;
    std::cout << "tr2.bigrade = (" << tr2.bigrade.first << "," << tr2.bigrade.second << ")" << std::endl;
    std::cout << "tr1.ordered_comp_label = " << cl1 << std::endl;
    std::cout << "tr2.ordered_comp_label = " << cl2 << std::endl;
    std::cout << "tr1.num_leaves = " << tr1.num_leaves << std::endl;
    std::cout << "tr2.num_leaves = " << tr2.num_leaves << std::endl;

    for (auto child_ptr_1 : tr1.children)
    {
        Time_root& child_1 = *child_ptr_1;
        std::string comp_label_1 = int_set_to_string(child_1.ordered_component_label);
        bool match_found = false;
        for (auto child_ptr_2 : tr2.children)
        {
            Time_root& child_2 = *child_ptr_2;
            std::string comp_label_2 = int_set_to_string(child_2.ordered_component_label);
            if (comp_label_1 == comp_label_2)
            {
                match_found = true;
                is_isomorphic = (is_isomorphic && is_dendrogram_isomorphic(child_1, child_2));
            }
        }
        if (!match_found)
            throw std::runtime_error("match not found");
    }
    return is_isomorphic;
}

std::set<int> PersistenceUpdater::populate_ordered_component_label(Time_root& tr)
{
    for (int v : tr.birth_label)
        tr.ordered_component_label.insert(v);
    for (auto tr_ptr : tr.children)
    {
        Time_root& child = *tr_ptr;
        std::set<int> child_comp_label = populate_ordered_component_label(child);
        for (int v : child_comp_label)
            tr.ordered_component_label.insert(v);
    }
    return tr.ordered_component_label;
}

std::string PersistenceUpdater::int_set_to_string(std::set<int> comp_label)
{
    std::stringstream ss;
    for (int v : comp_label)
    {
        ss << v << ",";
    }
    return ss.str();
}



Time_root PersistenceUpdater::compute_dendrogram_template(zeta_seq& zs)
{
    boost::unordered::unordered_map<double, double_bigrade> appearance_to_bigrade; // map from time of appearance to point in zeta seq
    boost::unordered::unordered_map<double_bigrade, double> bigrade_to_appearance; // inverse map
    std::vector<double_bigrade> zs_bigrades;

    int j = 0;
    for(std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it = zs.begin(); it != zs.end(); ++it)
    {
        auto entry = *it;
        double_bigrade big = double_bigrade(entry->x,entry->y);
        zs_bigrades.push_back(big);
        appearance_to_bigrade[j] = big;
        bigrade_to_appearance[big] = j;
        j += 1;
    }

    // 1-critical version SimplexSet vertices = bifiltration.get_ordered_simplices();
    // 1-critical version SimplexSet edges = bifiltration.get_ordered_high_simplices();
    boost::unordered::unordered_map<unsigned, double> vertex_appearance;
    boost::unordered::unordered_map<unsigned, std::pair<double,double>> vertex_bigrade;

    SimplexInfo* vertices_to_grades = bifiltration_data.getSimplices(0);
    SimplexInfo* edges_to_grades = bifiltration_data.getSimplices(1);

    std::cout << "vertexVec_grade info: " << std::endl;
    for (auto vertexVec_grade : *vertices_to_grades)
    {
        unsigned v = vertexVec_grade.first[0];
        AppearanceGrades grades = vertexVec_grade.second;
        for (Grade g : grades)
        {
            double_bigrade big = double_bigrade(g.x,g.y);
            std::vector<double_bigrade>::iterator it_lowerbound;
            it_lowerbound = bigrade_lower_bound(zs_bigrades.begin(),zs_bigrades.end(),big);
            double_bigrade lowerbound = *it_lowerbound;
            double birth = bigrade_to_appearance[lowerbound];


            std::cout << "(v,x,y,birth) = (" <<
                    v << " , " <<
                    g.x << " , " <<
                    g.y << " , " <<
                    birth << ")" << std::endl;

            if (vertex_appearance.find(v) == vertex_appearance.end())
            {
                vertex_appearance[v] = birth;
                vertex_bigrade[v] = lowerbound;
            }
            else
            {
                if (birth < vertex_appearance[v])
                {
                    vertex_appearance[v] = birth;
                    vertex_bigrade[v] = lowerbound;
                }
            }
        }
    }

    /* 1-critical version
    for (auto st_node : vertices)
    {
        unsigned v = st_node->get_vertex();
        unsigned x = st_node->grade_x();
        unsigned y = st_node->grade_y();
        //double birth = y;
        double_bigrade big = double_bigrade(x,y);
        std::vector<double_bigrade>::iterator it_lowerbound;
        it_lowerbound = bigrade_lower_bound(zs_bigrades.begin(),zs_bigrades.end(),big);
        double_bigrade lowerbound = *it_lowerbound;
        double birth = bigrade_to_appearance[lowerbound];
        /*it_lowerbound = std::lower_bound(zs_bigrades.begin(),zs_bigrades.end(),big,bigradeComparator);
        double_bigrade bound = *it_lowerbound;
        if (bound != big)
            bound = std::upper_bound(zs_bigrades.begin(),zs_bigrades.end(),big,bigradeComparator);
        double birth = bigrade_to_appearance[bound];*/

        /*
        debug() << "lowerbound = (" << lowerbound.first << "," << lowerbound.second << ")";

        debug() << "(v,x,y,birth) = (" <<
                v << " , " <<
                x << " , " <<
                y << " , " <<
                birth << ")";
        if (vertex_appearance.find(v) == vertex_appearance.end())
        {
            vertex_appearance[v] = birth;
            vertex_bigrade[v] = lowerbound;
        }
        else
        {
            if (birth < vertex_appearance[v])
            {
                vertex_appearance[v] = birth;
                vertex_bigrade[v] = lowerbound;
            }
        }
    }*/

    boost::unordered::unordered_map<std::pair<unsigned,unsigned>, double> edge_appearance;
    boost::unordered::unordered_map<std::pair<unsigned,unsigned>, std::pair<double,double>> edge_bigrade;

    for (auto edgeVec_grade : *edges_to_grades)
    {
        std::vector<int> edge_vector = edgeVec_grade.first;
        AppearanceGrades grades = edgeVec_grade.second;
        std::pair<unsigned,unsigned> edge (edge_vector[0], edge_vector[1]);
        for (Grade g : grades)
        {
            double_bigrade big = double_bigrade(g.x,g.y);
            std::vector<double_bigrade>::iterator it_lowerbound;
            it_lowerbound = bigrade_lower_bound(zs_bigrades.begin(),zs_bigrades.end(),big);
            if (it_lowerbound == zs_bigrades.end()) //case where edge is born after last element in zs
                continue;
            double_bigrade lowerbound = *it_lowerbound;
            double birth = bigrade_to_appearance[lowerbound];

            if (edge_appearance.find(edge) == edge_appearance.end())
            {
                edge_appearance[edge] = birth;
                edge_bigrade[edge] = lowerbound;
            }
            else
            {
                if (birth < edge_appearance[edge])
                {
                    edge_appearance[edge] = birth;
                    edge_bigrade[edge] = lowerbound;
                }
            }
        }
    }

    /* 1-critical version
    for (auto st_node : edges)
    {
        int global_index = st_node->global_index();
        std::vector<int> edge_vector = bifiltration.find_vertices(global_index);
        std::pair<unsigned,unsigned> edge (edge_vector[0], edge_vector[1]);
        unsigned x = st_node->grade_x();
        unsigned y = st_node->grade_y();
        //double birth = y;
        double_bigrade big = double_bigrade(x,y);
        std::vector<double_bigrade>::iterator it_lowerbound;
        it_lowerbound = bigrade_lower_bound(zs_bigrades.begin(),zs_bigrades.end(),big);
        if (it_lowerbound == zs_bigrades.end()) //case where edge is born after last element in zs
            continue;
        double_bigrade lowerbound = *it_lowerbound;
        double birth = bigrade_to_appearance[lowerbound];

        std::unordered_set<std::pair<int,int>> T_1 = cs->get_T_1();
        std::pair<int,int> big_int = std::pair<int,int>((int) big.first,(int) big.second);

        // case where edge does not live in T_1 \subset S
        /*if (T_1.find(big_int) == T_1.end())
        {
            continue;
        }

        if (edge_appearance.find(edge) == edge_appearance.end())
        {
            edge_appearance[edge] = birth;
            edge_bigrade[edge] = lowerbound;
        }
        else
        {
            if (birth < edge_appearance[edge])
            {
                edge_appearance[edge] = birth;
                edge_bigrade[edge] = lowerbound;
            }
        }
    }*/

    DenseGRAPH<Edge>* graph = new DenseGRAPH<Edge>(vertex_appearance.size(), true);
    //debug() << "print edges of graph";
    for (auto edge_wt : edge_appearance)
    {
        std::pair<unsigned,unsigned> e = edge_wt.first;
        unsigned v = e.first;
        unsigned w = e.second;
        graph-> insert(EdgePtr( new Edge(v,w,edge_wt.second,edge_bigrade[e])));
        /*debug() << "(v,w,edge_wt,edge_bigrade) = (" <<
                   v << " , " <<
                   w << " , " <<
                   edge_wt.second << " , (" <<
                   edge_bigrade[e].first << " , " << edge_bigrade[e].second << "))";
        */
    }

    std::vector<EdgePtr> mst;
    boost::unordered::unordered_map<std::pair<double,int>, Time_root> time_root_to_tr;

    //debug() << "START compute dendrogram template";
    std::vector<Time_root> upper_trs;
    Time_root tr;
    Dendrogram_data::compute_dendrogram(graph, mst, time_root_to_tr, tr, vertex_appearance, vertex_bigrade, appearance_to_bigrade, upper_trs);
    /*debug() << "upper_trs: ";
    Time_root::print_trs(upper_trs, true);
    for ( std::vector<EdgePtr>::const_iterator it = mst.begin(); it != mst.end(); ++it )
    {

        debug() << it->get()->v << " "
                << it->get()->w << " "
                << it->get()->wt << "\n";
    }*/
    if (upper_trs.size() > 1)
    {
        Time_root forest_connector;
        forest_connector.r = is_forest_connector;
        forest_connector.bigrade = bigrade(UINT_MAX,UINT_MAX);
        for (Time_root& upper_tr : upper_trs)
        {
            forest_connector.children.push_back(std::make_shared<Time_root>(upper_tr));
            forest_connector.num_leaves += upper_tr.num_leaves;
        }
        tr = forest_connector;
        std::cout << "forest_connector : " << std::endl;
        Time_root::print_tr(forest_connector);
        Time_root::print_children(forest_connector.children);
    }
    else // upper_trs.size() == 1
        tr = upper_trs.front();

    //debug() << "END compute dendrogram template";

    for (auto t_big : appearance_to_bigrade)
    {
        auto t = t_big.first;
        auto big = t_big.second;
        std::cout << "t = " << t << " , " << "big = ("
                  << big.first << "," << big.second << ")"
                  << std::endl;
    }
    return tr;
}

void PersistenceUpdater::relabel_dendrogram_with_oracle(Time_root& tr)
{
    if (tr.r != is_forest_connector)
    {
        int r = tr.r;
        auto ora = oracle[tr.bigrade];
        int h_r = Find(ora, r);
        tr.r = h_r;
    }
    for (auto child_ptr : tr.children)
        relabel_dendrogram_with_oracle(*child_ptr);
    return;
}

void PersistenceUpdater::store_dendrogram_templates(std::vector<std::shared_ptr<Halfedge>>& path, Progress& progress)
{
    // get underlying buffer
    //std::streambuf* orig_buf = std::cout.rdbuf();

    // set null
    //std::cout.rdbuf(NULL);

    // PART 1: INITIAL CELL COMPUTATION
    debug() << "PART 1: INITIAL CELL COMPUTATION";
    int zeta_seq_length = template_points_matrix.height();
    debug() << "zeta_seq_length = " << zeta_seq_length;
    zeta_seq cur_zs;

    debug() << "initial zs";
    for (int j = 0; j < zeta_seq_length; j++)
    {
        debug() << "j = " << j << " , ";
        if (template_points_matrix.get_row(j) != NULL)
        {
            const auto& entry = template_points_matrix.get_row(j);
            debug() << "(" << entry->x << "," << entry->y << ")";
            cur_zs.push_back(entry);
        }
    }
    debug() << " ========================= initial zs:";
    print_zeta_seq(cur_zs);

    //store the dendrogram template in the first cell
    debug() << "store the dendrogram template in the first cell";
    std::shared_ptr<Face> first_cell = arrangement.topleft->get_twin()->get_face();
    store_initial_dendrogram_template(first_cell, cur_zs);
    Time_root& tr_first = first_cell->get_dendrogram();
    std::cout << "tr_first before relabeling: " << std::endl;
    Time_root::recursively_print_tr(tr_first);
    //recursively_print_T_cur(tr_first);

    // relabel tr_first so that its .r labels correspond to oracle
    relabel_dendrogram_with_oracle(tr_first);
    std::cout << "tr_first after relabeling: " << std::endl;
    Time_root::recursively_print_tr(tr_first);

    /*Time_root T_cur = tr_first;

    std::cout << "extend dendrogram to last bigrade in cur_zs" << std::endl;
    if (get_bigrade(cur_zs.back()) != tr_first.bigrade)
    {
        T_cur.bigrade = get_bigrade(cur_zs.back());
        //T_cur.children.push_back(std::make_shared<Time_root>(tr_first));
    }*/
    Time_root T_cur = tr_first;
    //bigrade extended_bigrade = tr_first.bigrade; // used in step 4 + 5
    //tr_first.bigrade = get_bigrade(cur_zs.back()); //extend dendrogram template to last bigrade in cur_zs
    //Time_root::print_tr(tr_first);
    Dendrogram_viz::write_dendrogram_dot_file(tr_first, "../../../initial.dot", arrangement.x_grades, arrangement.y_grades);

    debug() << "first_cell = " << arrangement.FID(first_cell);

    // Time_root T_cur = tr_first;

    //truncate cur_zs to the minimum length without changing dendrogram
    //bigrade zs_truncation_pt = T_cur.bigrade; // we can truncate all bigrades past this point in cur_zs since they do not affect the dendrogram
    //std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_truncation_pt;
    /*for (auto entry : cur_zs)
    {
        if (get_bigrade(entry) == zs_truncation_pt)
    }*/
    /*int n = 0;
    for (auto it = cur_zs.begin(); it != cur_zs.end(); ++it)
    {
        n += 1;
        if (get_bigrade(*it) == zs_truncation_pt)
            break;
    }*/
    std::cout << " ========================= initial zs:";
    print_zeta_seq(cur_zs);

    /*cur_zs.resize(n);

    std::cout << " ========================= truncated initial zs:";
    print_zeta_seq(cur_zs);*/

    boost::unordered::unordered_map<std::shared_ptr<Face>, zeta_seq> face_to_zs;

    // PART 2: TRAVERSE PATH AND PERFORM DENDROGRAM UPDATES
    debug() << "PART 2: TRAVERSE PATH AND PERFORM DENDROGRAM UPDATES";



    for (unsigned j = 0; j < path.size(); j++)
    {
        std::cout << "---------------------------------------------------------" << std::endl;
        std::cout << "step j = " << j << std::endl;
        progress.progress(j); //update progress bar

        //determine which anchor is represented by this edge
        std::shared_ptr<Anchor> cur_anchor = (path[j])->get_anchor();
        std::shared_ptr<TemplatePointsMatrixEntry> at_anchor = cur_anchor->get_entry();

        std::cout << "entry at_anchor = (" << at_anchor->x << "," << at_anchor->y << ")" << std::endl;

        //get equivalence classes for this anchor
        std::shared_ptr<TemplatePointsMatrixEntry> down = at_anchor->down;
        std::shared_ptr<TemplatePointsMatrixEntry> left = at_anchor->left;

        // iterator for the changed zeta, z_i' in paper
        std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_prime;

        // compute next zeta sequence
        zeta_seq next_zs;
        next_zs = cur_zs;

        // Case: strict anchor
        if (down != nullptr && left != nullptr) //then this is a strict anchor and some simplices swap
        {
            std::cout << "CASE : strict anchor" << std::endl;
            std::cout << "down = (" << down->x << "," << down->y << ")" << std::endl;
            std::cout << "left = (" << left->x << "," << left->y << ")" << std::endl;
            if (verbosity >= 6) {
                debug() << "  step " << j << " of path: crossing (strict) anchor at (" << cur_anchor->get_x() << ", " << cur_anchor->get_y() << ") into cell " << arrangement.FID((path[j])->get_face()) << "; edge weight: " << cur_anchor->get_weight();
            }

            if (cur_anchor->is_above()) //then the anchor is crossed from below to above
            {
                std::cout << "anchor is crossed from below to above" << std::endl;
                // binary search for down in cur_zs
                std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_down;
                it_down = std::lower_bound(next_zs.begin(),next_zs.end(),down,zsComparator);
                std::cout << "*it_down before = (" << get_bigrade(*it_down).first << "," << get_bigrade(*it_down).second << ")" << std::endl;
                *it_down = left;
                std::cout << "*it_down after = (" << get_bigrade(*it_down).first << "," << get_bigrade(*it_down).second << ")" << std::endl;
                it_prime = it_down;
            }
            else //then anchor is crossed from above to below
            {
                std::cout << "anchor is crossed from above to below" << std::endl;
                // binary search for left in cur_zs
                std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_left;
                it_left = std::lower_bound(next_zs.begin(),next_zs.end(),left,zsComparator);
                std::cout << "*it_left before = (" << get_bigrade(*it_left).first << "," << get_bigrade(*it_left).second << ")" << std::endl;
                *it_left = down;
                std::cout << "*it_left after = (" << get_bigrade(*it_left).first << "," << get_bigrade(*it_left).second << ")" << std::endl;
                it_prime = it_left;
            }
        }

        // case : weak anchor
        else
        {
            std::cout << "CASE : weak anchor" << std::endl;
            bool neighbor_is_down = true;
            std::shared_ptr<TemplatePointsMatrixEntry> neighbor = at_anchor->down;
            if (neighbor == nullptr)
            {
                std::cout << "neighbor == nullptr" << std::endl;
                neighbor = at_anchor->left;
                neighbor_is_down = false;
            }
            std::cout << "neighbor = (" << neighbor->x << "," << neighbor->y << ")" << std::endl;
            std::cout << "neighbor_is_down = " << neighbor_is_down << std::endl;
            std::cout << "cur_anchor->is_above() = " << cur_anchor->is_above() << std::endl;
            // split
            if ((cur_anchor->is_above() && !neighbor_is_down) || (!cur_anchor->is_above() && neighbor_is_down))
            {
                std::cout << "split case" << std::endl;
                // binary search for neighbor in cur_zs
                std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_anchor;
                it_anchor = std::lower_bound(next_zs.begin(),next_zs.end(),at_anchor,zsComparator);
                std::cout << "*it_anchor = (" << get_bigrade(*it_anchor).first << "," << get_bigrade(*it_anchor).second << ")" << std::endl;
                next_zs.insert(it_anchor,neighbor);
                //std::cout << "========================= next_zs:";
                //print_zeta_seq(next_zs);
                it_prime = std::lower_bound(next_zs.begin(),next_zs.end(),neighbor,zsComparator);
                std::cout << "*it_prime = (" << get_bigrade(*it_prime).first << "," << get_bigrade(*it_prime).second << ")" << std::endl;
                //std::cout << "========================= next_zs:";
                //print_zeta_seq(next_zs);
            }

            // merge
            else
            {
                std::cout << "merge case" << std::endl;
                // binary search for neighbor in cur_zs
                std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_neighbor;
                it_neighbor = std::lower_bound(next_zs.begin(),next_zs.end(),neighbor,zsComparator);
                std::cout << "*it_neighbor = (" << get_bigrade(*it_neighbor).first << "," << get_bigrade(*it_neighbor).second << ")" << std::endl;
                it_prime = std::next(it_neighbor,1);
                std::cout << "*it_prime before = (" << get_bigrade(*it_prime).first << "," << get_bigrade(*it_prime).second << ")" << std::endl;
                next_zs.erase(it_neighbor);
                //std::cout << "========================= next_zs:";
                //print_zeta_seq(next_zs);
                it_prime = std::prev(it_prime,1); // need to increment iterator because we deleted the previous element
                std::cout << "*it_prime after = (" << get_bigrade(*it_prime).first << "," << get_bigrade(*it_prime).second << ")" << std::endl;
            }
        }

        std::cout << "========================= cur_zs:";
        print_zeta_seq(cur_zs);
        std::cout << "========================= next_zs:";
        print_zeta_seq(next_zs);

        // if we have already handled this face, move on to next face
        std::shared_ptr<Face> current_face = (path[j])->get_face();
        if (current_face->has_been_visited())
        {
            std::cout << "current_face->has_been_visited" << std::endl;
            std::cout << "face_to_zs[current_face] = ";
            print_zeta_seq(face_to_zs[current_face]);
            T_cur = current_face->get_dendrogram();
            //std::cout << "T_cur:" << std::endl;
            Time_root::recursively_print_tr(T_cur);
            cur_zs = next_zs;
            cur_anchor->toggle(); // we need to remember we crossed cur_anchor (as we do in step 8)
            continue;
        }


        //boost::unordered::unordered_map<int, Time_root> K;
        //enum V_type {L, i_minus_1, i_plus_1, U};
        //boost::unordered::unordered_map<V_type, Time_root>
        std::vector<Time_root> L;
        std::vector<Time_root> i_minus_1;
        std::vector<Time_root> i_plus_1;
        std::vector<Time_root> U;
        //std::vector<unsigned> R_i_minus_1;
        bool it_i_minus_1_exists = false;
        bool it_i_plus_1_exists = false;

        std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_i_minus_1;
        if (next_zs.begin() != it_prime)
        {
            std::cout << "it_i_minus_1 is defined" << std::endl;
            it_i_minus_1 = std::prev(it_prime,1);
            it_i_minus_1_exists = true;
        }
        std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_i_plus_1;
        if (std::prev(next_zs.end())!= it_prime)
        {
            std::cout << "it_i_plus_1 is defined" << std::endl;
            it_i_plus_1 = std::next(it_prime,1);
            it_i_plus_1_exists = true;
        }

        // step -1
        if (it_i_minus_1_exists)
        {
            if (T_cur.bigrade <= get_bigrade(*it_i_minus_1))
            {
                current_face->mark_as_visited();
                cur_anchor->toggle(); // we need to remember we crossed cur_anchor (as we do in step 8)
                face_to_zs[current_face] = next_zs;
                Time_root& den_template = current_face->get_dendrogram();
                Time_root template_tr = T_cur;
                //populate_num_leaves_and_birth_label_str(den_template);

                bool is_isomorphic = compare_template_to_naive_dendrogram(next_zs, template_tr);
                if (!is_isomorphic)
                    throw std::runtime_error("dendrograms are not isomorphic!");

                den_template = T_cur; // no need to reset T_cur since its same
                cur_zs = next_zs;
                continue;
            }
        }

        // step 0
        // I DONT KNOW WHETHER THE INITIALIZATION OF top_bigrade IS RIGHT
        bigrade top_bigrade; //bigrade of top of dendrogram
        if (it_i_plus_1_exists)
            top_bigrade = get_bigrade(*it_i_plus_1);
        else
            top_bigrade = get_bigrade(*it_prime);
        std::cout << "top_bigrade = (" << top_bigrade.first << "," << top_bigrade.second << ")" <<std::endl;

        int top_v;
        if (T_cur.r == is_forest_connector)
            top_v = is_forest_connector;
        else
        {
            if (!it_i_plus_1_exists)
            {
                auto oracle_prime = oracle[get_bigrade(*it_prime)];
                top_v = Find(oracle_prime, T_cur.r);
            }
            else
            {
                auto oracle_i_plus_1 = oracle[get_bigrade(*it_i_plus_1)];
                top_v = Find(oracle_i_plus_1, T_cur.r);
            }
        }



        // step 1
        // classify all nodes of T_cur by types L,i-1,i+1,U
        std::cout << "========================= STEP 1 " << std::endl;
        std::cout << "T_cur : " << std::endl;
        Time_root::print_tr(T_cur);
        std::cout << "it_prime = " << std::endl;
        print_iterator(it_prime);
        if (it_i_minus_1_exists)
        {
            std::cout << "it_i_minus_1 = " << std::endl;
            print_iterator(it_i_minus_1);
        }

        if (it_i_plus_1_exists)
        {
            std::cout << "it_i_plus_1 = " << std::endl;
            print_iterator(it_i_plus_1);
        }

        if (T_cur.bigrade == get_bigrade(*it_i_plus_1))
            i_plus_1.push_back(T_cur);

        std::cout << "Classify nodes: " << std::endl;
        std::cout << "it_i_minus_1_exists = " << it_i_minus_1_exists << std::endl;
        std::cout << "it_i_plus_1_exists = " << it_i_plus_1_exists << std::endl;

        std::cout << "true bigrade of it_prime:" << std::endl;
        print_true_bigrade(get_bigrade(*it_prime));

        // classify T_cur's top node (could restructure classify nodes to encapsulate this step)
        if (it_i_plus_1_exists)
        {
            std::cout << "true bigrade of it_i_plus_1:" << std::endl;
            print_true_bigrade(get_bigrade(*it_i_plus_1));
            if (T_cur.bigrade == get_bigrade(*it_i_plus_1))
                i_plus_1.push_back(T_cur);
        }
        if (it_i_minus_1_exists)
        {
            //std::cout << "reached classification of i_minus_1 and L" << std::endl;
            //std::cout << "child : " << std::endl;
            //Time_root::print_tr(child);
            std::cout << "true bigrade of it_i_minus_1:" << std::endl;
            print_true_bigrade(get_bigrade(*it_i_minus_1));
            if (T_cur.bigrade == get_bigrade(*it_i_minus_1))
                i_minus_1.push_back(T_cur);
        }
        // classify rest of the nodes in T_cur
        classify_nodes(L, i_minus_1, i_plus_1, U, T_cur, it_prime, it_i_minus_1_exists, it_i_plus_1_exists);
        std::cout << ">>>>>>>>>>>> L = " << std::endl;
        Time_root::print_trs(L, true);
        std::cout << ">>>>>>>>>>>> i_minus_1 = " << std::endl;
        Time_root::print_trs(i_minus_1, true);
        std::cout << ">>>>>>>>>>>> i_plus_1 = " << std::endl;
        Time_root::print_trs(i_plus_1, true);
        std::cout << ">>>>>>>>>>>> U = " << std::endl;
        Time_root::print_trs(U, true);

        if (it_i_minus_1_exists)
        {
            std::map<int, Subset_map>& oracle_i_minus_1 = oracle[get_bigrade(*it_i_minus_1)]; // H[z_i'] in paper
            std::cout << "oracle_i_minus_1:" << std::endl;
            print_oracle_image(oracle_i_minus_1);
        }

        // determine vertex representative of T_cur depending on whether it_i_plus_1 exists or not
        int T_cur_root;
        if (T_cur.r == is_forest_connector)
            T_cur_root = is_forest_connector;
        else
        {
            if (!it_i_plus_1_exists)
            {
                auto oracle_prime = oracle[get_bigrade(*it_prime)];
                T_cur_root = Find(oracle_prime, T_cur.r);
            }
            else
            {
                if (T_cur.bigrade == get_bigrade(*it_i_plus_1))
                {
                    auto oracle_i_plus_1 = oracle[get_bigrade(*it_i_plus_1)];
                    T_cur_root = Find(oracle_i_plus_1, T_cur.r);
                }
                else
                    T_cur_root = T_cur.r;
            }
        }


        std::cout << "T_cur_root = " << T_cur_root << std::endl;
        std::cout << "T_cur.r = " << T_cur.r << std::endl;

        boost::unordered::unordered_map<bigrade,
        boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>> T_next; //K in paper

        // Instantiate nodes in i - 1
        for (Time_root& tr : i_minus_1)
        {
            std::cout << "instantiate nodes in i - 1" << std::endl;
            Time_root tr_copy = tr;
            T_next[tr.bigrade][tr.r] = std::make_shared<Time_root>(tr_copy);
        }

        // Instantiate nodes in L plus their functor extensions
        // Need to copy actual roots in L? instead of using references
        for (Time_root& tr : L)
        {
            //Time_root tr_copy = tr;
            std::cout << "instantiate nodes in L" << std::endl;
            auto ora = oracle[get_bigrade(*it_i_minus_1)];
            int h_r = Find(ora, tr.r);
            Time_root tr_prime;
            //tr_prime.r = tr.r;
            tr_prime.r = h_r;
            tr_prime.bigrade = get_bigrade(*it_i_minus_1);
            tr_prime.children.push_back(std::make_shared<Time_root>(tr));
            //T_next[bigrade_root(tr_prime.bigrade,tr_prime.r)] = tr_prime;
            T_next[tr_prime.bigrade][tr_prime.r] = std::make_shared<Time_root>(tr_prime);
        }

        //std::cout << "Status of T_next at beginning of step 2:" << std::endl;
        //print_T_next(T_next);

        // step 2
        std::map<int, Subset_map>& oracle_it_prime = oracle[get_bigrade(*it_prime)]; // H[z_i'] in paper
        std::cout << "Oracle_it_prime:" << std::endl;
        print_oracle_image(oracle_it_prime);
        instantiate_level(oracle_it_prime, T_next, get_bigrade(*it_prime)); // instantiate K[z_i']
        std::cout << "Status of T_next after instantiate_level for step 2:" << std::endl;
        //print_T_next(T_next, true);
        if (it_i_minus_1_exists)
        {
            std::cout << "connect levels for step 2" << std::endl;
            connect_levels(oracle_it_prime, T_next, get_bigrade(*it_i_minus_1), get_bigrade(*it_prime));
            std::cout << "Status of T_next after connect_levels level for step 2:" << std::endl;
            print_T_next(T_next, true);
        }


        // step 3
        std::cout << "step 3" << std::endl;
        if (it_i_plus_1_exists)
        {
            std::map<int, Subset_map>& oracle_it_i_plus_1 = oracle[get_bigrade(*it_i_plus_1)]; // H[z_{i+1}] in paper
            std::cout << "Oracle_it_i_plus_1:" << std::endl;
            print_oracle_image(oracle_it_i_plus_1);
            instantiate_level(oracle_it_i_plus_1, T_next, get_bigrade(*it_i_plus_1));
            std::cout << "Status of T_next after instantiate_level for step 3:" << std::endl;
            //print_T_next(T_next, true);
            connect_levels(oracle_it_i_plus_1,T_next,get_bigrade(*it_prime),get_bigrade(*it_i_plus_1));
            std::cout << "Status of T_next after connect_levels level for step 3:" << std::endl;
            print_T_next(T_next, true);
        }


        // step 4 + 5
        std::cout << "step 4 + 5" << std::endl;
        std::vector<bigrade_root> parents_of_iplus1; //used in step 7 to remove redundant nodes (ACTUALLY SHOULD MAKE INTO UNORDERED_SET FOR UNIQUENESS)
        if (it_i_plus_1_exists)
        {
            if ((get_bigrade(*it_i_plus_1) < T_cur.bigrade) || (top_v == is_forest_connector)) // we need second condition because forest connectors are not expressed in the oracle so it has not been copied yet
            {
                std::map<int, Subset_map>& oracle_it_i_plus_1 = oracle[get_bigrade(*it_i_plus_1)]; // H[z_{i+1}] in paper
                std::cout << "we must copy_upper_dendrogram"<< std::endl;
                top_bigrade = T_cur.bigrade;
                top_v = T_cur.r;
                //std::map<int, Subset_map>& top_oracle = oracle[top_bigrade]; // H[z_{i+1}] in paper
                //top_v = T_cur.r;
                //top_v = Find(top_oracle, T_cur.r);
                copy_upper_dendrogram(T_cur, T_next, get_bigrade(*it_i_plus_1),parents_of_iplus1, oracle_it_i_plus_1, NULL);


                std::cout << "top_bigrade = (" << top_bigrade.first << "," << top_bigrade.second << ")" <<std::endl;
                std::cout << "top_v = " << top_v << std::endl;
            }
            /*if (get_bigrade(*it_i_plus_1) == extended_bigrade)
            {
                std:: cout << "get_bigrade(*it_i_plus_1) == extended_bigrade" << std::endl;
                std::cout << "there's no need to copy_upper_dendrogram" << std::endl;
                for (auto r_ptr : T_next[extended_bigrade])
                {
                    int r = r_ptr.first;
                    std::cout << "r = " << r << std::endl;
                    auto ptr = r_ptr.second;
                    std::cout << "*ptr = " << std::endl;
                    Time_root::print_tr(*ptr);
                    T_next[T_cur.bigrade][r] = ptr;
                    ptr->bigrade = T_cur.bigrade;
                    T_cur_root = r;
                }
                //T_next.erase(extended_bigrade);
            }
            else
            {
                if ((get_bigrade(*it_i_plus_1) < T_cur.bigrade) || (T_cur_root == is_forest_connector)) // we need second condition because forest connectors are not expressed in the oracle so it has not been copied yet
                {
                    std::map<int, Subset_map>& oracle_it_i_plus_1 = oracle[get_bigrade(*it_i_plus_1)]; // H[z_{i+1}] in paper
                    std::cout << "we must copy_upper_dendrogram"<< std::endl;
                    copy_upper_dendrogram(T_cur, T_next, get_bigrade(*it_i_plus_1),parents_of_iplus1, oracle_it_i_plus_1, NULL);
                }
            }*/
        }
        else
        {
            if (top_v == is_forest_connector)
            {
                std::cout << "T_cur_root == is_forest_connector" << std::endl;
                //std::map<int, Subset_map>& oracle_it_prime = oracle[get_bigrade(*it_prime)];
                //copy_upper_dendrogram(T_cur,T_next,get_bigrade(*it_prime),parents_of_iplus1,oracle_it_prime, NULL);
                Time_root tr;
                tr.r = top_v;
                tr.bigrade = T_cur.bigrade; // does this even matter if its a forest_connector?
                //tr.bigrade = bigrade(next_zs.back()->x,next_zs.back()->y);
                for (auto r_ptr : T_next[get_bigrade(*it_prime)])
                {
                    tr.children.push_back(r_ptr.second);
                }
                T_next[T_cur.bigrade][is_forest_connector] = std::make_shared<Time_root>(tr);
                top_bigrade = T_cur.bigrade;
                //top_v = is_forest_connector;
                std::cout << "top_bigrade = (" << top_bigrade.first << "," << top_bigrade.second << ")" <<std::endl;
                std::cout << "top_v = " << top_v;

                //T_next[tr.bigrade][is_forest_connector] = std::make_shared<Time_root>(tr);

            }
            else
                std::cout << "there's no need to copy_upper_dendrogram" << std::endl;
        }


        std::cout << "Status of T_next after copy_upper_dendrogram for step 4 + 5:" << std::endl;
        print_T_next(T_next,true);

        //std::cout << "about to recursively_print_T_cur" << std::endl;
        //std::cout << "T_cur.bigrade = (" << T_cur.bigrade.first << "," << T_cur.bigrade.second << ")" << std::endl;
        //std::cout << "T_cur_root = " << T_cur_root << std::endl;
        //Time_root::print_tr(*T_next[T_cur.bigrade][T_cur.r]);
        //recursively_print_T_cur(*T_next[T_cur.bigrade][T_cur_root]);

        // step 6
        std::cout << "step 6 begin" << std::endl;
        std::cout << "========================= next_zs:";
        print_zeta_seq(next_zs);
        if (it_i_minus_1_exists)
        {
            std::map<int, Subset_map>& oracle_i_minus_1 = oracle[get_bigrade(*it_i_minus_1)]; // H[z_i'] in paper
            std::cout << "oracle_i_minus_1:" << std::endl;
            //print_oracle_image(oracle_i_minus_1);
        }
        std::vector<bigrade> functor_bigrades;

        if (it_i_plus_1_exists)
        {
            functor_bigrades = {get_bigrade(*it_prime),
                                get_bigrade(*it_i_plus_1)};
        }
        else
        {
            functor_bigrades = {get_bigrade(*it_prime)};
        }

        if (it_i_minus_1_exists)
        {
            std::map<int, Subset_map>& oracle_it_i_minus_1 = oracle[get_bigrade(*it_i_minus_1)]; // H[z_{i-1}] in paper
            std::cout << "about to populate_component_labels" << std::endl;
            populate_component_labels(oracle_it_i_minus_1, T_next, get_bigrade(*it_i_minus_1));
            std::cout << "Status of T_next after populate_component_labels for step 6: " << std::endl;
            //print_T_next(T_next, true);
        }
        std::cout << "about to convert_component_to_birth_labels" << std::endl;
        convert_component_to_birth_labels(T_next, functor_bigrades);
        // still need to test convert_component_to_birth_labels with other examples
        std::cout << "Status of T_next after convert_component_to_birth_labels for step 6:" << std::endl;
        //print_T_next(T_next, true);

        // step 7
        std::cout << "parents_of_iplus1:" << std::endl;
        for (bigrade_root big_r : parents_of_iplus1)
        {
            std::cout << "big = (" << big_r.first.first << "," << big_r.first.second
                      << ") , r = " << big_r.second << std::endl;
        }
        for (bigrade_root big_r : parents_of_iplus1)
        {
            //Time_root& tr = *T_next[big_r.first][big_r.second];
            //remove_redundant_nodes(tr,get_bigrade(*it_i_minus_1),T_next);
            auto tr_ptr = T_next[big_r.first][big_r.second];
            if (it_i_minus_1_exists)
                remove_redundant_nodes(tr_ptr,get_bigrade(*it_i_minus_1),T_next);
            else
                remove_redundant_nodes(tr_ptr,get_bigrade(*it_prime),T_next);
        }
        if (parents_of_iplus1.empty())
        {
            std::cout << "parents_of_iplus1.empty()" << std::endl;
            std::cout << "top_bigrade = (" << top_bigrade.first << "," << top_bigrade.second << ")" << std::endl;
            std::cout << "T_cur_root = " << T_cur_root << std::endl;
            std::cout << "top_v = " << top_v;

            //remove_redundant_nodes(*T_next[T_cur.bigrade][T_cur.r],get_bigrade(*it_i_minus_1),T_next);

            if (it_i_minus_1_exists)
            {
                //remove_redundant_nodes(T_next[T_cur.bigrade][T_cur_root], get_bigrade(*it_i_minus_1),T_next);
                //remove_redundant_nodes(T_next[top_bigrade][T_cur_root], get_bigrade(*it_i_minus_1),T_next);
                remove_redundant_nodes(T_next[top_bigrade][top_v], get_bigrade(*it_i_minus_1),T_next);

            }
            else
            {
                //remove_redundant_nodes(T_next[T_cur.bigrade][T_cur_root], get_bigrade(*it_prime),T_next);
                //remove_redundant_nodes(T_next[top_bigrade][T_cur_root], get_bigrade(*it_prime),T_next);
                remove_redundant_nodes(T_next[top_bigrade][top_v], get_bigrade(*it_prime),T_next);
            }
        }
        std::cout << "Status of T_next after remove_redundant_nodes:" << std::endl;
        //print_T_next(T_next, true);

        // step 8
        //remember that we have crossed this anchor
        cur_anchor->toggle();
        std::cout << "========================= cur_zs:";
        print_zeta_seq(cur_zs);
        std::cout << "========================= next_zs:";
        print_zeta_seq(next_zs);

        //if this cell does not yet have a dendrogram template, then store it now
        std::shared_ptr<Face> cur_face = (path[j])->get_face();

        if (!cur_face->has_been_visited())
        {
            std::cout << "!cur_face->has_been_visited()" << std::endl;
            cur_face->mark_as_visited();
            face_to_zs[cur_face] = next_zs;
            //Time_root& tr = T_next[get_bigrade(*it_i_plus_1)].begin()->second;
            std::cout << "T_cur.bigrade = (" << T_cur.bigrade.first << "," << T_cur.bigrade.second << ")" << std::endl;
            std::cout << "T_cur_root = " << T_cur_root << std::endl;
            std::cout << "top_v = " << top_v << std::endl;
            // check to see if top node of dendrogram is redundant
            //if (T_next[T_cur.bigrade][T_cur_root]->children.size() == 1)
            if (T_next[top_bigrade][top_v]->children.size() == 1 &&
                T_next[top_bigrade][top_v]->birth_label.size() == 0)
            {
                std::cout << "top node of dendrogram is redundant" << std::endl;
                //cur_face->get_dendrogram() = *(T_next[T_cur.bigrade][T_cur_root]->children.front());
                cur_face->get_dendrogram() = *(T_next[top_bigrade][top_v]->children.front());
                top_bigrade = (T_next[top_bigrade][top_v]->children.front())->bigrade;
                std::cout << "top_bigrade = (" << top_bigrade.first << "," << top_bigrade.second << ")" <<std::endl;
            }
            else
            {
                std::cout << "top node of dendrogram is not redundant" << std::endl;
                //cur_face->get_dendrogram() = *T_next[T_cur.bigrade][T_cur_root];
                cur_face->get_dendrogram() = *T_next[top_bigrade][top_v];
            }

            // populate num_leaves
            Time_root& den_template = cur_face->get_dendrogram();
            populate_num_leaves_and_birth_label_str(den_template);

            /*
            // Testing correctness of algorithm by comparing naive vs templates
            std::cout << " %%%%%%%%%% testing correctness of algorithm %%%%%%%%%%" << std::endl;

            Time_root naive_tr = compute_dendrogram_template(next_zs);
            Time_root template_tr = cur_face->get_dendrogram();

            std::cout << "populate_ordered_component_label(naive_tr)" << std::endl;
            populate_ordered_component_label(naive_tr);
            std::cout << "populate_ordered_component_label(template_tr)" << std::endl;
            populate_ordered_component_label(template_tr);
            std::cout << "print naive_tr: " << std::endl;
            Time_root::recursively_print_tr(naive_tr);
            std::cout << "print template tr: " << std::endl;
            Time_root::recursively_print_tr(template_tr);*/

            Time_root template_tr = cur_face->get_dendrogram();
            bool is_isomorphic = compare_template_to_naive_dendrogram(next_zs, template_tr);

            //bool is_isomorphic = is_dendrogram_isomorphic(naive_tr,template_tr);

            std::cout << "is_isomorphic = " << is_isomorphic << std::endl;
            if (!is_isomorphic)
                throw std::runtime_error("dendrograms are not isomorphic!");
        }


        //Reset T_cur
        std::cout << "about to reset T_cur" << std::endl;
        //T_cur = *T_next[T_cur.bigrade][T_cur_root]; // Store T_cur as version with possibly redundant top node
        std::cout << "top_bigrade = (" << top_bigrade.first << "," << top_bigrade.second << ")" <<std::endl;

        T_cur = cur_face->get_dendrogram();
        std::cout << "T_cur: " << std::endl;
        Time_root::recursively_print_tr(T_cur);
        std::cout << "about to reset cur_zs" << std::endl;
        cur_zs = next_zs;
        debug() << "cur_face = " << arrangement.FID(cur_face);
    }

    // restore buffer
    //std::cout.rdbuf(orig_buf);
}

bool PersistenceUpdater::compare_template_to_naive_dendrogram(zeta_seq next_zs,
                                                             Time_root& template_tr)
{
    // Testing correctness of algorithm by comparing naive vs templates
    std::cout << " %%%%%%%%%% compare_template_to_naive_dendrogram %%%%%%%%%%" << std::endl;

    Time_root naive_tr = compute_dendrogram_template(next_zs);
    //Time_root template_tr = cur_face->get_dendrogram();

    std::cout << "populate_ordered_component_label(naive_tr)" << std::endl;
    populate_ordered_component_label(naive_tr);
    std::cout << "populate_ordered_component_label(template_tr)" << std::endl;
    populate_ordered_component_label(template_tr);
    std::cout << "print naive_tr: " << std::endl;
    Time_root::recursively_print_tr(naive_tr);
    std::cout << "print template tr: " << std::endl;
    Time_root::recursively_print_tr(template_tr);

    return is_dendrogram_isomorphic(naive_tr,template_tr);
}

// populate num_leaves and birth_label_str
int PersistenceUpdater::populate_num_leaves_and_birth_label_str(Time_root& tr)
{
    // first convert birth_label to birth_label_str
    std::stringstream ss;
    for (int x : tr.birth_label)
        ss << x << ",";
    tr.birth_label_str = ss.str();

    tr.num_leaves = 0;
    if (tr.children.size() == 0)
    {
        tr.num_leaves = 1;
        return 1;
    }
    int new_num_leaves = 0;
    for (auto child : tr.children)
    {
        new_num_leaves += populate_num_leaves_and_birth_label_str(*child);
    }
    tr.num_leaves = new_num_leaves;
    return new_num_leaves;
}

void PersistenceUpdater::recursively_print_T_cur(Time_root& tr)
{
    Time_root::print_tr(tr, true);
    for (auto child : tr.children)
    {
        recursively_print_T_cur(*child);
    }
    return;
}

// removes nodes with in and out degree 1 and null birth label
void PersistenceUpdater::remove_redundant_nodes(
        std::shared_ptr<Time_root> tr_ptr,
        bigrade lower_cutoff,
        boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next)
{
    //std::cout << " =========== remove_redundant_nodes" << std::endl;
    //std::cout << "tr:" << std::endl;
    Time_root& tr = *tr_ptr;
    //Time_root::print_tr(tr,true);
    //std::cout << "lower_cutoff: (" << lower_cutoff.first << "," << lower_cutoff.second << ")" << std::endl;
    if (tr.bigrade < lower_cutoff || tr.children.size() == 0)
    {
        //std::cout << "about to return" << std::endl;
        //std::cout << "tr.bigrade = (" << tr.bigrade.first << "," << tr.bigrade.second << ")" << std::endl;
        //std::cout << "tr.children.size() == " << tr.children.size() << std::endl;
        return;
    }

    std::queue<std::shared_ptr<Time_root>> unhandled_children;
    for (int i = 0; i < tr.children.size(); i++)
        unhandled_children.push(tr.children[i]);

    while (!unhandled_children.empty())
    {
        auto child_ptr = unhandled_children.front();
        unhandled_children.pop();
        Time_root& child = *child_ptr;
        //std::cout << "child:" << std::endl;
        //Time_root::print_tr(child);
        if (child.children.size() == 1 && child.birth_label.empty())
        {
            //std::cout << "child is redundant" << std::endl;
            auto unique_child_of_child = child.children.front();
            //std::cout << "unique_child_of_child:" << std::endl;
            //Time_root::print_tr(*unique_child_of_child);
            *child_ptr = *unique_child_of_child;
            //std::cout << "updated tr's children:" << std::endl;
            //Time_root::print_children(tr.children);
            //std::cout << "T_next[tr.bigrade][tr.r].children:" << std::endl;
            //Time_root::print_children(T_next[tr.bigrade][tr.r]->children);
            unhandled_children.push(child_ptr);
            //remove_redundant_nodes(child_ptr,lower_cutoff,T_next);
        }
        else
        {
            //std::cout << "child is not redundant" << std::endl;
            remove_redundant_nodes(child_ptr, lower_cutoff, T_next);
        }
    }
}

/*
// removes nodes with in and out degree 1 and null birth label
void PersistenceUpdater::remove_redundant_nodes(
        Time_root& tr,
        bigrade lower_cutoff,
        boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next)
{
    std::cout << " =========== remove_redundant_nodes" << std::endl;
    std::cout << "tr:" << std::endl;
    Time_root::print_tr(tr,true);
    std::cout << "lower_cutoff: (" << lower_cutoff.first << "," << lower_cutoff.second << ")" << std::endl;
    if (tr.bigrade < lower_cutoff || tr.children.size() == 0)
    {
        std::cout << "about to return" << std::endl;
        std::cout << "tr.bigrade = (" << tr.bigrade.first << "," << tr.bigrade.second << ")" << std::endl;
        std::cout << "tr.children.size() == " << tr.children.size() << std::endl;
        return;
    }
    for (int i = 0; i < tr.children.size(); i++)
    {
        Time_root& child = *tr.children[i];
        std::cout << "child:" << std::endl;
        Time_root::print_tr(child);
        if (child.children.size() == 1 && child.birth_label.empty())
        {
            std::cout << "child is redundant" << std::endl;
            Time_root& unique_child_of_child = *(child.children.front());
            std::cout << "unique_child_of_child:" << std::endl;
            Time_root::print_tr(unique_child_of_child);
            tr.children[i] = std::make_shared<Time_root>(unique_child_of_child); // check that this doesn't change the original child
            std::cout << "updated tr's children:" << std::endl;
            Time_root::print_children(tr.children);
            std::cout << "T_next[tr.bigrade][tr.r].children:" << std::endl;
            Time_root::print_children(T_next[tr.bigrade][tr.r]->children);

            remove_redundant_nodes(unique_child_of_child,lower_cutoff,T_next);
        }
        else
            remove_redundant_nodes(child, lower_cutoff,T_next);
    }
}*/

// convert component labels to birth_labels
void PersistenceUpdater::convert_component_to_birth_labels(
        boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
        std::vector<bigrade> functor_bigrades)
{
    //std::cout << "convert_component_to_birth_label" << std::endl;
    //std::cout << "############# beginning status of T_next:" << std::endl;
    //print_T_next(T_next,true);
    for (bigrade big: functor_bigrades)
    {
        for (auto& r_tr : T_next[big])
        {
            Time_root& tr = *r_tr.second;
            //std::cout << "tr:" << std::endl;
            //Time_root::print_tr(tr,true);
            for (int v : tr.component_label)
            {
                tr.birth_label.insert(v);
            }
            //tr.birth_label = tr.component_label;
        }
    }
    //std::cout << "############# middle status of T_next:" << std::endl;
    //print_T_next(T_next,true);
    //std::cout << "############# begin conversion" << std::endl;
    for (bigrade big: functor_bigrades)
    {
        for (auto& r_tr : T_next[big])
        {
            Time_root& tr = *r_tr.second;
            std::cout << "tr:" << std::endl;
            Time_root::print_tr(tr,true);
            for (std::shared_ptr<Time_root> child : tr.children)
            {
                std::cout << "child:" << std::endl;
                Time_root::print_tr(*child,true);
                std::cout << "component_label:";
                for (int v : child->component_label)
                {
                    std::cout << v << ",";
                    tr.birth_label.erase(v);
                }
              /*std::stringstream ss;
                for (int x : tr.birth_label)
                    ss << x << ",";
                tr.birth_label_str = ss.str();
                std::cout << std::endl;*/
            }
        }
    }
    //std::cout << "############# end status of T_next:" << std::endl;
    //print_T_next(T_next,true);
}

// copy upper part of T_cur
void PersistenceUpdater::copy_upper_dendrogram(Time_root& tr,
                                               boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                                               bigrade big_i_plus_1,
                                               std::vector<bigrade_root>& parents_of_iplus1,
                                               std::map<int, Subset_map>& oracle,
                                               std::shared_ptr<Time_root> parent_ptr)
{
    //std::cout << "copy_upper_dendrogram" << std::endl;
    //std::cout << "tr:" << std::endl;
    //Time_root::print_tr(tr, true);
    //std::cout << "children: " << std::endl;
    //Time_root::print_children(tr.children, true);
    //std::cout << "big_i_plus_1: (" << big_i_plus_1.first << "," << big_i_plus_1.second << ")" << std::endl;

    std::shared_ptr<Time_root> tr_copy_ptr = std::make_shared<Time_root>(tr);
    Time_root& tr_copy = *tr_copy_ptr;

    if (tr.bigrade < big_i_plus_1)
    {
        std::cout << "about to return" << std::endl;
        if (parent_ptr != NULL)
        {
            int h_r = Find(oracle, tr.r);
            //std::cout << "must update parent_ptr's children before returning" << std::endl;
            //std::cout << "*T_next[(" << big_i_plus_1.first << "," << big_i_plus_1.second << ")][" << h_r << "] = " << std::endl;
            Time_root::print_tr(*T_next[big_i_plus_1][h_r]);
            (*parent_ptr).children.push_back(T_next[big_i_plus_1][h_r]);
        }
        return; // is this return condition even the right one?
    }
    //Time_root tr_copy = tr;

    //std::cout << "tr_copy:" << std::endl;
    //Time_root::print_tr(tr_copy);
    tr_copy.children.clear();
    //std::cout << "tr_copy.children: " << std::endl;
    //Time_root::print_children(tr_copy.children);
    std::unordered_set<int> type_U_seen; // need because children of type U may have merged before i+1
    for (std::shared_ptr<Time_root> child_ptr : tr.children)
    {
        Time_root& child = *child_ptr;
        //std::cout << "child: " << std::endl;
        //Time_root::print_tr(child);
        bool child_is_at_cutoff = (child.bigrade == big_i_plus_1 ||
                                  (tr.bigrade > big_i_plus_1 && child.bigrade < big_i_plus_1));
        // child_is_at_cutoff iff child \in V(U) \cup V(i+1)
        if (child_is_at_cutoff)
        {
            //Time_root& child_twin = *T_next[child.bigrade][child.r]; //child_twin represents the same node in the functor dendrogram of T^f, image of child under XY in paper
            parents_of_iplus1.push_back(bigrade_root(tr.bigrade,tr.r));
            if (child.bigrade == big_i_plus_1)
            {
                //std::cout << "child is type i + 1" << std::endl;
                int h_r = Find(oracle, child.r);
                //std::cout << "h_r = " << h_r << std::endl;
                /*std::cout << "T_next[child.bigrade][child.r] = "
                          << "T_next[(" << child.bigrade.first << ","
                          << child.bigrade.second << ")]["
                          << h_r << "] = " << std::endl;
                Time_root::print_tr(*T_next[child.bigrade][h_r]);
                std::cout << "able to print_tr" << std::endl;*/
                tr_copy.children.push_back(T_next[child.bigrade][h_r]);
            }
            else
            {
                //std::cout << "child is type U" << std::endl;
                int h_r = Find(oracle, child.r);
                //std::cout << "h_r = " << h_r << std::endl;
                if (type_U_seen.find(h_r) == type_U_seen.end())
                {
                    tr_copy.children.push_back(T_next[big_i_plus_1][h_r]);
                    type_U_seen.insert(h_r);
                }
            }

        }
        else
        {
            //Time_root child_copy = child;
            std::shared_ptr<Time_root> child_copy_ptr = std::make_shared<Time_root>(child);
            //std::shared_ptr<Time_root> child_copy_ptr = std::make_shared<Time_root> (child_copy);
            //tr_copy.children.push_back(child_copy_ptr);
            copy_upper_dendrogram(*child_copy_ptr, T_next, big_i_plus_1, parents_of_iplus1, oracle, tr_copy_ptr);
        }
    }
    if (parent_ptr != NULL)
        (*parent_ptr).children.push_back(tr_copy_ptr);
    //std::cout << "about to assign T_next[(" << tr.bigrade.first << "," << tr.bigrade.second << ")][" << tr.r << "] = tr_copy_ptr" << std::endl;
    //T_next[tr.bigrade][tr.r] = std::make_shared<Time_root>(tr_copy);
    T_next[tr.bigrade][tr.r] = tr_copy_ptr;
    /*std::cout << "*tr_copy_ptr = " << std::endl;
    Time_root::print_tr(*tr_copy_ptr);
    std::cout << "*tr_copy_ptr.children = " << std::endl;
    Time_root::print_children((*tr_copy_ptr).children);
    std::cout << "sucessfully printed_children" << std::endl;*/

}

// connect different levels (big_1 and big_2) of T_next
void PersistenceUpdater::connect_levels(std::map<int, Subset_map>& oracle,
                                        boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                                        bigrade big_1,
                                        bigrade big_2)
{
    for (auto& v_tr : T_next[big_1])
    {
        int v = v_tr.first;
        Time_root& tr_1 = *v_tr.second;
        int h_v = Find(oracle,v);
        Time_root& tr_2 = *T_next[big_2][h_v];
        tr_2.children.push_back(v_tr.second);
    }
}

// instantiate a level (at the given bigrade) of T_next
void PersistenceUpdater::instantiate_level(std::map<int, Subset_map>& oracle,
                                           boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                                           std::pair<unsigned,unsigned> big)
{
    // can optimize this part by remembering another table "oracle_heads" to remember the heads
    std::unordered_set<int> heads_with_tr;
    for (auto kv : oracle)
    {
        int v = kv.first;
        int h_v = Find(oracle, v); //head of v (r_v in paper)

        if (heads_with_tr.find(h_v) == heads_with_tr.end())
        {
            heads_with_tr.insert(h_v);
            Time_root tr_h_v;
            tr_h_v.bigrade = big;
            tr_h_v.r = h_v;
            //T_next[bigrade_root(bigrade,h_v)] = tr_h_v;
            T_next[big][h_v] = std::make_shared<Time_root>(tr_h_v);
        }
        //T_next[bigrade_root(bigrade,h_v)].component_label.insert(v);
        T_next[big][h_v]->component_label.insert(v);
    }
}

// populate component labels of a level (at the given bigrade) in T_next
void PersistenceUpdater::populate_component_labels(std::map<int, Subset_map>& oracle,
                                                   boost::unordered::unordered_map<bigrade,boost::unordered::unordered_map<int, std::shared_ptr<Time_root>>>& T_next,
                                                   bigrade big)
{
    for (auto kv : oracle)
    {
        int v = kv.first;
        std::cout << "big = (" << big.first << "," << big.second << ")" << std::endl;
        std::cout << "v = " << v << std::endl;
        int h_v = Find(oracle, v); //head of v (r_v in paper)
        std::cout << "h_v = " << h_v << std::endl;
        T_next[big][h_v]->component_label.insert(v);
    }
}

//get bigrade from TemplatePointsMatrixEntry
std::pair<unsigned,unsigned> PersistenceUpdater::get_bigrade(std::shared_ptr<TemplatePointsMatrixEntry> u)
{
    return std::pair<unsigned,unsigned>(u->x,u->y);
}

//set oracle
void PersistenceUpdater::set_oracle(boost::unordered::unordered_map<std::pair<int,int>,std::map<int, Subset_map>>& ora)
{
    oracle = ora;
}

void PersistenceUpdater::set_cs(computing_s* comp_s)
{
    cs = comp_s;
}

//classify nodes of dendrogram during updates
void PersistenceUpdater::classify_nodes(std::vector<Time_root>& L,
                                        std::vector<Time_root>& i_minus_1,
                                        std::vector<Time_root>& i_plus_1,
                                        std::vector<Time_root>& U,
                                        Time_root& tr,
                                        std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_prime,
                                        bool it_i_minus_1_exists,
                                        bool it_i_plus_1_exists)
{
    //std::cout << ">>>>>>>>>>>>>>> classify_nodes" << std::endl;
    std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_i_minus_1;
    if (it_i_minus_1_exists)
    {
        it_i_minus_1 = std::prev(it_prime,1);
        //std::cout << "it_i_minus_1 : ";
        //print_iterator(it_i_minus_1);
    }
    std::vector<std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it_i_plus_1;
    if (it_i_plus_1_exists)
    {
        it_i_plus_1 = std::next(it_prime,1);
        //std::cout << "it_i_plus_1 : ";
        //print_iterator(it_i_plus_1);
    }

    //std::cout << "classify_nodes 2" << std::endl;
    //std::cout << "tr : " << std::endl;
    //Time_root::print_tr(tr);
    //Time_root::print_children(tr.children);
    if (tr.children.size() == 0)
        return;
    for (std::shared_ptr<Time_root> child_ptr : tr.children)
    {
        Time_root& child = *child_ptr;
        if (it_i_plus_1_exists)
        {
            if (child.bigrade == get_bigrade(*it_i_plus_1))
                i_plus_1.push_back(child);
            if (tr.bigrade > get_bigrade(*it_i_plus_1) && child.bigrade < get_bigrade(*it_i_plus_1))
                U.push_back(child);
        }
        if (it_i_minus_1_exists)
        {
            //std::cout << "reached classification of i_minus_1 and L" << std::endl;
            //std::cout << "child : " << std::endl;
            //Time_root::print_tr(child);
            if (child.bigrade == get_bigrade(*it_i_minus_1))
                i_minus_1.push_back(child);
            if (tr.bigrade > get_bigrade(*it_i_minus_1) && child.bigrade < get_bigrade(*it_i_minus_1))
                L.push_back(child);
        }
        classify_nodes(L,i_minus_1,i_plus_1,U,child,it_prime,it_i_minus_1_exists,it_i_plus_1_exists);
    }
}

//stores a dendrogram template in a 2-cell of the arrangement
void PersistenceUpdater::store_initial_dendrogram_template(std::shared_ptr<Face> cell, zeta_seq& zs)
{
    //mark this cell as visited

    cell->mark_as_visited();

    Time_root& tr = cell->get_dendrogram();

    // 1-critical version SimplexSet vertices = bifiltration.get_ordered_simplices();
    // 1-critical version SimplexSet edges = bifiltration.get_ordered_high_simplices();
    boost::unordered::unordered_map<unsigned, double> vertex_appearance;
    boost::unordered::unordered_map<unsigned, std::pair<double,double>> vertex_bigrade;
    boost::unordered::unordered_map<double, std::pair<double,double>> appearance_to_bigrade; // map from time of appearance to point in zeta seq

    SimplexInfo* vertices_to_grades = bifiltration_data.getSimplices(0);
    SimplexInfo* edges_to_grades = bifiltration_data.getSimplices(1);

    for (auto vertexVec_grade : *vertices_to_grades)
    {
        unsigned v = vertexVec_grade.first[0];
        AppearanceGrades grades = vertexVec_grade.second;
        for (Grade g : grades)
        {
            double birth = g.y;
            std::pair<double,double> big (g.x,g.y);
            debug() << "(v,x,y,birth) = (" <<
                    v << " , " <<
                    g.x << " , " <<
                    g.y << " , " <<
                    birth << ")";
            if (vertex_appearance.find(v) == vertex_appearance.end())
            {
                vertex_appearance[v] = birth;
                vertex_bigrade[v] = big;
            }
            else
            {
                if (birth < vertex_appearance[v])
                {
                    vertex_appearance[v] = birth;
                    vertex_bigrade[v] = big;
                }
            }
        }
    }

    for (auto entry : zs)
    {
        appearance_to_bigrade[entry->y] = std::pair<double,double>(entry->x,entry->y);
    }

    /* 1-critical version
    for (auto st_node : vertices)
    {
        unsigned v = st_node->get_vertex();
        double x = st_node->grade_x();
        double y = st_node->grade_y();
        double birth = y;
        std::pair<double,double> bigrade (x,y);

        /*
        debug() << "(v,x,y,birth) = (" <<
                v << " , " <<
                x << " , " <<
                y << " , " <<
                birth << ")";
        if (vertex_appearance.find(v) == vertex_appearance.end())
        {
            vertex_appearance[v] = birth;
            vertex_bigrade[v] = bigrade;
        }
        else
        {
            if (birth < vertex_appearance[v])
            {
                vertex_appearance[v] = birth;
                vertex_bigrade[v] = bigrade;
            }
        }
    }*/

    boost::unordered::unordered_map<std::pair<unsigned,unsigned>, double> edge_appearance;
    boost::unordered::unordered_map<std::pair<unsigned,unsigned>, std::pair<double,double>> edge_bigrade;

    for (auto edgeVec_grade : *edges_to_grades)
    {
        std::vector<int> edge_vector = edgeVec_grade.first;
        AppearanceGrades grades = edgeVec_grade.second;
        std::pair<unsigned,unsigned> edge (edge_vector[0], edge_vector[1]);
        for (Grade g : grades)
        {
            double birth = g.y;
            std::pair<double,double> big (g.x,g.y);
            if (edge_appearance.find(edge) == edge_appearance.end())
            {
                edge_appearance[edge] = birth;
                edge_bigrade[edge] = big;
            }
            else
            {
                if (birth < edge_appearance[edge])
                {
                    edge_appearance[edge] = birth;
                    edge_bigrade[edge] = big;
                }
            }
        }
    }

    /* 1-critical version
    for (auto st_node : edges)
    {
        int global_index = st_node->global_index();
        std::vector<int> edge_vector = bifiltration.find_vertices(global_index);
        std::pair<unsigned,unsigned> edge (edge_vector[0], edge_vector[1]);
        double x = st_node->grade_x();
        double y = st_node->grade_y();
        double birth = y;

        // case where edge does not live in T_1 \subset S
        if (appearance_to_bigrade.find(birth) == appearance_to_bigrade.end())
        {
            continue;
        }
        std::pair<double,double> bigrade (x,y);

        if (edge_appearance.find(edge) == edge_appearance.end())
        {
            edge_appearance[edge] = birth;
            edge_bigrade[edge] = bigrade;
        }
        else
        {
            if (birth < edge_appearance[edge])
            {
                edge_appearance[edge] = birth;
                edge_bigrade[edge] = bigrade;
            }
        }
    }*/

    DenseGRAPH<Edge>* graph = new DenseGRAPH<Edge>(vertex_appearance.size(), true);
    debug() << "print edges of graph";
    for (auto edge_wt : edge_appearance)
    {
        std::pair<unsigned,unsigned> e = edge_wt.first;
        unsigned v = e.first;
        unsigned w = e.second;
        graph-> insert(EdgePtr( new Edge(v,w,edge_wt.second,edge_bigrade[e])));
        /*debug() << "(v,w,edge_wt,edge_bigrade) = (" <<
                   v << " , " <<
                   w << " , " <<
                   edge_wt.second << " , (" <<
                   edge_bigrade[e].first << " , " << edge_bigrade[e].second << "))";*/

    }

    std::vector<EdgePtr> mst;
    boost::unordered::unordered_map<std::pair<double,int>, Time_root> time_root_to_tr;

    debug() << "START store_initial_dendrogram_template";
    std::vector<Time_root> upper_trs;
    Dendrogram_data::compute_dendrogram(graph, mst, time_root_to_tr, tr, vertex_appearance, vertex_bigrade, appearance_to_bigrade, upper_trs);
    debug() << "upper_trs: ";
    //Time_root::print_trs(upper_trs, true);
    for ( std::vector<EdgePtr>::const_iterator it = mst.begin(); it != mst.end(); ++it )
    {
        /*
        debug() << it->get()->v << " "
                << it->get()->w << " "
                << it->get()->wt << "\n";*/
    }
    if (upper_trs.size() > 1)
    {
        Time_root forest_connector;
        forest_connector.r = is_forest_connector;
        forest_connector.bigrade = bigrade(UINT_MAX,UINT_MAX);
        for (Time_root& upper_tr : upper_trs)
        {
            forest_connector.children.push_back(std::make_shared<Time_root>(upper_tr));
            forest_connector.num_leaves += upper_tr.num_leaves;
        }
        tr = forest_connector;
        std::cout << "forest_connector : " << std::endl;
        Time_root::print_tr(forest_connector);
        Time_root::print_children(forest_connector.children);
    }
    else // upper_trs.size() == 1
        tr = upper_trs.front();

    debug() << "END store_initial_dendrogram_template";
    for (auto t_big : appearance_to_bigrade)
    {
        auto t = t_big.first;
        auto big = t_big.second;
        /*
        std::cout << "t = " << t << " , " << "big = ("
                  << big.first << "," << big.second << ")"
                  << std::endl;*/
    }
    return;

}



//computes and stores a barcode template in each 2-cell of arrangement
//resets the matrices and does a standard persistence calculation for expensive crossings
void PersistenceUpdater::store_barcodes_with_reset(std::vector<std::shared_ptr<Halfedge>>& path, Progress& progress)
{

    // PART 1: GET THE BOUNDARY MATRICES WITH PROPER SIMPLEX ORDERING

    Timer timer;

    //initialize the lift map from simplex grades to LUB-indexes
    if (verbosity >= 10) {
        debug() << "  Mapping low simplices:";
    }
    std::cout << "break 1" << std::endl;
    IndexMatrix* ind_low = bifiltration.get_index_mx(dim); //can we improve this with something more efficient than IndexMatrix?
    std::cout << "break 2" << std::endl;
    store_multigrades(ind_low, true);
    std::cout << "break 3" << std::endl;

    if (verbosity >= 10) {
        debug() << "  Mapping high simplices:";
    }
    IndexMatrix* ind_high = bifiltration.get_index_mx(dim + 1); //again, could be improved?
    std::cout << "break 4" << std::endl;
    store_multigrades(ind_high, false);
    std::cout << "break 5" << std::endl;

    //get the proper simplex ordering
    std::vector<int> low_simplex_order; //this will be a map : dim_index --> order_index for dim-simplices; -1 indicates simplices not in the order
    unsigned num_low_simplices = build_simplex_order(ind_low, true, low_simplex_order);
    delete ind_low;
    std::cout << "break 5.1" << std::endl;

    std::vector<int> high_simplex_order; //this will be a map : dim_index --> order_index for (dim+1)-simplices; -1 indicates simplices not in the order
    unsigned num_high_simplices = build_simplex_order(ind_high, false, high_simplex_order);
    delete ind_high;
    std::cout << "break 5.2" << std::endl;

    //get boundary matrices (R) and identity matrices (U) for RU-decomposition
    R_low = bifiltration.get_boundary_mx(low_simplex_order, num_low_simplices);
    std::cout << "break 5.25" << std::endl;
    R_high = bifiltration.get_boundary_mx(low_simplex_order, num_low_simplices, high_simplex_order, num_high_simplices);
    std::cout << "break 5.3" << std::endl;
    //print runtime data
    if (verbosity >= 4) {
        debug() << "  --> computing initial order on simplices and building the boundary matrices took"
                << timer.elapsed() << "milliseconds";
    }

    //copy the boundary matrices (R) for fast reset later
    timer.restart();
    MapMatrix_Perm* R_low_initial = new MapMatrix_Perm(*R_low);
    MapMatrix_Perm* R_high_initial = new MapMatrix_Perm(*R_high);
    std::cout << "break 5.4" << std::endl;
    if (verbosity >= 4) {
        debug() << "  --> copying the boundary matrices took"
                << timer.elapsed() << "milliseconds";
    }

    //initialize the permutation vectors
    perm_low.resize(R_low->width());
    inv_perm_low.resize(R_low->width());
    perm_high.resize(R_high->width());
    inv_perm_high.resize(R_high->width());
    std::cout << "break 5.5" << std::endl;
    for (unsigned j = 0; j < perm_low.size(); j++) {
        perm_low[j] = j;
        inv_perm_low[j] = j;
    }
    for (unsigned j = 0; j < perm_high.size(); j++) {
        perm_high[j] = j;
        inv_perm_high[j] = j;
    }
    std::cout << "break 6" << std::endl;

    // PART 2: INITIAL PERSISTENCE COMPUTATION (RU-decomposition)

    timer.restart();

    //initial RU-decomposition
    U_low = R_low->decompose_RU();
    U_high = R_high->decompose_RU();

    unsigned total_time_for_resets = timer.elapsed();
    if (verbosity >= 4) {
        debug() << "  --> computing the RU decomposition took" << total_time_for_resets << "milliseconds";
    }

    //store the barcode template in the first cell
    std::shared_ptr<Face> first_cell = arrangement.topleft->get_twin()->get_face();
    store_barcode_template(first_cell);

    if (verbosity >= 4) {
        debug() << "Initial persistence computation in cell " << arrangement.FID(first_cell);
    }
    std::cout << "break 7" << std::endl;
    // PART 3: TRAVERSE THE PATH AND UPDATE PERSISTENCE AT EACH STEP

    if (verbosity >= 2) {
        debug() << "TRAVERSING THE PATH USING THE RESET ALGORITHM: path has" << path.size() << "steps";
    }

    //data members for analyzing the computation
    unsigned long total_transpositions = 0;
    unsigned total_time_for_transpositions = 0; //NEW
    unsigned number_of_resets = 1; //we count the initial RU-decomposition as the first reset
    int max_time = 0;

    // choose the initial value of the threshold intelligently
    unsigned long threshold;
    choose_initial_threshold(total_time_for_resets, total_transpositions, total_time_for_transpositions, threshold); 
        //if the number of swaps might exceed this threshold, then we will do a persistence calculation from scratch instead of vineyard updates
    if (verbosity >= 4) {
        debug() << "initial reset threshold set to" << threshold;
    }

    timer.restart();


    //traverse the path
    Timer steptimer;
    for (unsigned i = 0; i < path.size(); i++) {
        progress.progress(i); //update progress bar

        steptimer.restart(); //time update at each step of the path
        unsigned long num_trans = 0; //count of how many transpositions we will have to do if we do vineyard updates
        unsigned long swap_counter = 0; //count of how many transpositions we actually do

        //determine which anchor is represented by this edge
        std::shared_ptr<Anchor> cur_anchor = (path[i])->get_anchor();
        std::shared_ptr<TemplatePointsMatrixEntry> at_anchor = cur_anchor->get_entry();

        //get equivalence classes for this anchor
        std::shared_ptr<TemplatePointsMatrixEntry> down = at_anchor->down;
        std::shared_ptr<TemplatePointsMatrixEntry> left = at_anchor->left;

        //if this is a strict anchor, then swap simplices
        if (down != nullptr && left != nullptr) //then this is a strict anchor and some simplices swap
        {
            if (verbosity >= 6) {
                debug() << "  step " << i << " of path: crossing (strict) anchor at (" << cur_anchor->get_x() << ", " << cur_anchor->get_y() << ") into cell " << arrangement.FID((path[i])->get_face()) << "; edge weight: " << cur_anchor->get_weight();
            }

            //find out how many transpositions we will have to process if we do vineyard updates
            num_trans = count_transpositions(at_anchor, cur_anchor->is_above());

            if (cur_anchor->is_above()) //then the anchor is crossed from below to above
            {
                remove_lift_entries(at_anchor); //this block of the partition might become empty
                remove_lift_entries(down); //this block of the partition will move

                if (num_trans < threshold) //then do vineyard updates
                {
                    swap_counter += split_grade_lists(at_anchor, left, true); //move grades that come before left from anchor to left -- vineyard updates
                    swap_counter += move_columns(down, left, true); //swaps blocks of columns at down and at left -- vineyard updates
                } else //then reset the matrices
                {
                    split_grade_lists_no_vineyards(at_anchor, left, true); //only updates the xiSupportMatrix and permutation vectors; no vineyard updates
                    update_order_and_reset_matrices(down, left, true, R_low_initial, R_high_initial); //recompute the RU-decomposition
                }

                merge_grade_lists(at_anchor, down); //move all grades from down to anchor
                add_lift_entries(at_anchor); //this block of the partition might have previously been empty
                add_lift_entries(left); //this block of the partition moved
            }
            else //then anchor is crossed from above to below
            {
                remove_lift_entries(at_anchor); //this block of the partition might become empty
                remove_lift_entries(left); //this block of the partition will move

                if (num_trans < threshold) //then do vineyard updates
                {
                    swap_counter += split_grade_lists(at_anchor, down, false); //move grades that come before left from anchor to left -- vineyard updates
                    swap_counter += move_columns(left, down, false); //swaps blocks of columns at down and at left -- vineyard updates
                } else //then reset the matrices
                {
                    split_grade_lists_no_vineyards(at_anchor, down, false); //only updates the xiSupportMatrix and permutation vectors; no vineyard updates
                    update_order_and_reset_matrices(left, down, false, R_low_initial, R_high_initial); //recompute the RU-decomposition
                }

                merge_grade_lists(at_anchor, left); //move all grades from down to anchor
                add_lift_entries(at_anchor); //this block of the partition might have previously been empty
                add_lift_entries(down); //this block of the partition moved
            }
        }
        else //this is a non-strict anchor, and we just have to split or merge equivalence classes
        {
            if (verbosity >= 6) {
                debug() << "  step " << i << " of path: crossing (non-strict) anchor at (" << cur_anchor->get_x() << ", " << cur_anchor->get_y() << ") into cell " << arrangement.FID((path[i])->get_face()) << "; edge weight: " << cur_anchor->get_weight();
            }

            std::shared_ptr<TemplatePointsMatrixEntry> generator = at_anchor->down;
            if (generator == nullptr)
                generator = at_anchor->left;

            if ((cur_anchor->is_above() && generator == at_anchor->down) || (!cur_anchor->is_above() && generator == at_anchor->left))
            //then merge classes -- there will never be any transpositions in this case
            {
                remove_lift_entries(generator);
                merge_grade_lists(at_anchor, generator);
                add_lift_entries(at_anchor); //this is necessary in case the class was previously empty
            } else //then split classes
            {
                //find out how many transpositions we will have to process if we do vineyard updates
                unsigned junk = 0;
                bool horiz = (generator == at_anchor->left);
                count_transpositions_from_separations(at_anchor, generator, horiz, true, num_trans, junk);
                count_transpositions_from_separations(at_anchor, generator, horiz, false, num_trans, junk);

                //now do the updates
                remove_lift_entries(at_anchor); //this is necessary because the class corresponding to at_anchor might become empty

                if (num_trans < threshold) //then do vineyard updates
                    swap_counter += split_grade_lists(at_anchor, generator, horiz);
                else //then reset the matrices
                {
                    split_grade_lists_no_vineyards(at_anchor, generator, horiz); //only updates the xiSupportMatrix; no vineyard updates
                    update_order_and_reset_matrices(R_low_initial, R_high_initial); //recompute the RU-decomposition
                }

                add_lift_entries(at_anchor);
                add_lift_entries(generator);
            }
        }

        //remember that we have crossed this anchor
        cur_anchor->toggle();

        //if this cell does not yet have a barcode template, then store it now
        std::shared_ptr<Face> cur_face = (path[i])->get_face();
        if (!cur_face->has_been_visited())
            store_barcode_template(cur_face);

        //print/store data for analysis
        int step_time = steptimer.elapsed();

        if (num_trans < threshold) //then we did vineyard-updates
        {
            if (verbosity >= 6) {
                debug() << "  --> this step took" << step_time << "milliseconds and involved" << swap_counter << "transpositions; estimate was" << num_trans;
            }
            //TESTING: if (swap_counter != num_trans)
            //    debug() << "    ========>>> ERROR: transposition count doesn't match estimate!";

            if (swap_counter > 0) //don't track time for overhead that doesn't result in any transpositions
            {
                total_transpositions += swap_counter;
                total_time_for_transpositions += step_time;
            }
        } else {
            if (verbosity >= 6) {
                debug() << "  --> this step took" << step_time << "milliseconds; reset matrices to avoid" << num_trans << "transpositions";
            }
            //TESTING: if (swap_counter > 0)
            //    debug() << "    ========>>> ERROR: swaps occurred on a matrix reset!";
            number_of_resets++;
            total_time_for_resets += step_time;
        }

        if (step_time > max_time)
            max_time = step_time;

        //update the treshold
        if(swap_counter > 0 || num_trans >= threshold) {
            threshold = (unsigned long)(((double)total_transpositions / total_time_for_transpositions) * ((double)total_time_for_resets / number_of_resets));
            if (verbosity >= 6) {
                // debug() << "===>>> UPDATING THRESHOLD:";
                // debug() << "    total_trans: " << total_transpositions;
                // debug() << "    total_time_for_trans: " << total_time_for_transpositions;
                // debug() << "    total time for resets: " << total_time_for_resets;
                // debug() << "    number of resets:" << number_of_resets;
                debug() << "  -- new threshold:" << threshold;
            }
        }
    } //end path traversal

    //print runtime data
    if (verbosity >= 2) {
        debug() << "BARCODE TEMPLATE COMPUTATION COMPLETE: path traversal and persistence updates took" << timer.elapsed() << "milliseconds";
        if (verbosity >= 4) {
            debug() << "    max time per anchor crossing:" << max_time;
            debug() << "    total number of transpositions:" << total_transpositions;
            debug() << "    matrices were reset" << number_of_resets << "times when estimated number of transpositions exceeded" << threshold;
            if (number_of_resets > 0) {
                debug() << "    average time for reset:" << (total_time_for_resets / number_of_resets) << "milliseconds";
            }
        }
    }

    // PART 4: CLEAN UP

    //delete R_low;
    //delete R_high;
    delete U_low;
    delete U_high;

} //end store_barcodes_with_reset()

//function to set all edge weights to uniform value for dendrogram testing
void PersistenceUpdater::set_uniform_anchor_weights(std::vector<std::shared_ptr<Halfedge>>& path)
{
    for (unsigned i = 0; i < path.size(); i++)
    {
        //determine which anchor is represented by this edge
        std::shared_ptr<Anchor> cur_anchor = (path[i])->get_anchor();
        //std::shared_ptr<TemplatePointsMatrixEntry> at_anchor = cur_anchor->get_entry();

        if (verbosity >= 8) {
            debug() << "  step" << i << "of the short path: crossing anchor at (" << cur_anchor->get_x() << "," << cur_anchor->get_y() << ") into cell" << arrangement.FID((path[i])->get_face());
        }

        //store data
        cur_anchor->set_weight(1); //we expect that each separation produces a transposition about 25% of the time

        if (verbosity >= 8) {
            debug() << "     edge weight:" << cur_anchor->get_weight();
        }
    }
}


//function to set the "edge weights" for each anchor line
void PersistenceUpdater::set_anchor_weights(std::vector<std::shared_ptr<Halfedge>>& path)
{
    // PART 1: GET THE PROPER SIMPLEX ORDERING

    //initialize the lift map from simplex grades to LUB-indexes
    if (verbosity >= 10) {
        debug() << "  Mapping low simplices:";
    }
    IndexMatrix* ind_low = bifiltration.get_index_mx(dim); //can we improve this with something more efficient than IndexMatrix?
    store_multigrades(ind_low, true);
    delete ind_low;

    if (verbosity >= 10) {
        debug() << "  Mapping high simplices:";
    }
    IndexMatrix* ind_high = bifiltration.get_index_mx(dim + 1); //again, could be improved?
    store_multigrades(ind_high, false);
    delete ind_high;

    // PART 2: TRAVERSE THE PATH AND COUNT SWITCHES & SEPARATIONS AT EACH STEP

    for (unsigned i = 0; i < path.size(); i++) {
        unsigned long switches = 0;
        unsigned long separations = 0;

        //determine which anchor is represented by this edge
        std::shared_ptr<Anchor> cur_anchor = (path[i])->get_anchor();
        std::shared_ptr<TemplatePointsMatrixEntry> at_anchor = cur_anchor->get_entry();

        if (verbosity >= 8) {
            debug() << "  step" << i << "of the short path: crossing anchor at (" << cur_anchor->get_x() << "," << cur_anchor->get_y() << ") into cell" << arrangement.FID((path[i])->get_face());
        }

        //if this is a strict anchor, then there can be switches and separations
        if (at_anchor->down != nullptr && at_anchor->left != nullptr) //then this is a strict anchor
        {
            count_switches_and_separations(at_anchor, cur_anchor->is_above(), switches, separations);
        } else //this is a non-strict anchor, so there can be separations but not switches
        {
            std::shared_ptr<TemplatePointsMatrixEntry> generator = (at_anchor->down != nullptr) ? at_anchor->down : at_anchor->left;

            if ((cur_anchor->is_above() && generator == at_anchor->down) || (!cur_anchor->is_above() && generator == at_anchor->left))
            //then merge classes
            {
                separations += generator->low_count * at_anchor->low_count + generator->high_count * at_anchor->high_count;
                merge_grade_lists(at_anchor, generator);
            } else //then split classes
            {
                do_separations(at_anchor, generator, (at_anchor->y == generator->y));
                separations += generator->low_count * at_anchor->low_count + generator->high_count * at_anchor->high_count;
            }
        }

        //store data
        cur_anchor->set_weight(switches + separations / 4); //we expect that each separation produces a transposition about 25% of the time

        if (verbosity >= 8) {
            debug() << "     edge weight:" << cur_anchor->get_weight() << "; (" << switches << "," << separations << ")";
        }

    } //end path traversal
} //end set_anchor_weights()

//function to clear the levelset lists -- e.g., following the edge-weight calculation
void PersistenceUpdater::clear_levelsets()
{
    template_points_matrix.clear_grade_lists();
}

//stores multigrade info for the persistence computations (data structures prepared with respect to a near-vertical line positioned to the right of all \xi support points)
//  that is, this function creates the level sets of the lift map
//  low is true for simplices of dimension hom_dim, false for simplices of dimension hom_dim+1
//NOTE: this function has been updated for the new (unfactored) lift map of August 2015
void PersistenceUpdater::store_multigrades(IndexMatrix* ind, bool low)
{
    if (verbosity >= 8) {
        debug() << "STORING MULTIGRADES: low =" << low;
    }

    //initialize linked list to track the "frontier"
    typedef std::list<std::shared_ptr<TemplatePointsMatrixEntry>> Frontier;
    Frontier frontier;

    //loop through rows of TemplatePointsMatrix, from top to bottom
    for (unsigned y = ind->height(); y-- > 0;) //y counts down from (ind->height() - 1) to 0
    {
        //update the frontier for row y:
        //  if the last element of frontier has the same x-coord as cur, then replace that element with cur
        //  otherwise, append cur to the end of frontier
        std::shared_ptr<TemplatePointsMatrixEntry> cur = template_points_matrix.get_row(y);
        if (cur != nullptr) {
            Frontier::iterator it = frontier.end(); //the past-the-end element of frontier
            if (it != frontier.begin()) //then the frontier is not empty
            {
                --it; //the last element of frontier
                if ((*it)->x == cur->x) //then erase the last element of frontier
                    frontier.erase(it);
            }

            //append cur to the end of the frontier
            frontier.push_back(cur);
        }

        //store all multigrades and simplices whose y-grade is y
        Frontier::iterator it = frontier.begin();
        for (unsigned x = ind->width(); x-- > 0;) //x counts down from (ind->width() - 1) to 0
        {
            //get range of column indexes for simplices at multigrade (x,y)
            int last_col = ind->get(y, x); //arguments are row, then column
            int first_col = -1;
            if (x > 0)
                first_col = ind->get(y, x - 1);
            else if (y > 0)
                first_col = ind->get(y - 1, ind->width() - 1);

            //if there are any simplices at (x,y),
            //    and if x is not greater than the x-coordinate of the rightmost element of the frontier,
            //    then map multigrade (x,y) to the last element of the frontier such that x <= (*it)->x
            if (last_col > first_col && it != frontier.end() && x <= (*it)->x) {
                //advance the iterator to the first element of the frontier such that (*it)->x < x
                while (it != frontier.end() && (*it)->x >= x)
                    ++it;

                //back up one position, to the last element of the frontier such that (*it)->x >= x
                --it;

                //now map the multigrade to the xi support entry
                (*it)->add_multigrade(x, y, last_col - first_col, last_col, low);

                if (verbosity >= 10) {
                    debug() << "    simplices at (" << x << "," << y << "), in columns" << (first_col + 1) << "to" << last_col << ", mapped to TemplatePointsMatrixEntry at (" << (*it)->x << ", " << (*it)->y << ")";
                }
            }
        } //end x loop
    } //end y loop
} //end store_multigrades()

//finds the proper order of simplexes for the persistence calculation (with respect to a near-vertical line positioned to the right of all \xi support points)
//  NOTE: within each equivalence class, multigrades will occur in lexicographical order
//  PARAMETERS:
//    low is true for simplices of dimension hom_dim, false for simplices of dimension hom_dim+1
//    simplex_order will be filled with a map : dim_index --> order_index for simplices of the given dimension
//           If a simplex with dim_index i does not appear in the order (i.e. its grade is not less than the LUB of all xi support points), then simplex_order[i] = -1.
//  RETURN VALUE: the number of simplices in the order
unsigned PersistenceUpdater::build_simplex_order(IndexMatrix* ind, bool low, std::vector<int>& simplex_order)
{
    //count the number of simplices that will be in the order (i.e. simplices with grades less than the LUB of all xi support points)
    unsigned num_simplices = 0;
    for (unsigned row = 0; row < template_points_matrix.height(); row++) {
        std::shared_ptr<TemplatePointsMatrixEntry> cur = template_points_matrix.get_row(row);
        if (cur == nullptr)
            continue;
        std::list<std::shared_ptr<Multigrade>>* mgrades = (low) ? &(cur->low_simplices) : &(cur->high_simplices);
        for (std::list<std::shared_ptr<Multigrade>>::iterator it = mgrades->begin(); it != mgrades->end(); ++it) {
            num_simplices += (*it)->num_cols;
        }
    }

    //we will create the map starting by identifying the order index of each simplex, starting with the last simplex
    unsigned o_index = num_simplices - 1;

    //prepare the vector
    simplex_order.clear();
    simplex_order.resize(ind->last() + 1, -1); //all entries -1 by default

    //consider the rightmost TemplatePointsMatrixEntry in each row
    for (unsigned row = template_points_matrix.height(); row-- > 0;) //row counts down from (template_points_matrix->height() - 1) to 0
    {
        std::shared_ptr<TemplatePointsMatrixEntry> cur = template_points_matrix.get_row(row);
        if (cur == nullptr)
            continue;

        if (verbosity >= 10) {
            debug() << "----TemplatePointsMatrixEntry (" << cur->x << "," << cur->y << ")";
        }

        //store index of rightmost column that is mapped to this equivalence class
        auto cur_ind = (low) ? &(cur->low_index) : &(cur->high_index);
        *cur_ind = o_index;

        //get the multigrade list for this TemplatePointsMatrixEntry
        std::list<std::shared_ptr<Multigrade>>* mgrades = (low) ? &(cur->low_simplices) : &(cur->high_simplices);

        //sort the multigrades in lexicographical order
        mgrades->sort([](std::shared_ptr<Multigrade> a, std::shared_ptr<Multigrade> b) { return Multigrade::LexComparator(*a, *b); });

        //store map values for all simplices at these multigrades
        for (std::list<std::shared_ptr<Multigrade>>::iterator it = mgrades->begin(); it != mgrades->end(); ++it) {
            std::shared_ptr<Multigrade> mg = *it;
            if (verbosity >= 10) {
                debug() << "  multigrade (" << mg->x << "," << mg->y << ") has" << mg->num_cols << "simplices with last dim_index" << mg->simplex_index << "which will map to order_index" << o_index;
            }

            for (unsigned s = 0; s < mg->num_cols; s++) // simplex with dim_index (mg->simplex_index - s) has order_index o_index
            {
                simplex_order[mg->simplex_index - s] = o_index;
                o_index--;
            }
        }

        //if any simplices of the specified dimension were mapped to this equivalence class, then store information about this class
        if (*cur_ind != o_index) {
            if (low)
                lift_low.insert(std::pair<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>(*cur_ind, cur));
            else
                lift_high.insert(std::pair<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>(*cur_ind, cur));
        }
    } //end for(row > 0)

    return num_simplices;
} //end build_simplex_order()

//counts the number of transpositions that will happen if we cross an anchor and do vineyeard-updates
//this function DOES NOT MODIFY the xiSupportMatrix
unsigned long PersistenceUpdater::count_transpositions(std::shared_ptr<TemplatePointsMatrixEntry> anchor, bool from_below)
{
    //identify entries
    std::shared_ptr<TemplatePointsMatrixEntry> first = from_below ? anchor->down : anchor->left;
    std::shared_ptr<TemplatePointsMatrixEntry> second = from_below ? anchor->left : anchor->down;
    std::shared_ptr<TemplatePointsMatrixEntry> temp(new TemplatePointsMatrixEntry(anchor->left->x, anchor->down->y));

    //counters
    unsigned long count = 0;
    unsigned first_simplices_low = 0;
    unsigned first_simplices_high = 0;
    unsigned second_simplices_low = 0;
    unsigned second_simplices_high = 0;

    //count transpositions that occur when we separate out the grades that lift to the anchor from those that lift to second
    count_transpositions_from_separations(anchor, second, from_below, true, count, second_simplices_low);
    count_transpositions_from_separations(anchor, second, from_below, false, count, second_simplices_high);

    //count transpositions that occur when we separate out the grades that lift to GLB(first, second) from those that lift to first
    unsigned temp_simplices = 0;
    count_transpositions_from_separations(first, temp, from_below, true, count, temp_simplices);
    first_simplices_low = first->low_count - temp_simplices;
    temp_simplices = 0;
    count_transpositions_from_separations(first, temp, from_below, false, count, temp_simplices);
    first_simplices_high = first->high_count - temp_simplices;

    //count switches
    count += ((unsigned long)first_simplices_low) * ((unsigned long)second_simplices_low);
    count += ((unsigned long)first_simplices_high) * ((unsigned long)second_simplices_high);

    return count;
} //end count_transpositions()

//counts the number of transpositions that result from separations
//this function DOES NOT MODIFY the xiSupportMatrix
void PersistenceUpdater::count_transpositions_from_separations(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser, bool horiz, bool low, unsigned long& count_trans, unsigned& count_lesser)
{
    int gr_col = greater->low_index;
    int cur_col = gr_col;
    std::list<std::shared_ptr<Multigrade>> grades = low ? greater->low_simplices : greater->high_simplices;
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = grades.begin(); it != grades.end(); ++it) {
        std::shared_ptr<Multigrade> cur_grade = *it;
        if ((horiz && cur_grade->x > lesser->x) || (!horiz && cur_grade->y > lesser->y)) //then there will be transpositions from separations
        {
            count_trans += cur_grade->num_cols * (gr_col - cur_col);
            gr_col -= cur_grade->num_cols;
        } else
            count_lesser += cur_grade->num_cols;
        cur_col -= cur_grade->num_cols;
    }
} //end count_transpositions_from_separations()

//moves grades associated with TemplatePointsMatrixEntry greater, that come before TemplatePointsMatrixEntry lesser in R^2, so that they become associated with lesser
//  precondition: no grades lift to lesser (it has empty level sets under the lift map)
unsigned long PersistenceUpdater::split_grade_lists(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser, bool horiz)
{
    unsigned long swap_counter = 0;

    //low simplices
    int gr_col = greater->low_index;
    int cur_col = gr_col;
    std::list<std::shared_ptr<Multigrade>> grades = greater->low_simplices;
    greater->low_simplices.clear();
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = grades.begin(); it != grades.end(); ++it) {
        std::shared_ptr<Multigrade> cur_grade = *it;
        if ((horiz && cur_grade->x > lesser->x) || (!horiz && cur_grade->y > lesser->y)) //then this grade lifts to greater, so move columns to the right
        {
            if (cur_col != gr_col) //then we must move the columns
                swap_counter += move_low_columns(cur_col, cur_grade->num_cols, gr_col);

            greater->low_simplices.push_back(cur_grade);
            gr_col -= cur_grade->num_cols;
        } else //then this grade lifts to lesser, so update lift map, but no need to move columns
            lesser->low_simplices.push_back(cur_grade);

        cur_col -= cur_grade->num_cols;
    }
    lesser->low_index = gr_col;
    lesser->low_count = gr_col - cur_col;
    greater->low_count = greater->low_index - lesser->low_index;

    //high simplices
    gr_col = greater->high_index;
    cur_col = gr_col;
    grades = greater->high_simplices;
    greater->high_simplices.clear();
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = grades.begin(); it != grades.end(); ++it) {
        std::shared_ptr<Multigrade> cur_grade = *it;
        if ((horiz && cur_grade->x > lesser->x) || (!horiz && cur_grade->y > lesser->y)) //then this grade lifts to greater, so move columns to the right
        {
            if (cur_col != gr_col) //then we must move the columns
                swap_counter += move_high_columns(cur_col, cur_grade->num_cols, gr_col);

            greater->high_simplices.push_back(cur_grade);
            gr_col -= cur_grade->num_cols;
        } else //then this grade lifts to lesser, so update lift map, but no need to move columns
            lesser->high_simplices.push_back(cur_grade);

        cur_col -= cur_grade->num_cols;
    }
    lesser->high_index = gr_col;
    lesser->high_count = gr_col - cur_col;
    greater->high_count = greater->high_index - lesser->high_index;

    return swap_counter;
} //end split_grade_lists()

//splits grade lists and updates the permutation vectors, but does NOT do vineyard updates
void PersistenceUpdater::split_grade_lists_no_vineyards(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser, bool horiz)
{
    //STEP 1: update the lift map for all multigrades and store the current column index for each multigrade

    //first, low simpilices
    int gr_col = greater->low_index;
    int cur_col = gr_col;
    std::list<std::shared_ptr<Multigrade>> grades = greater->low_simplices;
    greater->low_simplices.clear(); ///this isn't so efficient...
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = grades.begin(); it != grades.end(); ++it) {
        std::shared_ptr<Multigrade> cur_grade = *it;
        cur_grade->simplex_index = cur_col;
        if ((horiz && cur_grade->x > lesser->x) || (!horiz && cur_grade->y > lesser->y)) //then this grade lifts to greater
        {
            greater->low_simplices.push_back(cur_grade);
            gr_col -= cur_grade->num_cols;
        } else //then this grade lifts to lesser
            lesser->low_simplices.push_back(cur_grade);

        cur_col -= cur_grade->num_cols;
    }
    lesser->low_index = gr_col;
    lesser->low_count = gr_col - cur_col;
    greater->low_count = greater->low_index - lesser->low_index;

    //now high simplices
    gr_col = greater->high_index;
    cur_col = gr_col;
    std::list<std::shared_ptr<Multigrade>> grades_h = greater->high_simplices;
    greater->high_simplices.clear(); ///this isn't so efficient...
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = grades_h.begin(); it != grades_h.end(); ++it) {
        std::shared_ptr<Multigrade> cur_grade = *it;
        cur_grade->simplex_index = cur_col;
        if ((horiz && cur_grade->x > lesser->x) || (!horiz && cur_grade->y > lesser->y)) //then this grade lifts to greater
        {
            greater->high_simplices.push_back(cur_grade);
            gr_col -= cur_grade->num_cols;
        } else //then this grade lifts to lesser
            lesser->high_simplices.push_back(cur_grade);

        cur_col -= cur_grade->num_cols;
    }
    lesser->high_index = gr_col;
    lesser->high_count = gr_col - cur_col;
    greater->high_count = greater->high_index - lesser->high_index;

    //STEP 2: traverse grades (backwards) in the new order and update the permutation vectors to reflect the new order on matrix columns

    //temporary data structures
    std::shared_ptr<TemplatePointsMatrixEntry> cur_entry = greater;
    int low_col = cur_entry->low_index;
    int high_col = cur_entry->high_index;

    //loop over xiMatrixEntrys
    while (true) //loop ends with break statement
    {
        //update positions of "low" simplices for this entry
        for (std::list<std::shared_ptr<Multigrade>>::iterator it = cur_entry->low_simplices.begin(); it != cur_entry->low_simplices.end(); ++it) {
            std::shared_ptr<Multigrade> cur_grade = *it;
            for (unsigned i = 0; i < cur_grade->num_cols; i++) {
                //column currently in position (cur_grade->simplex_index - i) has new position low_col
                unsigned original_position = inv_perm_low[cur_grade->simplex_index - i];
                perm_low[original_position] = low_col;
                low_col--;
            }
        }

        //update positions of "high" simplices for this entry
        for (std::list<std::shared_ptr<Multigrade>>::iterator it = cur_entry->high_simplices.begin(); it != cur_entry->high_simplices.end(); ++it) {
            std::shared_ptr<Multigrade> cur_grade = *it;
            for (unsigned i = 0; i < cur_grade->num_cols; i++) {
                //column currently in position (cur_grade->simplex_index - i) has new position high_col
                unsigned original_position = inv_perm_high[cur_grade->simplex_index - i];
                perm_high[original_position] = high_col;
                high_col--;
            }
        }

        //move to next entry
        if (cur_entry == greater)
            cur_entry = lesser;
        else
            break;
    } //end while

    //fix inverse permutation vectors -- is there a better way to do this?
    for (unsigned i = 0; i < perm_low.size(); i++)
        inv_perm_low[perm_low[i]] = i;
    for (unsigned i = 0; i < perm_high.size(); i++)
        inv_perm_high[perm_high[i]] = i;

} //end split_grade_lists_no_vineyards()

//moves all grades associated with TemplatePointsMatrixEntry lesser so that they become associated with TemplatePointsMatrixEntry greater
void PersistenceUpdater::merge_grade_lists(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser)
{
    //low simplices
    greater->low_simplices.splice(greater->low_simplices.end(), lesser->low_simplices);
    greater->low_count += lesser->low_count;
    lesser->low_count = 0;

    //high simplices
    greater->high_simplices.splice(greater->high_simplices.end(), lesser->high_simplices);
    greater->high_count += lesser->high_count;
    lesser->high_count = 0;
} //end merge_grade_lists()

//moves columns from an equivalence class given by std::shared_ptr<TemplatePointsMatrixEntry> first to their new positions after or among the columns in the equivalence class given by std::shared_ptr<TemplatePointsMatrixEntry> second
// the boolean argument indicates whether an anchor is being crossed from below (or from above)
// this version updates the permutation vectors required for the "reset" approach
unsigned long PersistenceUpdater::move_columns(std::shared_ptr<TemplatePointsMatrixEntry> first, std::shared_ptr<TemplatePointsMatrixEntry> second, bool from_below)
{
    //the following should never occur
    if (first->low_index + second->low_count != second->low_index || first->high_index + second->high_count != second->high_index) {
        throw std::runtime_error("PersistenceUpdater::move_columns(): swapping non-consecutive column blocks");
    }

    //get column indexes (so we know which columns to move)
    unsigned low_col = first->low_index; //rightmost column index of low simplices for the block that moves
    unsigned high_col = first->high_index; //rightmost column index of high simplices for the block that moves

    //set column indexes for the first class to their final position
    first->low_index = second->low_index;
    first->high_index = second->high_index;

    //initialize counter
    unsigned long swap_counter = 0;

    //move all "low" simplices for TemplatePointsMatrixEntry first (start with rightmost column, end with leftmost)
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = first->low_simplices.begin(); it != first->low_simplices.end();) //NOTE: iterator advances in loop
    {
        std::shared_ptr<Multigrade> cur_grade = *it;

        if ((from_below && cur_grade->x > second->x) || (!from_below && cur_grade->y > second->y))
        //then move columns at cur_grade past columns at TemplatePointsMatrixEntry second; lift map does not change ( lift : multigrades --> xiSupportElements )
        {
            swap_counter += move_low_columns(low_col, cur_grade->num_cols, second->low_index);
            second->low_index -= cur_grade->num_cols;
            ++it;
        } else //then cur_grade now lifts to TemplatePointsMatrixEntry second; columns don't move
        {
            //associate cur_grade with second
            second->insert_multigrade(cur_grade, true);
            it = first->low_simplices.erase(it); //NOTE: advances the iterator!!!

            //update column counts
            first->low_count -= cur_grade->num_cols;
            second->low_count += cur_grade->num_cols;
        }

        //update column index
        low_col -= cur_grade->num_cols;
    } //end "low" simplex loop

    //move all "high" simplices for TemplatePointsMatrixEntry first (start with rightmost column, end with leftmost)
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = first->high_simplices.begin(); it != first->high_simplices.end();) //NOTE: iterator advances in loop
    {
        std::shared_ptr<Multigrade> cur_grade = *it;

        // debug() << "  ====>>>> moving high simplices at grade (" << cur_grade->x << "," << cur_grade->y << ")";

        if ((from_below && cur_grade->x > second->x) || (!from_below && cur_grade->y > second->y))
        //then move columns at cur_grade past columns at TemplatePointsMatrixEntry second; lift map does not change ( lift : multigrades --> xiSupportElements )
        {
            swap_counter += move_high_columns(high_col, cur_grade->num_cols, second->high_index);
            second->high_index -= cur_grade->num_cols;
            ++it;
        } else //then cur_grade now lifts to TemplatePointsMatrixEntry second; columns don't move
        {
            // debug() << "====>>>> simplex at (" << cur_grade->x << "," << cur_grade->y << ") now lifts to (" << second->x << "," << second->y << ")";

            //associate cur_grade with second
            second->insert_multigrade(cur_grade, false);
            it = first->high_simplices.erase(it); //NOTE: advances the iterator!!!

            //update column counts
            first->high_count -= cur_grade->num_cols;
            second->high_count += cur_grade->num_cols;
        }

        //update column index
        high_col -= cur_grade->num_cols;
    } //end "high" simplex loop

    //the following should never occur
    if (second->low_index + first->low_count != first->low_index || second->high_index + first->high_count != first->high_index) {
        throw std::runtime_error("PersistenceUpdater::move_columns(): swap resulted in non-consecutive column blocks");
    }
    return swap_counter;
} //end move_columns()

//moves a block of n columns, the rightmost of which is column s, to a new position following column t (NOTE: assumes s <= t)
// this version maintains the permutation arrays required for the "reset" approach
unsigned long PersistenceUpdater::move_low_columns(int s, unsigned n, int t)
{
    // debug() << "   --Transpositions for low simplices: [" << s << "," << n << "," << t << "]:" << (n*(t-s)) << "total";

    //the following should never occur
    if (s > t) {
        throw std::runtime_error("PersistenceUpdater::move_low_columns(): illegal column move");
    }

    for (unsigned c = 0; c < n; c++) //move column that starts at s-c
    {
        for (int i = s; i < t; i++) {
            unsigned a = i - c;
            unsigned b = a + 1;

            //update the permutation vectors
            unsigned s = inv_perm_low[a];
            unsigned t = inv_perm_low[b];
            inv_perm_low[a] = t;
            inv_perm_low[b] = s;
            perm_low[t] = a;
            perm_low[s] = b;

            //now for the vineyards algorithm
            vineyard_update_low(a);

            /// TESTING ONLY - FOR CHECKING THAT D=RU
            //if (testing) {
            //    D_low->swap_columns(a, false);
            //    D_high->swap_rows(a, false);
            //}
        } //end for(i=...)
    } //end for(c=...)

    return n * (t - s);
} //end move_low_columns()

//moves a block of n columns, the rightmost of which is column s, to a new position following column t (NOTE: assumes s <= t)
// this version maintains the permutation arrays required for the "reset" approach
unsigned long PersistenceUpdater::move_high_columns(int s, unsigned n, int t)
{
    //    debug() << "   --Transpositions for high simplices: [" << s << "," << n << "," << t << "]:" << (n*(t-s)) << "total";

    //the following should never occur
    if (s > t) {
        throw std::runtime_error("PersistenceUpdater::move_high_columns(): illegal column move");
    }

    for (unsigned c = 0; c < n; c++) //move column that starts at s-c
    {
        for (int i = s; i < t; i++) {
            unsigned a = i - c;
            unsigned b = a + 1;

            //update the permutation vectors
            unsigned s = inv_perm_high[a];
            unsigned t = inv_perm_high[b];
            inv_perm_high[a] = t;
            inv_perm_high[b] = s;
            perm_high[t] = a;
            perm_high[s] = b;

            //now for the vineyards algorithm
            vineyard_update_high(a);

            /// TESTING ONLY - FOR CHECKING THAT D=RU
            //if (testing) {
            //    D_high->swap_columns(a, false);
            //}
        } //end for(i=...)
    } //end for(c=...)

    return n * (t - s);
} //end move_high_columns()

//performs a vineyard update corresponding to the transposition of columns a and (a + 1)
//  for LOW simplices
void PersistenceUpdater::vineyard_update_low(unsigned a)
{
    unsigned b = a + 1;

    bool a_pos = (R_low->low(a) == -1); //true iff simplex corresponding to column a is positive
    bool b_pos = (R_low->low(b) == -1); //true iff simplex corresponding to column b=a+1 is positive

    if (a_pos) //simplex a is positive (Vineyards paper - Cases 1 and 4)
    {
        if (b_pos) //simplex b is positive (Case 1)
        {
            //look for columns k and l in RH with low(k)=a, low(l)=b, and RH(a,l)=1 -- if these exist, then we must fix matrix RH following row/column swaps (Case 1.1)
            int k = R_high->find_low(a);
            int l = R_high->find_low(b);
            bool RHal = (l > -1 && R_high->entry(a, l)); //entry (a,l) in matrix RH

            //ensure that UL[a,b]=0
            U_low->clear(a, b);

            //transpose rows and columns (don't need to swap columns of RL, because these columns are zero)
            U_low->swap_columns(a);
            U_low->swap_rows(a);

            //swap rows, and fix RH if necessary
            if (k > -1 && RHal) //case 1.1
            {
                if (k < l) {
                    R_high->swap_rows(a, true); //in this case, low entries change
                    R_high->add_column(k, l);
                    U_high->add_row(l, k);
                } else {
                    R_high->swap_rows(a, false); //in this case, low entries do not change
                    R_high->add_column(l, k);
                    U_high->add_row(k, l);
                }
            } else
                R_high->swap_rows(a, !RHal); //in this case, only necessary to update low entries if RH(a,l)=0 or if column l does not exist
        } else //simplex b is negative (Case 4)
        {
            //ensure that UL[a,b]=0
            U_low->clear(a, b);

            //transpose rows and columns and update low arrays
            R_low->swap_columns(a, true);
            R_high->swap_rows(a, true);
            U_low->swap_columns(a);
            U_low->swap_rows(a);
        }
    } else //simplex a is negative (Vineyards paper - Cases 2 and 3)
    {
        if (b_pos) //simplex b is positive (Case 3)
        {
            //look for column l in RH with low(l)=b and RH(a,l)=1
            int l = R_high->find_low(b);
            bool RHal = (l > -1 && R_high->entry(a, l)); //entry (a,l) in matrix RH

            //transpose rows of R; update low array if necessary
            R_high->swap_rows(a, !RHal);

            if (U_low->entry(a, b)) //case 3.1 -- here, R = RWPW, so no further action required on R
            {
                U_low->add_row(b, a);
                U_low->swap_rows(a);
                U_low->add_row(b, a);
            } else //case 3.2
            {
                R_low->swap_columns(a, true);
                U_low->swap_rows(a);
            }
        } else //simplex b is negative (Case 2)
        {
            //transpose rows of R
            R_high->swap_rows(a, false); //neither of these rows contain lowest 1's in any column

            if (U_low->entry(a, b)) //case 2.1
            {
                U_low->add_row(b, a); //so that U will remain upper-triangular
                U_low->swap_rows(a); //swap rows of U

                if (R_low->low(a) < R_low->low(b)) //case 2.1.1
                {
                    R_low->add_column(a, b); //necessary due to the row addition on U; this doesn't change low entries
                    R_low->swap_columns(a, true); //now swap columns of R and update low entries
                } else //case 2.1.2
                {
                    R_low->add_column(a, b); //necessary due to the row addition on U; this doesn't change low entries
                    R_low->swap_columns(a, false); //now swap columns of R but DO NOT update low entries
                    R_low->add_column(a, b); //restore R to reduced form; low entries now same as they were initially
                    U_low->add_row(b, a); //necessary due to column addition on R
                }
            } else //case 2.2
            {
                R_low->swap_columns(a, true); //swap columns of R and update low entries
                U_low->swap_rows(a); //swap rows of U
            }
        }

        //finally, for cases 2 and 3, transpose columns of U
        U_low->swap_columns(a);
    }
} //end vineyard_update_low()

//performs a vineyard update corresponding to the transposition of columns a and (a + 1)
//  for HIGH simplices
void PersistenceUpdater::vineyard_update_high(unsigned a)
{
    unsigned b = a + 1;

    bool a_pos = (R_high->low(a) == -1); //true iff simplex corresponding to column a is positive
    bool b_pos = (R_high->low(b) == -1); //true iff simplex corresponding to column b is positive

    if (a_pos) //simplex a is positive, so its column is zero, and the fix is easy  (Vineyards paper - Cases 1 and 4)
    {
        if (!b_pos) //only have to swap columns of R if column b is nonzero
            R_high->swap_columns(a, true);

        //ensure that UL[a,b]=0
        U_high->clear(a, b);

        //transpose rows and columns of U
        U_high->swap_columns(a);
        U_high->swap_rows(a);

        //done -- we don't care about the ROWS corresponding to simplices a and b, because we don't care about the boundaries of (d+2)-simplices
    } else //simplex a is negative (Vineyards paper - Cases 2 and 3)
    {
        if (b_pos) //simplex b is positive (Case 3)
        {
            if (U_high->entry(a, b)) //case 3.1 -- here, R = RWPW, so no further action required on R
            {
                U_high->add_row(b, a);
                U_high->swap_rows(a);
                U_high->add_row(b, a);
            } else //case 3.2
            {
                R_high->swap_columns(a, true);
                U_high->swap_rows(a);
            }
        } else //simplex b is negative (Case 2)
        {
            if (U_high->entry(a, b)) //case 2.1
            {
                U_high->add_row(b, a); //so that U will remain upper-triangular
                U_high->swap_rows(a); //swap rows of U

                if (R_high->low(a) < R_high->low(b)) //case 2.1.1
                {
                    R_high->add_column(a, b); //necessary due to the row addition on U; this doesn't change low entries
                    R_high->swap_columns(a, true); //now swap columns of R and update low entries
                } else //case 2.1.2
                {
                    R_high->add_column(a, b); //necessary due to the row addition on U; this doesn't change low entries
                    R_high->swap_columns(a, false); //now swap columns of R but DO NOT update low entries
                    R_high->add_column(a, b); //restore R to reduced form; low entries now same as they were initially
                    U_high->add_row(b, a); //necessary due to column addition on R
                }
            } else //case 2.2
            {
                R_high->swap_columns(a, true); //swap columns and update low entries
                U_high->swap_rows(a); //swap rows of U
            }
        }

        //finally, for Cases 2 and 3, transpose columns of U
        U_high->swap_columns(a);
    }
} //end vineyard update_high()

//swaps two blocks of columns by updating the total order on columns, then rebuilding the matrices and computing a new RU-decomposition
void PersistenceUpdater::update_order_and_reset_matrices(std::shared_ptr<TemplatePointsMatrixEntry> first, std::shared_ptr<TemplatePointsMatrixEntry> second, bool from_below, MapMatrix_Perm* RL_initial, MapMatrix_Perm* RH_initial)
{
    //STEP 1: update the lift map for all multigrades and store the current column index for each multigrade

    //store current column index for each multigrade that lifts to TemplatePointsMatrixEntry second
    int low_col = second->low_index;
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = second->low_simplices.begin(); it != second->low_simplices.end(); ++it) //starts with rightmost column, ends with leftmost
    {
        (*it)->simplex_index = low_col;
        low_col -= (*it)->num_cols;
    }
    int high_col = second->high_index;
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = second->high_simplices.begin(); it != second->high_simplices.end(); ++it) //starts with rightmost column, ends with leftmost
    {
        (*it)->simplex_index = high_col;
        high_col -= (*it)->num_cols;
    }

    //get column indexes for the first equivalence class
    low_col = first->low_index; //rightmost column index of low simplices for the equivalence class to move
    high_col = first->high_index; //rightmost column index of high simplices for the equivalence class to move

    //set column indexes for the first class to their final position
    first->low_index = second->low_index;
    first->high_index = second->high_index;

    //store current column index and update the lift map for each multigrade that lifts to TemplatePointsMatrixEntry first
    //"low" simplices (start with rightmost column, end with leftmost)
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = first->low_simplices.begin(); it != first->low_simplices.end();) //NOTE: iterator advances in loop
    {
        std::shared_ptr<Multigrade> cur_grade = *it;

        //remember current position of this grade
        cur_grade->simplex_index = low_col;

        if ((from_below && cur_grade->x > second->x) || (!from_below && cur_grade->y > second->y))
        //then move columns at cur_grade past columns at TemplatePointsMatrixEntry second; lift map does not change ( lift : multigrades --> xiSupportElements )
        {
            second->low_index -= cur_grade->num_cols;
            ++it;
        } else //then cur_grade now lifts to TemplatePointsMatrixEntry second; columns don't move
        {
            //associate cur_grade with second
            second->insert_multigrade(cur_grade, true);
            it = first->low_simplices.erase(it); //NOTE: advances the iterator!!!

            //update column counts
            first->low_count -= cur_grade->num_cols;
            second->low_count += cur_grade->num_cols;
        }

        //update column index
        low_col -= cur_grade->num_cols;
    } //end "low" simplex loop

    //"high" simplices (start with rightmost column, end with leftmost)
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = first->high_simplices.begin(); it != first->high_simplices.end();) //NOTE: iterator advances in loop
    {
        std::shared_ptr<Multigrade> cur_grade = *it;

        //remember current position of this grade
        cur_grade->simplex_index = high_col;

        if ((from_below && cur_grade->x > second->x) || (!from_below && cur_grade->y > second->y))
        //then move columns at cur_grade past columns at TemplatePointsMatrixEntry second; lift map does not change ( lift : multigrades --> xiSupportElements )
        {
            second->high_index -= cur_grade->num_cols;
            ++it;
        } else //then cur_grade now lifts to TemplatePointsMatrixEntry second; columns don't move
        {
            //associate cur_grade with target
            second->insert_multigrade(cur_grade, false);
            it = first->high_simplices.erase(it); //NOTE: advances the iterator!!!

            //update column counts
            first->high_count -= cur_grade->num_cols;
            second->high_count += cur_grade->num_cols;
        }

        //update column index
        high_col -= cur_grade->num_cols;
    } //end "high" simplex loop

    //STEP 2: traverse grades (backwards) in the new order and update the permutation vectors to reflect the new order on matrix columns

    //temporary data structures
    std::shared_ptr<TemplatePointsMatrixEntry> cur_entry = first;
    low_col = first->low_index;
    high_col = first->high_index;

    //loop over xiMatrixEntrys
    while (first != second) {
        //update positions of "low" simplices for this entry
        for (std::list<std::shared_ptr<Multigrade>>::iterator it = cur_entry->low_simplices.begin(); it != cur_entry->low_simplices.end(); ++it) {
            std::shared_ptr<Multigrade> cur_grade = *it;
            for (unsigned i = 0; i < cur_grade->num_cols; i++) {
                //column currently in position (cur_grade->simplex_index - i) has new position low_col
                unsigned original_position = inv_perm_low[cur_grade->simplex_index - i];
                perm_low[original_position] = low_col;
                low_col--;
            }
        }

        //update positions of "high" simplices for this entry
        for (std::list<std::shared_ptr<Multigrade>>::iterator it = cur_entry->high_simplices.begin(); it != cur_entry->high_simplices.end(); ++it) {
            std::shared_ptr<Multigrade> cur_grade = *it;
            for (unsigned i = 0; i < cur_grade->num_cols; i++) {
                //column currently in position (cur_grade->simplex_index - i) has new position high_col
                unsigned original_position = inv_perm_high[cur_grade->simplex_index - i];
                perm_high[original_position] = high_col;
                high_col--;
            }
        }

        //move to next entry
        if (cur_entry == first)
            cur_entry = second;
        else
            break;
    } //end while

    //fix inverse permutation vectors -- is there a better way to do this?
    for (unsigned i = 0; i < perm_low.size(); i++)
        inv_perm_low[perm_low[i]] = i;
    for (unsigned i = 0; i < perm_high.size(); i++)
        inv_perm_high[perm_high[i]] = i;

    //STEP 3: re-build the matrix R based on the new order

    R_low->rebuild(RL_initial, perm_low);
    R_high->rebuild(RH_initial, perm_high, perm_low);

    //STEP 4: compute the new RU-decomposition

    ///TODO: should I avoid deleting and reallocating matrix U?
    delete U_low;
    U_low = R_low->decompose_RU();
    delete U_high;
    U_high = R_high->decompose_RU();

} //end update_order_and_reset_matrices()

//updates the total order on columns, rebuilds the matrices, and computing a new RU-decomposition for a NON-STRICT anchor
void PersistenceUpdater::update_order_and_reset_matrices(MapMatrix_Perm* RL_initial, MapMatrix_Perm* RH_initial)
{
    //anything to do here?????

    //re-build the matrix R based on the new order
    R_low->rebuild(RL_initial, perm_low);
    R_high->rebuild(RH_initial, perm_high, perm_low);

    //compute the new RU-decomposition
    ///TODO: should I avoid deleting and reallocating matrix U?
    delete U_low;
    U_low = R_low->decompose_RU();
    delete U_high;
    U_high = R_high->decompose_RU();

} //end update_order_and_reset_matrices()

//swaps two blocks of simplices in the total order, and returns the number of transpositions that would be performed on the matrix columns if we were doing vineyard updates
void PersistenceUpdater::count_switches_and_separations(std::shared_ptr<TemplatePointsMatrixEntry> at_anchor, bool from_below, unsigned long& switches, unsigned long& seps)
{
    //identify entries
    std::shared_ptr<TemplatePointsMatrixEntry> first = from_below ? at_anchor->down : at_anchor->left;
    std::shared_ptr<TemplatePointsMatrixEntry> second = from_below ? at_anchor->left : at_anchor->down;
    std::shared_ptr<TemplatePointsMatrixEntry> temp(new TemplatePointsMatrixEntry(at_anchor->left->x, at_anchor->down->y)); //temporary entry for holding grades that come before BOTH first and second

    //separate out the grades that lift to anchor from those that lift to second
    do_separations(at_anchor, second, from_below);
    seps += second->low_count * at_anchor->low_count + second->high_count * at_anchor->high_count;
    do_separations(first, temp, from_below);
    seps += temp->low_count * first->low_count + temp->high_count * first->high_count;

    //count switches
    int i = first->low_index;
    first->low_index = second->low_index;
    second->low_index = i;
    i = first->high_index;
    first->high_index = second->high_index;
    second->high_index = i;
    switches += first->low_count * second->low_count + first->high_count * second->high_count;

    //count final separations
    seps += first->low_count * at_anchor->low_count + first->high_count * at_anchor->high_count;
    merge_grade_lists(at_anchor, first);
    seps += temp->low_count * second->low_count + temp->high_count * second->high_count;
    merge_grade_lists(second, temp);
} //end count_switches_and_separations()

//used by the previous function to split grade lists at each anchor crossing
void PersistenceUpdater::do_separations(std::shared_ptr<TemplatePointsMatrixEntry> greater, std::shared_ptr<TemplatePointsMatrixEntry> lesser, bool horiz)
{
    //first, low simpilicse
    int gr_col = greater->low_index;
    int cur_col = gr_col;
    std::list<std::shared_ptr<Multigrade>> grades = greater->low_simplices;
    greater->low_simplices.clear(); ///this isn't so efficient...
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = grades.begin(); it != grades.end(); ++it) {
        std::shared_ptr<Multigrade> cur_grade = *it;
        if ((horiz && cur_grade->x > lesser->x) || (!horiz && cur_grade->y > lesser->y)) //then this grade lifts to greater
        {
            greater->low_simplices.push_back(cur_grade);
            gr_col -= cur_grade->num_cols;
        } else //then this grade lifts to lesser
            lesser->low_simplices.push_back(cur_grade);

        cur_col -= cur_grade->num_cols;
    }
    lesser->low_index = gr_col;
    lesser->low_count = gr_col - cur_col;
    greater->low_count = greater->low_index - lesser->low_index;

    //now high simplices
    gr_col = greater->high_index;
    cur_col = gr_col;
    std::list<std::shared_ptr<Multigrade>> grades_h = greater->high_simplices;
    greater->high_simplices.clear(); ///this isn't so efficient...
    for (std::list<std::shared_ptr<Multigrade>>::iterator it = grades_h.begin(); it != grades_h.end(); ++it) {
        std::shared_ptr<Multigrade> cur_grade = *it;
        if ((horiz && cur_grade->x > lesser->x) || (!horiz && cur_grade->y > lesser->y)) //then this grade lifts to greater
        {
            greater->high_simplices.push_back(cur_grade);
            gr_col -= cur_grade->num_cols;
        } else //then this grade lifts to lesser
            lesser->high_simplices.push_back(cur_grade);

        cur_col -= cur_grade->num_cols;
    }
    lesser->high_index = gr_col;
    lesser->high_count = gr_col - cur_col;
    greater->high_count = greater->high_index - lesser->high_index;
} //end do_separations

//removes entries corresponding to TemplatePointsMatrixEntry head from lift_low and lift_high
void PersistenceUpdater::remove_lift_entries(std::shared_ptr<TemplatePointsMatrixEntry> entry)
{
    if (verbosity >= 10) {
        debug() << "    ----removing partition entries for TemplatePointsMatrixEntry" << entry->index << "(" << entry->low_index << ";" << entry->high_index << ")";
    }

    //low simplices
    std::map<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it1 = lift_low.find(entry->low_index);
    if (it1 != lift_low.end() && it1->second == entry)
        lift_low.erase(it1);

    //high simplices
    std::map<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it2 = lift_high.find(entry->high_index);
    if (it2 != lift_high.end() && it2->second == entry)
        lift_high.erase(it2);

} //end remove_lift_entries()

//if the equivalence class corresponding to TemplatePointsMatrixEntry head has nonempty sets of "low" or "high" simplices, then this function creates the appropriate entries in lift_low and lift_high
void PersistenceUpdater::add_lift_entries(std::shared_ptr<TemplatePointsMatrixEntry> entry)
{
    if (verbosity >= 10) {
        debug() << "    ----adding partition entries for TemplatePointsMatrixEntry" << entry->index << "(" << entry->low_index << ";" << entry->high_index << ")";
    }

    //low simplices
    if (entry->low_count > 0)
        lift_low.insert(std::pair<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>(entry->low_index, entry));

    //high simplices
    if (entry->high_count > 0)
        lift_high.insert(std::pair<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>(entry->high_index, entry));
} //end add_lift_entries()

//stores a barcode template in a 2-cell of the arrangement
///TODO: IMPROVE THIS!!! (store previous barcode at the simplicial level, and only examine columns that were modified in the recent update)
/// Is there a better way to handle endpoints at infinity?
void PersistenceUpdater::store_barcode_template(std::shared_ptr<Face> cell)
{
    Debug qd = debug(true);
    if (verbosity >= 6) {
        qd << "  -- barcode: ";
    }

    //mark this cell as visited
    cell->mark_as_visited();

    //get a reference to the barcode template object
    BarcodeTemplate& dbc = cell->get_barcode();

    //loop over all zero-columns in matrix R_low
    for (unsigned c = 0; c < R_low->width(); c++) {
        if (R_low->col_is_empty(c)) //then simplex corresponding to column c is positive
        {
            //find index of template point corresponding to simplex c
            std::map<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>::iterator tp1 = lift_low.lower_bound(c);
            unsigned a = (tp1 != lift_low.end()) ? tp1->second->index : -1; //index is -1 iff the simplex maps to infinity

            //is simplex s paired?
            int s = R_high->find_low(c);
            if (s != -1) //then simplex c is paired with negative simplex s
            {
                //find index of xi support point corresponding to simplex s
                std::map<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>::iterator tp2 = lift_high.lower_bound(s);
                unsigned b = (tp2 != lift_high.end()) ? tp2->second->index : -1; //index is -1 iff the simplex maps to infinity

                if (a != b) //then we have a bar of positive length
                {
                    dbc.add_bar(a, b);
                    if (verbosity >= 6) {
                        qd << "(" << c << "," << s << ")-->(" << a << "," << b << ") ";
                    }
                }
            } else //then simplex c generates an essential cycle
            {
                dbc.add_bar(a, -1); //b = -1 = MAX_UNSIGNED indicates this is an essential cycle
                if (verbosity >= 6) {
                    qd << c << "-->" << a << " ";
                }
            }
        }
    }
    if (verbosity >= 6) {
        qd << "\n ";
    }
} //end store_barcode_template()

//chooses an initial threshold by timing vineyard updates corresponding to random transpositions
void PersistenceUpdater::choose_initial_threshold(unsigned decomp_time, unsigned long & num_trans, unsigned & trans_time, unsigned long & threshold)
{
    if (verbosity >= 4) {
        debug() << "RANDOM VINEYARD UPDATES TO CHOOSE THE INITIAL THRESHOLD";
    }

    //other data structures
    unsigned num_cols = R_low->width() + R_high->width();
    std::list<unsigned> trans_list;

    //avoid trivial cases
    if (num_cols <= 3) { //if neither the low or high matrix has at least 2 columns, then we can't do transpositions
        threshold = 1000;
        return;
    }

    //determine the time for which we will do transpositions
    int runtime = decomp_time / 20;
    if (runtime < 100)
        runtime = 100; //run for at least 100 milliseconds

    //start the timer
    Timer timer;

    //do transpositions
    if (verbosity >= 8) {
        debug() << "  -->Doing some random vineyard updates...";
    }
    while ( (timer.elapsed() < runtime || trans_list.size() == 0) &&
        (timer.elapsed() < 5 || trans_list.size() < 5000) ) //do a transposition
    {
        unsigned rand_col = rand() % (num_cols - 1); //random integer in {0, 1, ..., num_cols - 2}

        if (rand_col + 1 < R_low->width()) //then transpose LOW columns rand_col and (rand_col + 1)
        {
            vineyard_update_low(rand_col);
            trans_list.push_back(rand_col);
        } else if (R_low->width() <= rand_col) //then transpose HIGH columns rand_col and (rand_col + 1)
        {
            vineyard_update_high(rand_col - R_low->width());
            trans_list.push_back(rand_col);
        }
        //note that if (rand_col + 1 == R_low->width()), then we don't do anything, since rand_col is a "low" simplex but rand_col + 1 is a "high" simplex, and these cannot swap
    }

    //do the inverse transpositions
    if (verbosity >= 8) {
        debug() << "  -->Undoing the random vineyard updates...";
    }
    for (auto rit = trans_list.rbegin(); rit != trans_list.rend(); ++rit) {
        auto col = *rit;
        if (col < R_low->width()) {
            vineyard_update_low(col);
        } else {
            vineyard_update_high(col - R_low->width());
        }
    }

    //record the time and number of transpositions
    trans_time = timer.elapsed();
    num_trans = 2 * trans_list.size();

    if (verbosity >= 8) {
        debug() << "  -->Did" << num_trans << "vineyard updates in" << trans_time << "milliseconds.";
    }

    //compute the threshold
    threshold = (unsigned long)(((double)num_trans / (double)trans_time) * decomp_time);
} //end choose_initial_threshold()

///TESTING ONLY
/// functions to check that D=RU
// void PersistenceUpdater::check_low_matrix(MapMatrix_Perm* RL, MapMatrix_RowPriority_Perm* UL)
// {
//     bool err_low = false;
//     for (unsigned row = 0; row < D_low->height(); row++) {
//         for (unsigned col = 0; col < D_low->width(); col++) {
//             bool temp = false;
//             for (unsigned e = 0; e < D_low->width(); e++)
//                 temp = (temp != (RL->entry(row, e) && UL->entry(e, col)));
//             if (temp != D_low->entry(row, col))
//                 err_low = true;
//         }
//     }
//     if (err_low) {
//         debug() << "====>>>> MATRIX ERROR (low) AT THIS STEP!";
//         //        debug() << "  Reduced matrix for low simplices:";
//         //        RL->print();
//         //        debug() << "  Matrix U for low simplices:";
//         //        UL->print();
//         //        debug() << "  Matrix D for low simplices:";
//         //        D_low->print();
//     } else
//         debug() << "low matrix ok";
// }

// void PersistenceUpdater::check_high_matrix(MapMatrix_Perm* RH, MapMatrix_RowPriority_Perm* UH)
// {
//     bool err_high = false;
//     for (unsigned row = 0; row < D_high->height(); row++) {
//         for (unsigned col = 0; col < D_high->width(); col++) {
//             bool temp = false;
//             for (unsigned e = 0; e < D_high->width(); e++)
//                 temp = (temp != (RH->entry(row, e) && UH->entry(e, col)));
//             if (temp != D_high->entry(row, col))
//                 err_high = true;
//         }
//     }
//     if (err_high)
//         debug() << "====>>>> MATRIX ERROR (high) AT THIS STEP!\n";
//     else
//         debug() << "high matrix ok";
// }

void PersistenceUpdater::print_perms(Perm& per, Perm& inv)
{
    debug(true) << "  permutation: ";
    for (unsigned i = 0; i < per.size(); i++)
        debug(true) << per[i] << " ";
    debug(true) << "\n  inverse permutation: ";
    for (unsigned i = 0; i < inv.size(); i++)
        debug(true) << inv[i] << " ";
}

void PersistenceUpdater::print_high_partition()
{
    debug(true) << "  high partition: ";
    for (std::map<unsigned, std::shared_ptr<TemplatePointsMatrixEntry>>::iterator it = lift_high.begin(); it != lift_high.end(); ++it)
        debug(true) << it->first << "->" << it->second->index << ", ";
}
