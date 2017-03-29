//#include "Graph.h"
//#include <iostream>
//#include <QtWidgets>
//#include <QApplication>
//#include "qnemainwindow.h"
//#include "dendrogram_data.h"
#include "Graph.h"
#include "time_root.h"
//#include <unordered_map>
//#include <unordered_set>
//#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>


#ifndef DENDROGRAM_DATA_H
#define DENDROGRAM_DATA_H

class Dendrogram_data
{
public:
    // A structure to represent a subset for union-find

    struct subset
    {
        int parent;
        int rank;
    };
    typedef subset Subset;

    /*Dendrogram_data::Dendrogram_data() :
    { }

    Dendrogram_data::~Dendrogram_data()
    { }*/

    // A utility function to find set of an element i
    // (uses path compression technique)
    static int Find( Subset subsets[], int i )
    {
        // find root and make root as parent of i (path compression)
        if (subsets[i].parent != i)
            subsets[i].parent = Find(subsets, subsets[i].parent);

        return subsets[i].parent;
    }

    // A function that does union of two sets of x and y by rank
    // returns (lower root, upper root)
    static std::pair<int,int> Union( Subset subsets[], int x, int y )
    {
        int xroot = Find( subsets, x );
        int yroot = Find( subsets, y );

        // Attach smaller rank tree under root of high rank tree
        // (Union by Rank)
        if ( subsets[xroot].rank < subsets[yroot].rank )
        {
            subsets[xroot].parent = yroot;
            return std::pair<int,int>(xroot,yroot);
        }
        else if (subsets[xroot].rank > subsets[yroot].rank)
        {
            subsets[yroot].parent = xroot;
            return std::pair<int,int>(yroot,xroot);
        }
        else
        {
            // If ranks are same, then make one as root and increment
            // its rank by one
            /*
            subsets[yroot].parent = xroot;
            subsets[xroot].rank++;
            return std::pair<int,int>(yroot,xroot);*/


            std::cout << "ranks are the same" << std::endl;
            int big = std::max(xroot, yroot);
            int small = std::min(xroot, yroot);
            /*
            subsets[yroot].parent = xroot;
            subsets[xroot].rank++;
            return std::pair<int,int>(yroot,xroot);*/
            subsets[small].parent = big;
            subsets[big].rank++;
            return std::pair<int,int>(small,big);


        }
    }

    static void print_upper_roots(std::unordered_set<int> upper_roots)
    {
        std::stringstream ss;
        for (int x : upper_roots)
        {
            ss << x << ",";
        }
        std::cout << ss.str();
    }

    static void print_time_root_to_tr(boost::unordered::unordered_map<std::pair<double,int>, Time_root>& time_root_to_tr,
                                          bool with_comp_labels)
    {
        for(auto t_r_tr : time_root_to_tr)
        {
            std::cout << "===============" << std::endl;
            Time_root::print_tr(t_r_tr.second, with_comp_labels);
            std::cout << "children: " << std::endl;
            Time_root::print_children(t_r_tr.second.children, with_comp_labels);
        }
    }

    // The main function to compute dendrogram using Kruskal's algorithm
    static void compute_dendrogram( DenseGRAPH<Edge>* graph,
                     std::vector<EdgePtr>& mst,
                     boost::unordered::unordered_map<std::pair<double,int>, Time_root>& time_root_to_tr,
                     Time_root& last_upper_tr,
                     boost::unordered::unordered_map<unsigned, double>& vertex_appearance,
                     boost::unordered::unordered_map<unsigned, std::pair<double,double>>& vertex_bigrade,
                     boost::unordered::unordered_map<double, std::pair<double,double>>& appearance_to_bigrade,
                     std::vector<Time_root>& upper_trs)
    {
        const int V = graph->V();

        // Allocate memory for creating V ssubsets
        std::unique_ptr<subset[]> subsets( new subset[ V ]() );
        //std::unordered_map<int, double> root_to_prev_time; // maps a root to the last time a node was made for it

        // Create V subsets with single elements
        std::cout << "Create V subsets with single elements" << std::endl;
        for (auto v_grade : vertex_appearance)
        {
            int v = v_grade.first;
            subsets[v].parent = v;
            subsets[v].rank = 0;
            graph->insert( EdgePtr( new Edge( v, v, vertex_appearance[v], vertex_bigrade[v] ) ) );
            std::cout << "(v, vertex_appearance) = (" << v << " , "
                      << vertex_appearance[v] << " , ("
                      << vertex_bigrade[v].first << " , " << vertex_bigrade[v].second << "))"
                      << std::endl;
        }
        // insert dummy edge with inifinite weight at the end to update last wave new_heads
        // MAGIC NUMBER
        //graph->insert( EdgePtr( new Edge( 0, 0, std::numeric_limits<double>::infinity() ) ) );
        // Sort edges in non-decreasing order of their weight
        std::vector<EdgePtr> edges;
        graph->edges( edges );
        std::sort(edges.begin(), edges.end(), EdgeSort() );

        std::cout << "print edges" << std::endl;
        for ( std::vector<EdgePtr>::iterator it = edges.begin(); it != edges.end(); ++it )
        {
            EdgePtr edgePtr = *it;
            std::cout << "(v,w,wt,bigrade) = (" << it->get()->v << " , "
                      << it->get()->w << " , "
                      << it->get()->wt << " , ("
                      << appearance_to_bigrade[it->get()->wt].first << ","
                      << appearance_to_bigrade[it->get()->wt].second << "))"
                      << std::endl;
        }

        double current_time = edges.front()->wt;
        //std::pair<double,double> current_bigrade = edges.front()->bigrade;
        std::pair<double,double> current_bigrade = appearance_to_bigrade[current_time];
        std::unordered_set<int> new_heads;
        std::vector<EdgePtr> new_edges;
        std::vector<int> lower_roots;
        std::unordered_set<int> upper_roots;
        //boost::unordered::unordered_map<int, std::vector<int>> new_head_to_children;
        //boost::unordered::unordered_map<int, Time_root> V_to_prev_tr; // maps a 0-simplex v to the last time root that represents v's component
        boost::unordered::unordered_map<int, std::shared_ptr<std::pair<double,int>>> V_to_prev_tr;


        //std::vector<int> all_roots_at_next_time; // keep track of all time roots formed at next_time for num_leaves updating

        std::vector<EdgePtr>::size_type size = edges.size();

        // for ( std::vector<EdgePtr>::iterator it = edges.begin(); it != edges.end(); ++it )
        for (std::vector<EdgePtr>::size_type i = 0; i < size; ++i)
        {
            //EdgePtr edgePtr = *it;
            EdgePtr edgePtr = edges[i];

            //std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
            /*std::cout << edgePtr->v << " "
                      << edgePtr->w << " "
                      << edgePtr->wt << std::endl;*/
            //bool is_last_edge = (it + 1 == edges.end());
            bool is_last_edge = (i == size - 1);
            //std::cout << "is_last_edge = " << is_last_edge << std::endl;
            if (is_last_edge && edgePtr->wt != std::numeric_limits<double>::infinity())
            {
                int v = edgePtr->v;
                int w = edgePtr->w;
                EdgePtr dummy_edge ( new Edge( v, w, std::numeric_limits<double>::infinity() ) );
                edges.push_back(dummy_edge);
                ++size;
                //bool edges_has_one_more = (it + 1 != edges.end());
                //std::cout << "edges_has_one_more = " << edges_has_one_more << std::endl;
            }

            if (edgePtr->wt > current_time)
            {
                std::cout << "wt > current_time" << std::endl;
                std::cout << "wt = " << edgePtr->wt << std::endl;
                std::cout << "current_time = " << current_time << std::endl;
                std::cout << "current_bigrade = (" << current_bigrade.first << "," << current_bigrade.second << ")" << std::endl;
                std::cout << "start new_heads: " << std::endl;
                std::cout << "upper_roots before erase" << std::endl;
                print_upper_roots(upper_roots);
                std::cout << std::endl;
                for (int v : lower_roots)
                {
                    new_heads.erase(v);

                    if (upper_roots.find(v) != upper_roots.end())
                    {
                        upper_roots.erase(v);
                    }
                }
                std::cout << "upper_roots after erase:" << std::endl;
                print_upper_roots(upper_roots);
                std::cout << std::endl;
                for ( int v : new_heads )
                {
                    std::cout << v << std::endl;
                    Time_root tr_v;
                    tr_v.t = current_time;
                    tr_v.r = v;
                    tr_v.bigrade = current_bigrade;
                    tr_v.num_children = 0;
                    tr_v.num_leaves = 0;
                    time_root_to_tr[std::pair<double,int>(current_time,v)] = tr_v;
                }
                std::cout << "end new heads" << std::endl;
                std::cout << " ^^^^^^^^^^^^^^^^ upper_roots: " << std::endl;
                print_upper_roots(upper_roots);

                std::cout << "------ updating new_edges ------" << std::endl;
                std::unordered_set<int> V_that_have_become_children;
                for ( std::vector<EdgePtr>::iterator it_new = new_edges.begin(); it_new != new_edges.end(); ++it_new )
                {
                    std::cout << "*************************" << std::endl;
                    EdgePtr edgePtr = *it_new;
                    int v = edgePtr->v;
                    int w = edgePtr->w;
                    /*std::cout << it_new->get()->v << " "
                              << it_new->get()->w << " "
                              << it_new->get()->wt << std::endl;*/
                    //int head_v = Find (subsets.get(), edgePtr->v);
                    int h_v = Find (subsets.get(), edgePtr->v);
                    int h_w = Find (subsets.get(), edgePtr->w);
                    std::cout << "(v,w) = (" << v << "," << w << ")" << std::endl;
                    std::cout << "(h_v,h_w) = (" << h_v << "," << h_w << std::endl;
                    Time_root& head_tr = time_root_to_tr[std::pair<double,int>(current_time,h_v)];

                    //if (root_to_prev_time.find(v) != root_to_prev_time.end()
                    //        && root_to_prev_time[v] < current_time)
                    if (V_to_prev_tr.find(v) != V_to_prev_tr.end())
                    {
                        std::cout << "case 1" << std::endl;
                        std::cout << "head_tr: ";
                        Time_root::print_tr(head_tr);
                        //Time_root v_tr = time_root_to_tr[std::pair<double,int>(root_to_prev_time[v],v)];
                        //Time_root& v_tr = time_root_to_tr[*V_to_prev_tr[v]];
                        // follow pointers in V_to_prev_tr until you get to current root
                        // WARNING : MAY INCUR LOG(N) COMPLEXITY BUT ONLY DONE FOR INITIAL TEMPLATE SO PROBABLY OKAY
                        int r = V_to_prev_tr[v]->second;
                        while (V_to_prev_tr[r]->second != r)
                        {
                            r = V_to_prev_tr[r]->second;
                        }
                        r = V_to_prev_tr[r]->second;
                        Time_root& v_tr = time_root_to_tr[*V_to_prev_tr[r]];
                        std::cout << "v_tr : ";
                        Time_root::print_tr(v_tr);
                        if (V_that_have_become_children.find(v_tr.r) == V_that_have_become_children.end())
                        {
                            std::cout << "update head_tr" << std::endl;
                            std::shared_ptr<Time_root> v_tr_ptr = std::make_shared<Time_root>(v_tr);
                            head_tr.children.push_back(v_tr_ptr);
                            head_tr.num_children += v_tr.num_children;
                            head_tr.num_leaves += v_tr.num_leaves;
                            V_that_have_become_children.insert(v_tr.r);

                        }

                    }
                    else
                    {
                        std::cout << "case 2" << std::endl;
                        head_tr.birth_label.insert(v);
                        V_that_have_become_children.insert(v);
                    }
                    //if (root_to_prev_time.find(w) != root_to_prev_time.end()
                    //        && root_to_prev_time[w] < current_time)
                    if (V_to_prev_tr.find(w) != V_to_prev_tr.end())
                    {
                        std::cout << "case 3" << std::endl;
                        std::cout << "head_tr: ";
                        Time_root::print_tr(head_tr);
                        //Time_root w_tr = time_root_to_tr[std::pair<double,int>(root_to_prev_time[w],w)];
                        //Time_root& w_tr = time_root_to_tr[*V_to_prev_tr[w]];
                        int r = V_to_prev_tr[w]->second;
                        while (V_to_prev_tr[r]->second != r)
                        {
                            r = V_to_prev_tr[r]->second;
                        }
                        r = V_to_prev_tr[r]->second;
                        Time_root& w_tr = time_root_to_tr[*V_to_prev_tr[r]];

                        std::cout << "w_tr : ";
                        Time_root::print_tr(w_tr);
                        if (V_that_have_become_children.find(w_tr.r) == V_that_have_become_children.end())
                        {
                            std::cout << "update head_tr" << std::endl;
                            std::shared_ptr<Time_root> w_tr_ptr = std::make_shared<Time_root>(w_tr);
                            head_tr.children.push_back(w_tr_ptr);
                            head_tr.num_children += w_tr.num_children;
                            head_tr.num_leaves += w_tr.num_leaves;
                            V_that_have_become_children.insert(w_tr.r);
                        }
                    }
                    else
                    {
                        std::cout << "case 4" << std::endl;
                        head_tr.birth_label.insert(w);
                        V_that_have_become_children.insert(w);
                    }
                    //V_to_prev_tr[v] = head_tr;
                    //V_to_prev_tr[w] = head_tr;
                    //root_to_prev_time[head] = current_time;
                }

                std::cout << "end updating new_edges" << std::endl;
                std::cout << "upper_roots: " << std::endl;
                print_upper_roots(upper_roots);

                // figure out which of new_heads are leaves, add birth labels to singleton leaves, and update V_to_prev_tr
                std::cout << "figure out which of new_heads are leaves, add birth labels to singleton leaves, and update V_to_prev_tr" << std::endl;
                for (int v : new_heads)
                {
                    std::cout << "v = " << v << std::endl;
                    std::pair<double,int> time_root (current_time,v);
                    Time_root& tr = time_root_to_tr[time_root];
                    std::cout << "V_to_prev_tr[v] = ";
                    Time_root::print_tr(tr,true);
                    if (tr.children.size() == 0)
                    {
                        tr.num_leaves = 1;
                        tr.num_children = 1;
                        if (tr.birth_label.empty())
                            tr.birth_label.insert(tr.r);
                    }
                    std::stringstream ss;
                    for (int x : tr.birth_label)
                        ss << x << ",";
                    tr.birth_label_str = ss.str();
                    // V_to_prev_tr[v] = tr;
                    if (V_to_prev_tr.find(v) == V_to_prev_tr.end())
                    {
                        std::cout << "kase 1" << std::endl;
                        V_to_prev_tr[v] = std::make_shared<std::pair<double,int>>(time_root);
                    }
                    else
                    {
                        std::cout << "kase 2" << std::endl;
                        std::cout << "time_root.first = " << time_root.first << std::endl;
                        std::cout << "time_root.second = " << time_root.second << std::endl;
                        V_to_prev_tr[v]->first = time_root.first;
                        V_to_prev_tr[v]->second = time_root.second;
                    }
                }

                std::cout << "update V_to_prev_tr for V_that_have_become_children" << std::endl;
                for (int x : V_that_have_become_children)
                {
                    std::cout << "x = " << x << std::endl;
                    int h_x = Find(subsets.get(), x);
                    std::cout << "h_x = " << h_x << std::endl;
                    //V_to_prev_tr[x] = V_to_prev_tr[h_x];
                    if (V_to_prev_tr.find(x) == V_to_prev_tr.end())
                    {
                        std::cout << "case 1" << std::endl;
                        V_to_prev_tr[x] = V_to_prev_tr[h_x];
                    }
                    else
                    {
                        std::cout << "case 2" << std::endl;
                        V_to_prev_tr[x]->first = V_to_prev_tr[h_x]->first;
                        V_to_prev_tr[x]->second = V_to_prev_tr[h_x]->second;
                        std::cout << "V_to_prev_tr[h_x]->first = " << V_to_prev_tr[h_x]->first << std::endl;
                        std::cout << "V_to_prev_tr[h_x]->second = " << V_to_prev_tr[h_x]->second << std::endl;
                    }
                }

                std::cout << "------ time_root_to_tr ------ :" << std::endl;
                print_time_root_to_tr(time_root_to_tr,true);
                /*
                std::cout << "------ V_to_prev_tr ------ :" << std::endl;
                for (auto kv : V_to_prev_tr)
                {
                    std::cout << "#########################" << std::endl;
                    std::cout << "root = " << kv.first << " , "
                              << "prev_tr = (" << kv.second.t << "," << kv.second.r << ") , ";
                    std::stringstream ss;
                    for (int x : kv.second.birth_label)
                        ss << x << ",";
                    std::cout << "bl = " << ss.str() << std::endl;
                }*/

                // update last_upper_tr to be new_heads.begin() (when there is only one component left this is going to be the right assignment)
                /* std::cout << "update last_upper_tr" << std::endl;
                if (new_heads.size() == 1)
                {
                    int final_head = *new_heads.begin();
                    Time_root final_tr = time_root_to_tr[std::pair<double,int>(current_time,final_head)];
                    last_upper_tr = final_tr;
                }*/

                // update upper_trs
                /*
                std::cout << "update upper_trs" << std::endl;
                std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                std::cout << edgePtr->v << " "
                          << edgePtr->w << " "
                          << edgePtr->wt << std::endl;*/

                if (edgePtr->wt == std::numeric_limits<double>::infinity())
                {
                    //std::cout << "inside loop" << std::endl;
                    for (int v : upper_roots)
                    {
                        std::cout << "v = " << v << std::endl;
                        if (V_to_prev_tr.find(v) == V_to_prev_tr.end())
                            std::cout << "V_to_prev_tr[v] does not exist";
                        std::cout << "*V_to_prev_tr[v] : ";
                        std::cout << V_to_prev_tr[v]->first << " , "
                                  << V_to_prev_tr[v]->second << std::endl;

                        //upper_trs.push_back(V_to_prev_tr[v]);
                        upper_trs.push_back(time_root_to_tr[*V_to_prev_tr[v]]);
                    }
                }
                std::cout << "done updating upper_trs" << std::endl;


                if (edgePtr->wt == std::numeric_limits<double>::infinity())
                {
                    std::cout << "about to return from compute_dendrogram" << std::endl;
                    return;
                }

                std::cout << "start resetting" << std::endl;
                new_edges.clear();
                new_heads.clear();
                lower_roots.clear();
                V_that_have_become_children.clear();
                current_time = edgePtr->wt;
                //current_bigrade = edgePtr->bigrade;
                current_bigrade = appearance_to_bigrade[current_time];
                std::cout << "finished resetting" << std::endl;

                /*
                std::cout << "------ After V_to_prev_tr ------ :" << std::endl;
                for (auto kv : V_to_prev_tr)
                {
                    std::cout << "root = " << kv.first << " , "
                              << " prev_tr = (" << kv.second.t << "," << kv.second.r << ")"
                              << std::endl;
                }

                std::cout << "------ time_root_to_tr ------- : " << std::endl;
                for (auto kv : time_root_to_tr)
                {
                    std::cout << "tr = (" << kv.first.first << " , " << kv.first.second << ")"  << std::endl;
                    std::cout << "tr.children = " << std::endl;
                    //Time_root::print_vector_of_tr(kv.second.children);
                }
                std::cout << "end time_root_to_tr" << std::endl;*/
                // current_heads = new_heads;
                //continue;
            }
            if (edgePtr->v == edgePtr->w)
            {
                //std::cout << "edgePtr->v == edgePtr->w" << std::endl;
                new_heads.insert(edgePtr->v);
                upper_roots.insert(edgePtr->v);
                /*
                if (root_to_prev_time.find(edgePtr->v) == root_to_prev_time.end())
                {
                    new_heads.insert(edgePtr->v);
                    root_to_prev_time[edgePtr->v] = current_time;
                }*/
            }
            else
            {
                //std::cout << "edgePtr->v != edgePtr->w" << std::endl;
                int r_v = Find( subsets.get(), edgePtr->v );
                int r_w = Find( subsets.get(), edgePtr->w );
                //std::cout << "(r_v, r_w) = (" << r_v << "," << r_w << ")" << std::endl;

                // If including this edge does't cause cycle, include it
                // in result and increment the index of result for next edge
                if ( r_v != r_w )
                {

                    mst.push_back(edgePtr);
                    new_edges.push_back(edgePtr);

                    std::pair<int,int> roots = Union( subsets.get(), r_v, r_w );
                    int lower_root = roots.first;
                    int upper_root = roots.second;

                    // keep track of all top roots of trees in MSF
                    upper_roots.insert(upper_root); // insert will automatically ensure upper_root is unique
                    /*std::cout << "upper_roots before erase" << std::endl;
                    print_upper_roots(upper_roots);
                    if (upper_roots.find(lower_root) != upper_roots.end())
                    {
                        upper_roots.erase(lower_root);
                        std::cout << "upper_roots after erase:" << std::endl;
                        print_upper_roots(upper_roots);
                    }*/

                    lower_roots.push_back(lower_root);
                    new_heads.insert(upper_root);
                    /* POSSIBLY SWITCH LOCATIONS
                    if (V_to_prev_tr.find(r_v) != V_to_prev_tr.end())
                        V_to_prev_tr[edgePtr->v] = V_to_prev_tr[r_v];
                    if (V_to_prev_tr.find(r_w) != V_to_prev_tr.end())
                        V_to_prev_tr[edgePtr->w] = V_to_prev_tr[r_w];
                    */
                    //new_heads_to_children[upper_root].push_back(lower_root);
                    //new_heads_to_children[upper_root].push_back(upper_root);
                    //new_heads.erase(lower_root);

                    std::cout << "lower_root = " << lower_root << std::endl;
                    std::cout << "upper_root = " << upper_root << std::endl;

                    //std::cout << "MST updated: " << std::endl;
                    /*for ( std::vector<EdgePtr>::const_iterator it = mst.begin(); it != mst.end(); ++it )
                    {
                        std::cout << it->get()->v << " "
                                  << it->get()->w << " "
                                  << it->get()->wt << std::endl;
                    }*/
                }
            }
            /*std::cout << "check if last_upper_tr must be updated" << std::endl;
            if (is_last_edge)
            {
                std::cout << "1" << std::endl;
                int final_head = *new_heads.begin();
                std::cout << "2" << std::endl;
                Time_root final_tr = time_root_to_tr[std::pair<double,int>(current_time,final_head)];
                std::cout << "3" << std::endl;
                last_upper_tr = final_tr;
                std::cout << "4" << std::endl;
            }*/
        }
        std::cout << "about to normally return from compute_dendrogram" << std::endl;
        return;
    }
};

#endif // DENDROGRAM_DATA_H
