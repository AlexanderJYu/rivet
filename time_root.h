#ifndef TIME_ROOT_H
#define TIME_ROOT_H

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <string>
#include <set>

namespace std
{
    template <> struct hash<std::pair<int, int>> {
        inline size_t operator()(const std::pair<int, int> &v) const {
            std::hash<int> int_hasher;
            return int_hasher(v.first) ^ int_hasher(v.second);
        }
    };
}
/*
namespace std
{
    template <> struct hash<Time_root> {
        inline size_t operator()(const Time_root &tr) const {
            std::hash<int> int_hasher;
            return int_hasher(tr.t) ^ int_hasher(tr.r);
        }
    };
}*/


class Time_root
{
public:
    double t; // 1-dim filtration parameter
    int r; //0-simplex label (representative from UF)
    int num_children; // I believe num_leaves replaced the function of num_children (used for drawing dendrograms in a planar fashion)
    int num_leaves; //number of descendants that are leaves
    std::vector<std::shared_ptr<Time_root>> children;
    //std::vector<int> comp_vertices;
    std::unordered_set<int> birth_label;
    std::string birth_label_str; // string version of birth_label
    std::pair<unsigned,unsigned> bigrade; //2-dim filtration paramter (possibly should combine t and bigrade into one parameter)
    std::unordered_set<int> component_label;
    std::set<int> ordered_component_label; // sorted version of component_label (used for comparing templates to naive algorithm)
    Time_root(double t = -1, int r = -1, int num_children = 0, int num_leaves = 0) :
        t(t), r(r), num_children(num_children), num_leaves(num_leaves) { }

    template <class Archive>
    void serialize(Archive& ar, const unsigned int /*version*/)
    {
        ar& t& r& num_leaves& children& birth_label_str& bigrade;
    }

    // various functions for printing data structures and testing

    static void print_tr(Time_root& tr, bool with_comp_label = false)
    {
        std::stringstream ss_bl;
        ss_bl << "(";
        for (auto v : tr.birth_label)
            ss_bl << v << ",";
        ss_bl << ")";
        std::cout << "r = " << tr.r << " , ";
        std::cout << "b = (" << tr.bigrade.first << "," << tr.bigrade.second << ") , ";
        std::cout << "bl = " << ss_bl.str() << std::endl;
        if (with_comp_label)
        {
            std::stringstream ss_cl;
            ss_cl << "(";
            for (auto v : tr.component_label)
                ss_cl << v << ",";
            ss_cl << ")";
            std::cout << "cl = " << ss_cl.str() << std::endl;
        }
    }

    /*Time_root(Time_root& tr) :
        t(tr.t), r(tr.r), num_children(tr.num_children), num_leaves(tr.num_leaves),
        children(tr.children), birth_label(tr.birth_label), bigrade(tr.bigrade),
        component_label(tr.component_label) { } */

    static void print_tr2(const Time_root& tr, bool use_t)
    {
        std::cout << "r = " << tr.r << " , ";
        if (use_t)
            std::cout << "t = " << tr.t << " , ";
        std::cout << "b = (" << tr.bigrade.first << "," << tr.bigrade.second << ") , ";
        std::cout << "bl = " << tr.birth_label_str << " , ";
        std::cout << "cl = " << int_set_to_string(tr.ordered_component_label) << ",";
        std::cout << "num_leaves = " << tr.num_leaves << std::endl;
    }

    static void recursively_print_tr(const Time_root& tr, bool use_t = false)
    {
        std::cout << "=========================" << std::endl;
        print_tr2(tr, use_t);
        std::cout << "children: " << std::endl;
        for (auto child : tr.children)
            print_tr2(*child,use_t);
        for (auto child : tr.children)
        {
            recursively_print_tr(*child, use_t);
        }
        return;
    }

    static std::string int_set_to_string(std::set<int> comp_label)
    {
        std::stringstream ss;
        for (int v : comp_label)
        {
            ss << v << ",";
        }
        return ss.str();
    }

    static void print_trs(std::vector<Time_root>& trs, bool with_comp_label = false)
    {
        for (Time_root& tr: trs)
        {
            if (with_comp_label)
                print_tr(tr, true);
            else
                print_tr(tr);
        }
    }

    static void print_children(std::vector<std::shared_ptr<Time_root>>& trs, bool with_comp_label = false)
    {
        for (auto tr_ptr: trs)
        {
            if (with_comp_label)
                print_tr(*tr_ptr, true);
            else
                print_tr(*tr_ptr);
        }
    }

    /*
    static void print_vector_of_tr(std::vector<Time_root> children)
    {
        for(Time_root child : children)
        {
            std::cout << "child.t = " << child.t << " , "
                      << "child.r = " << child.r << " , "
                      << "child.num_children = " << child.num_children << std::endl;
        }
    }*/
    /*
    static void print_vector_of_comp(std::vector<int> comp_vertices)
    {
        std::stringstream ss;
        ss.str("");
        for (int v : comp_vertices)
        {
            ss << v << ",";
        }
        std::cout << "comp_vertices = " << ss.str() << std::endl;
    }*/
    /*
    static void print_map(std::unordered_map<std::pair<int,int>, Time_root> time_root_to_tr)
    {
        for(auto kv : time_root_to_tr)
        {
            std::cout << "time = " << kv.first.first << " , "
                      << "root = " << kv.first.second << " , "
                      << "num_children = " << kv.second.num_children << " , "
                      << "children = " << std::endl;
            print_vector_of_tr(kv.second.children);
            print_vector_of_comp(kv.second.comp_vertices);
            std::cout << "###################################" << std::endl;
        }
    }*/
};

#endif // TIME_ROOT_H
