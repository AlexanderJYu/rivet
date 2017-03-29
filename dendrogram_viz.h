#ifndef DENDROGRAM_VIZ_H
#define DENDROGRAM_VIZ_H

//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

#include <boost/config.hpp>
#include <iostream>                         // for std::cout
#include <utility>                          // for std::pair
#include <algorithm>                        // for std::for_each
#include <boost/utility.hpp>                // for boost::tie
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>
//#include "graph.h"
#include "time_root.h"

using namespace boost;

class Dendrogram_viz
{
public:

    template <class Graph> struct exercise_vertex
    {
      exercise_vertex(Graph& g_, const char name_[]) : g(g_),name(name_) { }
      typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
      void operator()(const Vertex& v) const
      {
        using namespace boost;
        typename property_map<Graph, vertex_index_t>::type
          vertex_id = get(vertex_index, g);
        std::cout << "vertex: " << name[get(vertex_id, v)] << std::endl;

        // Write out the outgoing edges
        std::cout << "\tout-edges: ";
        typename graph_traits<Graph>::out_edge_iterator out_i, out_end;
        typename graph_traits<Graph>::edge_descriptor e;
        for (boost::tie(out_i, out_end) = out_edges(v, g);
             out_i != out_end; ++out_i)
        {
          e = *out_i;
          Vertex src = source(e, g), targ = target(e, g);
          std::cout << "(" << name[get(vertex_id, src)]
                    << "," << name[get(vertex_id, targ)] << ") ";
        }
        std::cout << std::endl;

        // Write out the incoming edges
        std::cout << "\tin-edges: ";
        typename graph_traits<Graph>::in_edge_iterator in_i, in_end;
        for (boost::tie(in_i, in_end) = in_edges(v, g); in_i != in_end; ++in_i)
        {
          e = *in_i;
          Vertex src = source(e, g), targ = target(e, g);
          std::cout << "(" << name[get(vertex_id, src)]
                    << "," << name[get(vertex_id, targ)] << ") ";
        }
        std::cout << std::endl;

        // Write out all adjacent vertices
        std::cout << "\tadjacent vertices: ";
        typename graph_traits<Graph>::adjacency_iterator ai, ai_end;
        for (boost::tie(ai,ai_end) = adjacent_vertices(v, g);  ai != ai_end; ++ai)
          std::cout << name[get(vertex_id, *ai)] <<  " ";
        std::cout << std::endl;
      }
      Graph& g;
      const char *name;
    };

    template <class RootMap,class BigradeMap,class BirthlabelMap>
    class vertex_writer {
    public:
      vertex_writer(RootMap r, BigradeMap b, BirthlabelMap bl) : rm(r),bm(b), blm(bl) {}
      template <class Vertex>
      void operator()(std::ostream &out, const Vertex& v) const {
        out << "[label=\""
            << rm[v] << "\\n"
            << bm[v] << "\\n"
            << blm[v] << "\"]";
      }
    private:
      RootMap rm;
      BigradeMap bm;
      BirthlabelMap blm;
    };

    template <class RootMap, class BigradeMap, class BirthlabelMap>
    inline vertex_writer<RootMap,BigradeMap,BirthlabelMap>
    static make_vertex_writer(RootMap r,BigradeMap b,BirthlabelMap bl)
    {
        return vertex_writer<RootMap,BigradeMap,BirthlabelMap>(r,b,bl);
    };

    struct vert
    {
       std::string root;
       std::string bigrade;
       std::string birthlabel;
    };

    //typedef std::map<std::string, std::string> GraphvizAttributes;

    typedef adjacency_list<vecS, vecS, bidirectionalS, vert> Graph;
    /*typedef subgraph <
      adjacency_list< vecS, vecS, directedS, vert,
      property<edge_index_t, int, property<edge_attribute_t, GraphvizAttributes> >,
      property<graph_name_t, std::string,
            property<graph_graph_attribute_t,  GraphvizAttributes,
            property<graph_vertex_attribute_t, GraphvizAttributes,
            property<graph_edge_attribute_t,   GraphvizAttributes>>>>
        >
    > Graph;*/
    /*typedef subgraph< adjacency_list< vecS, vecS, bidirectionalS,
      vert, property< edge_index_t, int > > > Graph;*/
    typedef boost::graph_traits<Graph>::vertex_descriptor vert_id;

    static void add_top_vertex(Graph& g,
                        Time_root& tr,
                        const std::vector<double>& x_grades,
                        const std::vector<double>& y_grades)
    {
        vert prop;
        std::stringstream ss_r;
        ss_r << tr.r;
        prop.root = ss_r.str();
        std::stringstream ss_b;
        ss_b << "( " << x_grades[tr.bigrade.first] << " , " << y_grades[tr.bigrade.second] << ")";
        prop.bigrade = ss_b.str();
        std::stringstream ss_bl;
        for (auto u : tr.birth_label)
            ss_bl << u << ",";
        prop.birthlabel = ss_bl.str();
        vert_id v = add_vertex(prop,g);
        for (std::shared_ptr<Time_root> child_tr : tr.children)
        {
            recursively_define_g(g,*child_tr,v,x_grades,y_grades);
        }
    }

    static void write_dendrogram_dot_file(Time_root& last_upper_tr,
                                          std::string file_name,
                                          const std::vector<double>& x_grades,
                                          const std::vector<double>& y_grades)
    {
        // create a typedef for the Graph type
        Graph g;
        //Graph& g_0 = g.create_subgraph();
        //typedef std::pair<int,int> Edge;
        //Edge edge_array[];
        //int v_id = 0;
        //std::unordered_map<vert_id, vert> id_to_vert;
        //std::string name[V];
        enum {is_forest_connector = -2};
        if (last_upper_tr.r != is_forest_connector)
        {
            add_top_vertex(g,last_upper_tr,x_grades,y_grades);
        }
        else
        {
            for (std::shared_ptr<Time_root> child_tr : last_upper_tr.children)
            {
                add_top_vertex(g,*child_tr,x_grades,y_grades);
            }
        }

        /*
        std::stringstream ss_label;
        ss_label << ss_r << "\\n" << ss_b << "\\n" << ss_bl;
        //auto v = add_vertex(prop,g_0);
        vert_id v = add_vertex(g_0);
        put(v, g_0, v, prop);*/
        //name[v] = ss_label.str();


        std::map<std::string,std::string> graph_attr, vertex_attr, edge_attr;
        graph_attr["size"] = "3,3";
        graph_attr["rankdir"] = "LR";
        graph_attr["ratio"] = "fill";
        vertex_attr["shape"] = "circle";

        std::filebuf fb;
        fb.open (file_name,std::ios::out);
        std::ostream os(&fb);

        boost::write_graphviz(os, g,
                              make_vertex_writer(boost::get(&vert::root,g),
                                                 boost::get(&vert::bigrade,g),
                                                 boost::get(&vert::birthlabel,g)),
                              default_writer(),
                              make_graph_attributes_writer(graph_attr, vertex_attr,
                                                           edge_attr));
        /*boost::write_graphviz(os, g,
                              make_label_writer(name),
                              default_writer(),
                              make_graph_attributes_writer(graph_attr, vertex_attr,
                                                           edge_attr));*/
    }

    static void recursively_define_g( Graph& g,
                                      Time_root& current_tr,
                                      vert_id& parent,
                                      const std::vector<double>& x_grades,
                                      const std::vector<double>& y_grades)
    {
        vert prop;
        std::stringstream ss_r;
        ss_r << current_tr.r;
        prop.root = ss_r.str();
        std::stringstream ss_b;
        ss_b << "( " << x_grades[current_tr.bigrade.first] << " , " << y_grades[current_tr.bigrade.second] << ")";
        prop.bigrade = ss_b.str();
        std::stringstream ss_bl;
        for (auto u : current_tr.birth_label)
            ss_bl << u << ",";
        prop.birthlabel = ss_bl.str();
        //vert_id v = add_vertex(prop,g);
        /*
        std::stringstream ss_label;
        ss_label << ss_r << "\\n" << ss_b << "\\n" << ss_bl;
        //auto v = add_vertex(prop,g_0);
        vert_id v = add_vertex(g);
        put(prop, g, v);
        //name[v] = ss_label.str();
        add_edge(v, parent, g);*/


        auto v = add_vertex(prop, g);
        add_edge(v, parent, g);

        if (current_tr.children.size() == 0)
            return;

        for (std::shared_ptr<Time_root> child_tr : current_tr.children)
        {
            recursively_define_g(g,*child_tr,v,x_grades,y_grades);
        }
    }
};

/*
int main(int,char*[])
{

  // create a typedef for the Graph type
  typedef adjacency_list<vecS, vecS, bidirectionalS, vert> Graph;

  // Make convenient labels for the vertices
  enum { A, B, C, D, E, F, G, N };
  const int num_vertices = N;
  const char name[] = "ABCDEFG";

  // writing out the edges in the graph
  typedef std::pair<int,int> Edge;
  Edge edge_array[] =
  { Edge(A,B), Edge(A,C), Edge(B,D), Edge(B,E), Edge(C,F), Edge(C,G)};
  const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);

  Graph g;
  //property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
  for (std::size_t j = 0; j < num_vertices; ++j)
  {
    vert prop;
    prop.root = name[j];
    std::stringstream ss;
    ss << "( " << j << " , " << j + 1 << ")";
    prop.bigrade = ss.str();
    add_vertex(prop,g);
  }

  for (std::size_t j = 0; j < num_edges; ++j)
  {
    add_edge(edge_array[j].first,edge_array[j].second, g);
  }

  /*
  boost::property_map<Graph, vertex_index_t>::type
    vertex_id = get(vertex_index, g);

  std::cout << "vertices(g) = ";
  typedef graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> vp;
  for (vp = vertices(g); vp.first != vp.second; ++vp.first)
    std::cout << name[get(vertex_id, *vp.first)] <<  " ";
  std::cout << std::endl;

  std::cout << "edges(g) = ";
  graph_traits<Graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei)
    std::cout << "(" << name[get(vertex_id, source(*ei, g))]
              << "," << name[get(vertex_id, target(*ei, g))] << ") ";
  std::cout << std::endl;

  std::for_each(vertices(g).first, vertices(g).second,
                exercise_vertex<Graph>(g, name));

  std::map<std::string,std::string> graph_attr, vertex_attr, edge_attr;
  graph_attr["size"] = "3,3";
  graph_attr["rankdir"] = "LR";
  graph_attr["ratio"] = "fill";
  vertex_attr["shape"] = "circle";

  std::filebuf fb;
  fb.open ("1.dot",std::ios::out);
  std::ostream os(&fb);

  boost::write_graphviz(os, g,
                        make_vertex_writer(boost::get(&vert::root,g),boost::get(&vert::bigrade,g)),
                        default_writer(),
                        make_graph_attributes_writer(graph_attr, vertex_attr,
                                                     edge_attr));



  return 0;
}*/


#endif // DENDROGRAM_VIZ_H
