#include "dendrogram_data.h"
#include "subset.h"
#include "math/simplex_tree.h"
#include "math/bifiltration_data.h"
#include <numerics.h>

#ifndef COMPUTING_S_H
#define COMPUTING_S_H

namespace std
{
    template <> struct hash<std::pair<unsigned, unsigned>> {
        inline size_t operator()(const std::pair<unsigned, unsigned> &v) const {
            std::hash<unsigned> unsigned_hasher;
            return unsigned_hasher(v.first) ^ unsigned_hasher(v.second);
        }
    };
}

// class used to compute S = T_0 \cup T_1 (section 4 of paper)
class computing_s
{
private:
public:
    //Dendrogram_data::Subset oracle[][][]; // used later for computing dendrogram templates
    boost::unordered::unordered_map<std::pair<int,int>,std::map<int, Subset_map>> oracle;

    //Subset UF[]; // is reset at each new column; UF[x_0][v] is the instance of v in column x = x_0
    // 1-critical version SimplexTree* bifiltration;
    BifiltrationData* bifiltration;
    std::unordered_set<std::pair<int,int>> T_0; // every bigrade where a 0-simplex is born
    std::unordered_set<std::pair<int,int>> T_1; // Supp(xi_1) I think
    int m_x; //largest x-coordinate of bifiltration
    int m_y; //largest_y-coordinate of bifiltration
    std::vector<exact> x_grades;
    std::vector<exact> y_grades;
    //int V; //number of 0-simplices in bifiltration
    boost::unordered::unordered_map<int, std::unordered_set<unsigned>>  MSF_vertices; //called Q in paper
    boost::unordered::unordered_map<int, std::vector<std::pair<unsigned,unsigned>>> MSF_edges; //called Q in paper
    boost::unordered::unordered_map<std::pair<int,int>, std::unordered_set<unsigned>> bigrade_to_vertices; // map from a bigrade to the 0-simplices born there
    boost::unordered::unordered_map<std::pair<int,int>, std::unordered_set<std::pair<unsigned,unsigned>>> bigrade_to_edges; // map from a bigrade to the 1-simplices born there
    std::unordered_set<int> x_coordinates; // set of all x such that there exists a bigrade of the form (x,_) in the bifiltration
    std::unordered_set<int> y_coordinates; // set of all y such that there exists a bigrade of the form (_,y) in the bifiltration

    computing_s(/*SimplexTree* bif*/BifiltrationData* bif, int max_x, int max_y, std::vector<exact> x_exact, std::vector<exact> y_exact) :
        bifiltration(bif),
        m_x(max_x),
        m_y(max_y),
        x_grades(x_exact),
        y_grades(y_exact)
    {
        /* 1-critical version
        SimplexSet vertices = bifiltration->get_ordered_simplices();
        SimplexSet edges = bifiltration->get_ordered_high_simplices();
        V = vertices.size();*/

        SimplexInfo* vertices_to_grades = bifiltration->getSimplices(0);
        SimplexInfo* edges_to_grades = bifiltration ->getSimplices(1);

        for (auto vertexVec_grade : *vertices_to_grades)
        {
            unsigned v = vertexVec_grade.first[0];
            AppearanceGrades grades = vertexVec_grade.second;
            for (Grade g : grades)
            {
                bigrade_to_vertices[std::pair<int,int>(g.x,g.y)].insert(v);
                T_0.insert(std::pair<int,int>(g.x,g.y));
                x_coordinates.insert(g.x);
                y_coordinates.insert(g.y);
            }
        }

        for (auto edgeVec_grade : *edges_to_grades)
        {
            std::vector<int> edge_vector = edgeVec_grade.first;
            AppearanceGrades grades = edgeVec_grade.second;
            std::pair<unsigned,unsigned> edge (edge_vector[0], edge_vector[1]);
            for (Grade g : grades)
            {
                bigrade_to_edges[std::pair<int,int>(g.x,g.y)].insert(edge);
                x_coordinates.insert(g.x);
                y_coordinates.insert(g.y);
            }
        }

        /* 1-critical version
        for (auto st_node : vertices)
        {
            unsigned v = st_node->get_vertex();
            int x = st_node->grade_x();
            int y = st_node->grade_y();
            bigrade_to_vertices[std::pair<int,int>(x,y)].insert(v);
            T_0.insert(std::pair<int,int>(x,y));
            x_coordinates.insert(x);
            y_coordinates.insert(y);
        }

        for (auto st_node : edges)
        {
            int global_index = st_node->global_index();
            std::vector<int> edge_vector = bifiltration->find_vertices(global_index);
            std::pair<unsigned,unsigned> edge (edge_vector[0], edge_vector[1]);
            int x = st_node->grade_x();
            int y = st_node->grade_y();
            bigrade_to_edges[std::pair<int,int>(x,y)].insert(edge);
            x_coordinates.insert(x);
            y_coordinates.insert(y);
        }*/
    }

    std::unordered_set<std::pair<int,int>> get_T_0 () { return T_0; }
    std::unordered_set<std::pair<int,int>> get_T_1 () { return T_1; }
    boost::unordered::unordered_map<std::pair<int,int>,std::map<int, Subset_map>> get_oracle ()
    {
        return oracle;
    }

    //stores the S support points in lexicographical order
    void store_S_points(std::vector<TemplatePoint>& tpts)
    {
        for(int i = 0; i < m_x; i++)
        {
            for(int j = 0; j < m_y; j++)
            {
                std::pair<int,int> current (i,j);
                if (T_0.find(current) != T_0.end() || T_1.find(current) != T_1.end())
                    tpts.push_back( TemplatePoint(i,j,0,0,0)); //should I make the latter 3 arguments 0's?
            }
        }
    }//end store_support_points()

    // various functions for printing data structures and testing
    void print_T_0()
    {
        std::cout << "T_0:" << std::endl;
        for (auto grade : T_0)
        {
            std::cout << "(" << grade.first << "," << grade.second << ")" << " === "
                      << "(" << x_grades[grade.first] << "," << y_grades[grade.second] << ")" << std::endl;

        }
    }

    void print_T_1()
    {
        std::cout << "T_1:" << std::endl;
        for (auto grade : T_1)
        {
            std::cout << "(" << grade.first << "," << grade.second << ")" << " === "
                      << "(" << x_grades[grade.first] << "," << y_grades[grade.second] << ")" << std::endl;
        }
    }

    void print_oracle()
    {
        std::cout << "oracle:" << std::endl;
        for (auto grade_UF : oracle)
        {
            std::pair<int,int> grade = grade_UF.first;
            std::map<int, Subset_map> UF = grade_UF.second;
            std::cout << "-----------------" << std::endl;
            std::cout << "grade = (" << grade.first << "," << grade.second << ")" << std::endl;
            for (auto v_subset : UF)
            {
                std::cout << "v = " << v_subset.first << std::endl;
                std::cout << "v.parent = " << v_subset.second.parent << std::endl;
                std::cout << "v.rank = " << v_subset.second.rank << std::endl;
            }
        }
    }

    // main function to compute S
    void compute_s()
    {
        // handle column x = 0
        // Allocate memory for creating V ssubsets
        std::map<int, Subset_map> UF;
        std::unordered_set<unsigned> unique_vertices; // keep track of vertices currently in the UF structure

        for (int y = 0; y <= m_y; y++)
        {
            std::cout << "y = " << y << std::endl;
            std::pair<int,int> g = std::pair<int,int>(0,y);
            for (unsigned v : bigrade_to_vertices[g])
            {
                std::cout << "v = " << v << std::endl;
                if (unique_vertices.find(v) == unique_vertices.end())
                {
                    std::cout << "v is not in unique_vertices" << std::endl;
                    UF[v].parent = v;
                    UF[v].rank = 0;
                    unique_vertices.insert(v);
                    MSF_vertices[y].insert(v);
                }
            }

            bool merged = false;
            for (std::pair<unsigned,unsigned> e : bigrade_to_edges[g])
            {
                std::cout << " ======== " << std::endl;
                std::cout << "e = (" << e.first << "," << e.second << ")" << std::endl;
                int a = Find( UF, e.first );
                int b = Find( UF, e.second );
                std::cout << "(a,b) = (" << a << "," << b << ")" << std::endl;

                if ( a != b )
                {
                    std::cout << "a != b" << std::endl;
                    Union( UF, a, b );
                    std::cout << "Find(UF, a) = " << Find(UF, a) << std::endl;
                    std::cout << "Find(UF, b) = " << Find(UF, b) << std::endl;
                    //MSF_edges[y].insert(e);
                    MSF_edges[y].push_back(e);
                    /*if (bigrade_to_vertices[g].find(e.first) == bigrade_to_vertices[g].end() &&
                        bigrade_to_vertices[g].find(e.second) == bigrade_to_vertices[g].end())*/
                    if (bigrade_to_vertices[g].find(a) == bigrade_to_vertices[g].end() &&
                        bigrade_to_vertices[g].find(b) == bigrade_to_vertices[g].end())
                    {
                        std::cout << "merged" << std::endl;
                        merged = true;
                    }
                }
            }
            if (merged)
                T_1.insert(g);

            //update oracle
            if (x_coordinates.find(g.first) != x_coordinates.end() &&
                y_coordinates.find(g.second) != y_coordinates.end())
            {
                oracle[g] = UF;
            }
        }
        //std::cout << "got past first column" << std::endl;
        // handle columns 1 <= x <= m_x
        for (int x = 1; x <= m_x; x++)
        {
            for (int y = 0; y <= m_y; y++)
            {
                // step 1
                std::cout << "(x,y) = (" << x << "," << y << ")" << std::endl;
                std::pair<int,int> g = std::pair<int,int>(x,y);
                if (y == 0)
                {
                    unique_vertices.clear();
                    UF.clear();
                }

                // step 2
                std::cout << "====== STEP 2 ======" << std::endl;
                bool vertex_identified = false; // one of the conditions for g being in T_1
                std::cout << "MSF_vertices[y].size() = " << MSF_vertices[y].size() << std::endl;
                std::vector<unsigned> verts_to_erase;
                for (unsigned v : MSF_vertices[y])
                {
                    std::cout << "v = " << v << std::endl;
                    if (unique_vertices.find(v) == unique_vertices.end())
                    {
                        std::cout << "v is not in unique_vertices" << std::endl;
                        UF[v].parent = v;
                        UF[v].rank = 0;
                        unique_vertices.insert(v);
                    }
                    else
                    {
                        std::cout << "v is in unique_vertices" << std::endl;
                        //MSF_vertices[y].erase(v);
                        verts_to_erase.push_back(v);
                        vertex_identified = true;
                    }
                }
                for (unsigned v : verts_to_erase)
                {
                    MSF_vertices[y].erase(v);
                }
                if (vertex_identified)
                {
                    std::cout << "vertex identified" << std::endl;
                    T_1.insert(g);
                }
                //for (std::pair<unsigned,unsigned> e : MSF_edges[y])
                for (std::vector<std::pair<unsigned,unsigned>>::iterator it = MSF_edges[y].begin() ; it != MSF_edges[y].end(); ++it)
                {
                    std::pair<unsigned,unsigned> e = *it;
                    std::cout << " ======== " << std::endl;
                    std::cout << "e = (" << e.first << "," << e.second << ")" << std::endl;
                    int a = Find( UF, e.first );
                    int b = Find( UF, e.second );
                    std::cout << "(a,b) = (" << a << "," << b << ")" << std::endl;

                    if ( a != b )
                    {
                        std::cout << "a != b" << std::endl;
                        Union( UF, a, b );
                    }
                }

                // step 3 and 4
                std::cout << "====== STEP 3 + 4 ======" << std::endl;
                //bool has_relevant_simplices = false;
                bool merged = false;
                for (unsigned v : bigrade_to_vertices[g])
                {
                    if (unique_vertices.find(v) == unique_vertices.end())
                    {
                        UF[v].parent = v;
                        UF[v].rank = 0;
                        unique_vertices.insert(v);
                        MSF_vertices[y].insert(v);
                        //has_relevant_simplices = true;
                        //std::cout << "has_relevant_simplices = true" << std::endl;
                    }
                }
                for (std::pair<unsigned,unsigned> e : bigrade_to_edges[g])
                {
                    std::cout << " ======== " << std::endl;
                    std::cout << "e = (" << e.first << "," << e.second << ")" << std::endl;
                    int a = Find( UF, e.first );
                    int b = Find( UF, e.second );
                    std::cout << "(a,b) = (" << a << "," << b << ")" << std::endl;

                    if ( a != b )
                    {
                        std::cout << "a != b" << std::endl;
                        Union( UF, a, b );
                        //MSF_edges[y].insert(e);
                        MSF_edges[y].push_back(e);
                        /*if (bigrade_to_vertices[g].find(e.first) == bigrade_to_vertices[g].end() &&
                            bigrade_to_vertices[g].find(e.second) == bigrade_to_vertices[g].end())*/
                        if (bigrade_to_vertices[g].find(a) == bigrade_to_vertices[g].end() &&
                            bigrade_to_vertices[g].find(b) == bigrade_to_vertices[g].end())
                        {
                            std::cout << "merged" << std::endl;
                            merged = true;
                        }
                    }
                }
                if (merged)
                {
                    std::cout << "T_1.insert(g)" << std::endl;
                    T_1.insert(g);
                }

                //update oracle
                if (x_coordinates.find(g.first) != x_coordinates.end() &&
                    y_coordinates.find(g.second) != y_coordinates.end())
                {
                    oracle[g] = UF;
                }
            }
        }
    }
};

#endif // COMPUTING_S_H













