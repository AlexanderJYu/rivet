//
// Created by Bryn Keller on 7/6/16.
//

#ifndef RIVET_CONSOLE_SERIALIZATION_H
#define RIVET_CONSOLE_SERIALIZATION_H

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
//#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>

#include "numerics.h"
#include "anchor.h"
#include "dcel.h"
#include "math/xi_point.h"
#include "math/xi_support_matrix.h"
#include "barcode_template.h"

//namespace boost {
//    namespace multiprecision {
//        template<class Archive>
//        void load(Archive &ar, exact &num, const unsigned int version) {
//            std::string rep;
//            ar &rep;
//            num.assign(rep;
//        }
//
//        template<class Archive>
//        void save(Archive &ar, exact const &num , const unsigned int version) {
//            std::stringstream ss;
//            ss & num;
//            ar &ss.str();
//        }
//    }
//}

BOOST_CLASS_EXPORT(BarcodeTemplate);
BOOST_CLASS_EXPORT(BarTemplate);
BOOST_CLASS_EXPORT(xiMatrixEntry);
BOOST_CLASS_EXPORT(Anchor);
BOOST_CLASS_EXPORT(Face);
BOOST_CLASS_EXPORT(Halfedge);
BOOST_CLASS_EXPORT(Mesh);
BOOST_CLASS_EXPORT(Vertex);


template <class Archive>
void Vertex::serialize(Archive &ar, const unsigned int version) {
    ar &incident_edge & x & y;
}

template <class Archive>
void Halfedge::serialize(Archive &ar, const unsigned int version) {
    ar &origin & twin & next & prev & face & anchor;
}

template<class Archive>
void Face::serialize(Archive & ar, const unsigned int version) {
    ar &boundary & dbc & visited;
}

template <class Archive>
void Anchor::serialize(Archive &ar, const unsigned int version) {
    ar &x_coord & y_coord & entry & dual_line & position & above_line & weight;
}

template <class Archive>
void serialize(Archive &ar, xiMatrixEntry &x, const unsigned int version) {
    ar &x.x & x.y & x.index & x.down & x.left & x.low_simplices & x.high_simplices & x.low_count & x.high_count &
       x.low_index & x.high_index;
}

template <class Archive>
void serialize(Archive &ar, Multigrade &m, const unsigned int version) {
    ar &m.num_cols & m.simplex_index & m.x & m.y;
}

template <class Archive>
void Mesh::serialize(Archive &ar, const unsigned int version) {
ar & x_exact
    & y_exact
    & x_grades
    & y_grades;
    std::cout << "Processing faces" << std::endl;
      ar & topleft
      & topright
      & bottomleft
      & bottomright;

    std::cout << "Processing collections" << std::endl;
    ar
    & all_anchors
    & halfedges
    & faces
    & verbosity
    & vertical_line_query_list
    & vertices;
}

#include <boost/serialization/split_free.hpp>
BOOST_SERIALIZATION_SPLIT_FREE(unsigned_matrix)

namespace boost {


    template<class Archive>
    void save(Archive &ar, unsigned_matrix const &mat, const unsigned int &version) {
        assert(mat.num_dimensions() == 2);
        std::vector<unsigned> dims(mat.shape(), mat.shape() + mat.num_dimensions());
        std::vector<unsigned> data(mat.data(), mat.data() + mat.num_elements());
        ar &dims &  data;
    }

    template <class Archive>
    void load(Archive & ar, unsigned_matrix &mat, const unsigned int &version) {
        std::vector<unsigned> dims;
        std::vector<unsigned> data;
        ar &dims & data;
        unsigned_matrix::extent_gen extents;
        auto size = extents[dims[0]][dims[1]];
        mat.resize(size);
        std::memcpy(data.data(), mat.data(), data.size() * sizeof(unsigned));
    }
}

#endif //RIVET_CONSOLE_SERIALIZATION_H
