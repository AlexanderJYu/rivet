//
// Created by Bryn Keller on 6/21/16.
//

#ifndef RIVET_CONSOLE_NUMERICS_H
#define RIVET_CONSOLE_NUMERICS_H

#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/multiprecision/cpp_int.hpp>

typedef boost::multiprecision::cpp_rational exact;

namespace rivet
{
    namespace numeric {
        exact str_to_exact(const std::string& str);
        bool is_number(const std::string& str);
    }

}


#endif //RIVET_CONSOLE_NUMERICS_H
