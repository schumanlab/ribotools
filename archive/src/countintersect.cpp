#include <iostream>
#include <string>
#include <sstream>

#include "version.hpp"
#include "bedrecord.hpp"

int countintersect(int arg, char const *argv[])
{
    for (std::string line; std::getline(std::cin, line);) {
        auto ss = std::stringstream(line);
        auto ref = BedRecord();
        auto qry = BedRecord();
        int dist = 0;
        ss >> ref;
        ss >> qry;
        ss >> dist;
        std::cout << line << std::endl;
        std::cout << ref.name << " " << qry.name << " " << dist << std::endl;
        break;
    }
    return 0;
}