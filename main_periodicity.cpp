#include <iostream>

#include "argumentparser.h"
#include "bedio.h"
#include "version.h"

int main_periodicity(int argc, const char *argv[])
{
    std::string fileBed;

    auto p = ArgumentParser("periodicity", std::string(VERSION), "creates P-site periodicity histogram");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");

    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    auto hAnn = BedIO(fileBed);
    while (hAnn.next()) {
        //std::cout << hBed.line() << std::endl;
        std::cout << hAnn.bed().name() << "\t" << hAnn.bed().span << std::endl;
    }


    return 0;
}
