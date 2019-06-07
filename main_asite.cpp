#include <iostream>
#include <fstream>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parserargv.h"

int main_asite(int argc, const char *argv[])
{
    std::string fileBed;
    //std::vector<BamHandle> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::asite::error, provide BED file." << std::endl;
        return 1;
    }

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            //BamHandle handle;
            //handle.name = fileNameNext;
            //handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::asite::error, provide BAM file." << std::endl;
        return 1;
    }

    /*
    if (openBamHandles(handlesBam) > 0) {
        std::cerr << "ribotools::asite::error, failed to open BAM handles." << std::endl;
        return 1;
    }

    std::cout << "BED: " << fileBed << std::endl;
    std::cout << "BAMS: " << handlesBam.size() << std::endl;

    if (closeBamHandles(handlesBam) > 0) {
        std::cerr << "ribotools::asite::error, failed to close BAM handles." << std::endl;
        return 1;
    }
    */
    return 0;
}