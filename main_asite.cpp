#include <iostream>
#include <fstream>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"

int main_asite(int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<BamHandle*> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::asite::error, provide BED file." << std::endl;
        return 1;
    }

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            BamHandle *handle = new BamHandle(fileNameNext);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::asite::error, provide BAM file." << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::asite::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // loop over BED records
    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;
        bed.parseExons();
        int codonsCDS = bed.cdsSpan / 3;

        // calculate coverage per file
        for (auto handle : handlesBam) {

            //std::vector<int> coverage(codonsCDS, 0);

        }

    }

    // destructors
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}