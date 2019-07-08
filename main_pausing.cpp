#include <iostream>
#include <fstream>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"

int main_pausing(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::vector <BamHandle*> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);

    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::pausing::error, provide BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::pausing::error, provide FASTA file." << std::endl;
        return 1;
    }

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            BamHandle *handle = new BamHandle(fileNameNext, 255, 0);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::pausing::error, provide BAM file." << std::endl;
        return 1;
    }


    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::pausing::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::pausing::error, failed to load fasta reference " << fileFasta << std::endl;
        return 1;
    }

    // loop over BED record
    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;
    }

    // destructors
    fhBed.close();
    fai_destroy(fhFai);
    for (auto handle : handlesBam)
        delete handle;


    return 0;
}