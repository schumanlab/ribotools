#include <iostream>
#include <fstream>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parserargv.hpp"

struct BamHandle {
    std::string name;
    htsFile *bam;
    hts_idx_t *bai;
    bam_hdr_t *header;
};

int openBamHandles(std::vector<BamHandle> &handlesBam);
int closeBamHandles(std::vector<BamHandle> &handlesBam);

int main_asite(int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<BamHandle> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::asite::error, provide BED file." << std::endl;
        return 1;
    }

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            BamHandle handle;
            handle.name = fileNameNext;
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::asite::error, provide BAM file." << std::endl;
        return 1;
    }

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
    return 0;
}


int openBamHandles(std::vector<BamHandle> &handlesBam)
{
    for (auto &handle: handlesBam) {

        // open BAM file
        handle.bam = hts_open(handle.name.c_str(), "r");
        if (!handle.bam) {
            std::cerr << "ribotools::asite::error, failed to open BAM file " << handle.name << std::endl;
            return 1;
        }

        // load BAI file
        handle.bai = hts_idx_load(handle.name.c_str(), HTS_FMT_BAI);
        if (!handle.bai) {
            std::cerr << "ribotools::asite::error, failed to load BAM index " << handle.name << std::endl;
            return 1;
        }

        // read BAM header
        handle.header = sam_hdr_read(handle.bam);
        if (!handle.header) {
            std::cerr << "ribotools::asite::error, failed to read BAM header " << handle.name << std::endl;
            return 1;
        }
    }
    return 0;
}


int closeBamHandles(std::vector<BamHandle> &handlesBam)
{
    for (auto &handle: handlesBam) {

        // destroy header
        if (handle.header)
            bam_hdr_destroy(handle.header);
        
        // destory index
        if (handle.bai)
            hts_idx_destroy(handle.bai);

        // close BAM file
        if (handle.bam)
            hts_close(handle.bam);
    }
    return 0;
}