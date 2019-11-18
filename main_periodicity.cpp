#include <iostream>
#include <map>

#include "argumentparser.h"
#include "bamio.h"
#include "bedio.h"
#include "seqio.h"
#include "version.h"


struct Offset {
    uint32_t readLength;
    uint32_t fpSite;
    uint32_t tpSite;
};

void calculatePSiteOffset(std::map<Offset, uint32_t> &poMap, BedIO &hBed, BamIO &hBam);


int main_periodicity(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileBam;

    auto p = ArgumentParser("periodicity", std::string(VERSION), "creates P-site periodicity histogram");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentRequired("alignment").setKeyShort("-b").setKeyLong("--bam").setHelp("BAM file for alignment");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");


    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        fileBam = p.get<std::string>("alignment");
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }



    auto hBed = BedIO(fileBed);
    if (!hBed.isOpen()) {
        std::cerr << "ribotools::" + hBed.error() << std::endl;
        return EXIT_FAILURE;
    }


    auto hBam = BamIO(fileBam, 255);
    if (!hBam.isOpen()) {
        std::cerr << "ribotools::" + hBam.error() << std::endl;
        return EXIT_FAILURE;
    }


    std::map<Offset, uint32_t> map_pSiteOffset;

    calculatePSiteOffset(map_pSiteOffset, hBed, hBam);


    return 0;
}


void calculatePSiteOffset(std::map<Offset, uint32_t> &poMap, BedIO &hBed, BamIO &hBam)
{

    while (hBed.next()) {

        hBam.query(hBed.bed().name(), hBed.bed().orfStart(), hBed.bed().orfStart() + 1);
        while (hBam.next()) {
            int32_t readStart = hBam.readStart();
            int32_t readLength = hBam.readLength();
            std::cout << hBed.bed().name(2) << "\t" << readStart << "\t" << readLength << "\t" << hBed.bed().orfStart() << std::endl;
        }

    }
}
