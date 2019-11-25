#include <iostream>
#include <map>

#include "argumentparser.h"
#include "bamio.h"
#include "bedio.h"
#include "seqio.h"
#include "version.h"


struct Offset {
    int32_t span, fp, tp;

    bool operator==(const Offset &qry) const {
        return (span == qry.span) && (fp == qry.fp) && (tp == qry.tp);
    }

    bool operator<(const Offset &qry) const {
        return (span < qry.span) || (fp < qry.fp) || (tp < qry.tp);
    }
};

void calculatePSiteOffset(std::map<Offset, uint32_t> &posMap, BedIO &hBed, BamIO &hBam){

    while (hBed.next()) {

        hBam.query(hBed.bed().name(1), hBed.bed().orfStart(), hBed.bed().orfStart() + 1);
        while (hBam.next()) {
            int32_t readStart = hBam.readStart();
            int32_t readEnd = hBam.readEnd();
            int32_t readLength = hBam.readLength();

            Offset key = {readLength, hBed.bed().orfStart() - readStart, readEnd - hBed.bed().orfStart()};

            posMap[key]++;
            //std::cout << hBed.bed().name(2) << "\t" << readStart << "\t" << readEnd << "\t" << hBed.bed().orfStart() << std::endl;
        }

    }
}


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

    for (auto const &[key, value] : map_pSiteOffset) {
        std::cout << key.span << "::" << key.fp << "-" << key.tp << "\t" << value << std::endl;
    }


    return 0;
}



