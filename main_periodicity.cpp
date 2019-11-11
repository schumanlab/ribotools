#include <iostream>

#include "argumentparser.h"
#include "bedio.h"
#include "seqio.h"
#include "version.h"

int main_periodicity(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::string fileBam;
    std::vector<std::string> filesBams;

    auto p = ArgumentParser("periodicity", std::string(VERSION), "creates P-site periodicity histogram");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentRequired("sequence").setKeyShort("-f").setKeyLong("--fasta").setHelp("FASTA file containing transcript sequence");
    p.addArgumentRequired("alignment").setKeyShort("-b").setKeyLong("--bam").setHelp("BAM file for alignment");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");


    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        fileFasta = p.get<std::string>("sequence");
        fileBam = p.get<std::string>("alignment");
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }



    auto hAnn = BedIO(fileBed);
    if (!hAnn.isOpen()) {
        std::cerr << "ribotools::" + hAnn.error() << std::endl;
        return EXIT_FAILURE;
    }

    auto hSeq = SeqIO(fileFasta);
    if (!hSeq.isOpen()) {
        std::cerr << "ribotools::" + hSeq.error() << std::endl;
        return EXIT_FAILURE;
    }




    while (hAnn.next()) {
        //std::cout << hBed.line() << std::endl;
        std::cout << hAnn.bed().name() << "\t" << hAnn.bed().span() << std::endl;
        hSeq.fetch(hAnn.bed().name(), hAnn.bed().orfStart(), hAnn.bed().orfEnd());
        std::cout << hSeq.sequence() << std::endl;
    }


    return 0;
}
