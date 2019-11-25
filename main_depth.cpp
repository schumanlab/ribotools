#include <iostream>
#include <memory>

#include "argumentparser.h"
#include "bedio.h"
#include "bamio.h"
#include "version.h"

int main_depth(const int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<std::string> filesBam;

    auto p = ArgumentParser("depth", std::string(VERSION), "coverage per ORF from BAM files");
    p.addArgumentRequired("annotation").setKeyShort("-b").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentPositional("alignment").setCount(-1).setHelp("list of RiboSeq BAM files");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        filesBam = p.get<std::vector<std::string>>("alignment");
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }


    // open BED file
    auto hBed = BedIO(fileBed);
    if (!hBed.isOpen()) {
        std::cerr << "ribotools::" << hBed.error() << std::endl;
        return EXIT_FAILURE;
    }

    // open BAM files
    auto hBam = BamIO(filesBam);
    if (!hBam.isOpen()) {
        std::cerr << "ribotools::" << hBam.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << hBam.aux.size() << std::endl;
    for (auto handle : hBam.aux) {
        std::cout << handle->name() << "\t" << handle->count() << std::endl;
    }









    return 0;
}
