#include <iostream>
#include <memory>

#include "argumentparser.h"
#include "bedio.h"
#include "bamio.h"
#include "version.h"

std::vector<int> reduceToCodons(std::vector<int> &depth, std::size_t size) {
    std::vector<int> codons(size, 0);
    std::vector<int>::iterator ov = codons.begin();
    for (std::vector<int>::const_iterator it = depth.begin();
         it != depth.end(); it += 3) {

        if (ov != codons.end()) {
            *ov = *it + *(it + 1) + *(it + 2);
            ++ov;
        }
    }
    return codons;
}



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

    std::size_t sizeColumns = hBam.aux.size();

    // loop over BED records
    while (hBed.next()) {

        // size in codons
        std::size_t sizeRows = static_cast<std::size_t>(hBed.bed().orfSpan() / 3);
        if (sizeRows < 50) continue;


        // calculate counts
        std::size_t useFlag = 0;
        for (auto handle : hBam.aux) {
                handle->query(hBed.bed().name(1), hBed.bed().orfStart(), hBed.bed().orfEnd());
                int count = handle->count();
                float avgDepth = static_cast<float>(count) / hBed.bed().orfSpan();
                if (avgDepth >= 0.1f) useFlag++;
        }

        // filter based on minimum ratio in all replica
        if (useFlag < 4) continue;

        // depth per replica
        std::vector<int> vDepth(sizeRows, 0);
        std::vector<std::vector<int>> vData(sizeColumns, vDepth);
        std::vector<std::vector<int>>::iterator iv = vData.begin();
        for (auto handle : hBam.aux) {
            std::vector<int> depthNucleotides = handle->depth(hBed.bed().name(1), hBed.bed().orfStart(), hBed.bed().orfEnd());
            std::vector<int> depthCodons = reduceToCodons(depthNucleotides, sizeRows);
            if (iv != vData.end()) {
                *iv = depthCodons;
                ++iv;
            }
        }

        // print loop
        for (std::size_t r = 0; r < sizeRows; ++r) {

            std::cout << hBed.bed().name(2) << "\t" << r;
            for (std::size_t c = 0; c < sizeColumns; ++c) {
                std::cout << "\t" << vData[c][r];
            }
            std::cout << std::endl;

        }

    }


    return 0;
}
