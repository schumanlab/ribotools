#include <iostream>
#include <memory>
#include <fstream>
#include <numeric>

#include "argumentparser.h"
#include "bamhandle.h"
#include "bedrecord.h"
#include "version.h"

int main_runoff(const int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<std::string> filesBam;
    std::vector<std::shared_ptr<BamHandle>> handlesBam;
    int sizeHistogram;
    int sizeBackground;

    auto p = ArgumentParser("runoff", std::string(VERSION), "project runoff footprints normalised to transcript end");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentRequired("footprints").setKeyShort("-b").setKeyLong("--bam").setHelp("BAM files of footprint alignment").setCount(-1);
    p.addArgumentOptional("bins").setKeyShort("-n").setKeyLong("--bins").setHelp("number of histogram bins").setDefaultValue<int>(400);
    p.addArgumentOptional("background").setKeyShort("-w").setKeyLong("--window").setHelp("size of background window in 3-prime end").setDefaultValue<int>(50);
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        filesBam = p.get<std::vector<std::string>>("footprints");
        sizeHistogram = p.get<int>("bins");
        sizeBackground = p.get<int>("background");
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // open BAM handles
    for (const auto &fileName : filesBam) {
        std::shared_ptr<BamHandle> handle = std::make_shared<BamHandle>(fileName, 255, 0);
        handlesBam.emplace_back(handle);
    }


    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::pausing::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // loop over BED record
    std::string line;

    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::stringstream iss(line);
        iss >> bed;

        // filter based on length
        int sizeNucleotides = bed.span;
        int sizeCodons = bed.cdsSpan / 3;
        if ((sizeHistogram + sizeBackground + 10) >= sizeCodons) continue;

        for(auto handle : handlesBam) {

            // pileup nucleotides
            std::vector<int> fcNucleotides(static_cast<std::size_t>(sizeNucleotides), 0);
            handle->calculateFootprintCoverage(fcNucleotides, bed.transcript, 0, bed.span);

            std::vector<int>::const_iterator itCdsStart = fcNucleotides.begin() + bed.cdsStart;
            std::vector<int>::const_iterator itCdsEnd = fcNucleotides.begin() + bed.cdsEnd;

            int background = std::accumulate(fcNucleotides.begin() + bed.cdsEnd - sizeHistogram - 10,
                                             fcNucleotides.begin() + bed.cdsEnd - 10, 0);

            if (background < 50) continue;

            std::cout << bed.gene << "\t" << handle->name() << "\t" << background;
            int codonSum = 0;
            for (std::vector<int>::const_iterator it = itCdsStart; it != itCdsEnd; ++it) {
                codonSum += (*it);
                long int codonIndex = (it - itCdsStart + 1);
                if ((codonIndex % 3 == 0) && (codonIndex/3 <= sizeHistogram)) {
                    std::cout << "\t" << codonSum;
                    codonSum = 0;
                }
            }
            std::cout << std::endl;

        }

    }

    // print header
    /*
    std::cout << "# index";
    for (auto handle : handlesBam) {
        std::cout << "\t" << handle->name();
    }
    std::cout << std::endl;

    // print histogram
    for (std::vector<double>::const_iterator it = buffer.begin(); it != (buffer.begin() + sizeHistogram); ++it) {

        std::cout << (it - buffer.begin());
        for (std::size_t k = 0; k < handlesBam.size(); ++k) {
            std::vector<double>::const_iterator itNext = it + static_cast<int>(k) * sizeHistogram;
            std::cout << "\t" << (*itNext);
        }
        std::cout << std::endl;

    }
    */

    // destructors
    fhBed.close();

    return 0;
}
