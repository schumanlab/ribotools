#include <iostream>
#include <fstream>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"

int main_count(const int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<BamHandle *> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);

    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::metagene::error, provide BED file." << std::endl;
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
        std::cerr << "ribotools::metagene::error, provide BAM file." << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::metagene::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // counting variables
    std::size_t count_files = handlesBam.size();
    int count_genes = 0;
    bam1_t *alignment = bam_init1();

    // print header
    std::cout << "# name" << "\t" << "length";
    for (auto handle : handlesBam) {
        std::cout << "\t" << handle->name();
    }
    std::cout << std::endl;

    // count total reads per file
    std::cout << "total" << "\t" << -1;
    for (auto handle : handlesBam) {
        int count_reads = 0;
        while (handle->readBam(alignment) > 0) {
            count_reads++;
        }
        std::cout << "\t" << count_reads;
    }
    std::cout << std::endl;

    // loop over BED records
    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // count reads per file
        int count_total = 0;
        std::vector<int> count_each(count_files, 0);
        int i = 0;
        for (auto handle : handlesBam) {
            handle->query(bed.transcript, bed.cdsStart, bed.cdsEnd);
            int count_reads = 0;
            while (handle->readBam(alignment) > 0) {
                count_reads++;
            }
            count_each.at(i++) = count_reads;
            count_total += count_reads;
        }

        // print result
        if (count_total > 0) {
            std::cout << bed.name << "\t" << bed.cdsSpan;
            for (auto val : count_each) {
                std::cout << "\t" << val;
            }
            std::cout << std::endl;
            count_genes++;
        }
    }

    if (alignment)
        bam_destroy1(alignment);
    std::cerr << "Used genes: " << count_genes << std::endl;

    // destructors
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}
