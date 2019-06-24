#include <iostream>
#include <fstream>
#include <map>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"

int main_poffset(int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<BamHandle *> handlesBam;
    

    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::poffset::error, provide a BED file." << std::endl;
        return 1;
    }

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            auto handle = new BamHandle(fileNameNext, 255, 0);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::poffset::error, provide BAM file." << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::poffset::error, failed to open BED reference " << fileBed << std::endl;
        return 1;
    }

    // read bed file
    int counter = 0;
    std::string line;
    std::map<std::string, std::map<int, std::map<int,int>>> offset_map;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;
        bed.parseExons();

        for (auto handle : handlesBam) {
            handle->query(bed.transcript, bed.cdsStart - 1, bed.cdsStart);
            int reads = 0;
            
            bam1_t *b = bam_init1();
            while (handle->readBam(b) > 0) {
                int readStart = b->core.pos;
                int readLength = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
                int readEnd = readStart + readLength;
                int offsetStart = readStart - bed.cdsStart;
                int offsetEnd = readEnd - bed.cdsStart;
                offset_map[handle->name()][readLength][offsetStart]++;
                offset_map[handle->name()][readLength][offsetEnd]++;
                reads++;
            }
            bam_destroy1(b);
        }

        counter++;

    }

    for (auto name : offset_map)
    {
        for (auto length : name.second) {
            int reads = 0;
            int best_min = 0;
            int best_max = 0;
            int value_min = 0;
            int value_max = 0;
            for (auto offset : length.second) {
                if (offset.second > value_min && offset.first < 0) {
                    value_min = offset.second;
                    best_min = offset.first;
                }

                if (offset.second > value_max && offset.first > 0) {
                    value_max = offset.second;
                    best_max = offset.first;
                }

                reads += offset.second;
            }

            std::cout << name.first << "\t" << length.first << "\t" << reads << "\t" << best_min << "\t" << best_max << "\t" << value_min << "\t" << value_max << std::endl;
        }


    }

    std::cerr << "BedRecords: " << counter << std::endl;


    // destructors
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}