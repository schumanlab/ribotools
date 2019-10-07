#include <iostream>
#include <fstream>
#include <unordered_map>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"

struct PSite {
    std::string name;
    int readLength;
    int offset5p;
    int offset3p;

    bool operator==(const PSite &other) const {
        return (name == other.name &&
                readLength == other.readLength &&
                offset5p == other.offset5p &&
                offset3p == other.offset3p);
    }
};

struct PSiteHasher {
    std::size_t operator()(const PSite &other) const {
        size_t value = 17;
        value = value * 31 + std::hash<std::string>()(other.name);
        value = value * 31 + std::hash<int>()(other.readLength);
        value = value * 31 + std::hash<int>()(other.offset5p);
        value = value * 31 + std::hash<int>()(other.offset3p);
        return value;
    }
};



int main_poffset(const int argc, const char *argv[])
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
    std::unordered_map<PSite, int, PSiteHasher> psite_map;
    
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        for (auto handle : handlesBam) {
            std::string name = handle->name();
            handle->query(bed.transcript, bed.cdsStart - 1, bed.cdsStart);
            int reads = 0;
            
            bam1_t *b = bam_init1();
            while (handle->readBam(b) > 0) {
                int readStart = b->core.pos;
                int readLength = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
                int readEnd = readStart + readLength;
                int offsetStart = readStart - bed.cdsStart;
                int offsetEnd = readEnd - bed.cdsStart;
                
                PSite next = {name, readLength, offsetStart, offsetEnd};
                auto node = psite_map.find(next);
                if (node == psite_map.end())
                    psite_map[next] = 1;

                psite_map[next]++;


                reads++;
            }
            bam_destroy1(b);
        }

        counter++;

    }

    for (auto node : psite_map) {
        std::cout << node.first.name << "\t" 
                  << node.first.readLength << "\t"
                  << node.first.offset5p << "\t"
                  << node.first.offset3p << "\t"
                  << node.second << std::endl;
    }


    std::cerr << "BedRecords: " << counter << std::endl;


    // destructors
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}
