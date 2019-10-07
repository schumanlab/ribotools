#include <iostream>
//#include <fstream>
#include <memory>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"

int main_gcratio(const int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<std::shared_ptr<BamHandle>> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);

    /*
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::translate::error, provide a BED file." << std::endl;
        return 1;
    }
    */

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            std::shared_ptr<BamHandle> handle = std::make_shared<BamHandle>(fileNameNext, 255, 0);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::gcratio::error, provide BAM file." << std::endl;
        return 1;
    }

    // open BED file
    /*
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::translate::error, failed to open BED reference " << fileBed << std::endl;
        return 1;
    }
    */




    // calculate gc-content per BAM file
    //std::cout << "# name\tread.count\tgc-ratio_mean\tgc-ratio_std" << std::endl;
    for (auto handle : handlesBam) {
        //double gc_mean = 0.0;
        //double gc_M2 = 0.0;
        //double gc_var = 0.0;
        //int readCount = 0;
        uint64_t baseCount = 0;
        uint64_t gcCount = 0;
        handle->calculateGCbases(baseCount, gcCount);
        std::cout << handle->name() << "\t" << baseCount << "\t" << gcCount << std::endl;

        //handle->calculateGCcontent(gc_mean, gc_M2, gc_var, readCount);

        /*
        // loop over bed file
        std::string line;
        while (std::getline(fhBed, line)) {

            // parse bed line
            auto bed = BedRecord();
            std::istringstream iss(line);
            iss >> bed;

            handle->query(bed.transcript, bed.cdsStart + 30, bed.cdsEnd - 30);
            handle->calculateGCperORF(gc_mean, gc_M2, gc_var, readCount);
        }

        // rewind file
        fhBed.clear();
        fhBed.seekg(0);
        */

        //double gc_std = std::sqrt(gc_var);

        //std::cout << handle->name() << "\t" << readCount << "\t" << gc_mean << "\t" << gc_std << std::endl;


    }

    // close bed file
    //fhBed.close();

    return 0;
}
