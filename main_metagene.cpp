#include <iostream>
#include <fstream>
#include <numeric>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"

int main_metagene(int argc, const char *argv[])
{
    std::string fileBed;
    std::vector <BamHandle*> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);

    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::metagene::error, provide BED file." << std::endl;
        return 1;
    }

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            BamHandle *handle = new BamHandle(fileNameNext);
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

    // prepare metagene
    int bufferSize = 100 * handlesBam.size();
    std::vector<double> metagene(bufferSize, 0.0);

    // loop over BED records
    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;
        bed.parseExons();

        //int codonsGene = bed.
        //int codons5pUTR = BedRecord::nextCodon(bed.cdsStart);

        std::cout << bed.gene << std::endl;
        
        for (int i = 38; i < 47; ++i) {
            std::cout << i << "\t" << bed.prevCodon(i) << "\t" << bed.nextCodon(i) << std::endl;
        }
         

        /*

        // calculate coverage per file
        for (auto handle : handlesBam) {
            
            std::vector<double> depth(codonsGene, 0.0);
            handle->codonDepth(depth, bed.transcript, bed.span, bed.cdsStart, codons5pUTR);
            //double orfAverage = std::accumulate(depth.begin() + 14, depth.end() - 10, 0.0);

            int c = 0;
            for (auto value : depth) {
                std::cout << c++ << "\t" << value << std::endl;
            }

        }
         */
        
    }



    // destructors
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}