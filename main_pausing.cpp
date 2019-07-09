#include <iostream>
#include <fstream>
#include <numeric>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"

int main_pausing(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::vector <BamHandle*> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);

    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::pausing::error, provide BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::pausing::error, provide FASTA file." << std::endl;
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
        std::cerr << "ribotools::pausing::error, provide BAM file." << std::endl;
        return 1;
    }


    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::pausing::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::pausing::error, failed to load fasta reference " << fileFasta << std::endl;
        return 1;
    }

    // loop over BED record
    std::string line;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // calculate CDS coverage
        std::vector<int> fc(bed.cdsSpan/3, 0);
        for (auto handle : handlesBam)
            handle->calculateASiteCoverage(fc, bed.transcript, 0, bed.span, bed.cdsStart);
        
        for(int i = 0; i < fc.size(); ++i)
            std::cout << i << "\t" << fc[i] << std::endl;

        int offset = std::min(150, bed.cdsSpan/3);
        double fcs = static_cast<double>(std::accumulate(fc.begin(), fc.begin() + offset, 0)) / offset;
        std::cout << "SUM " << fcs << std::endl;
        
        // retrieve sequence
        //char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);



        //if (sequence)
        //    free(sequence);
    }

    // destructors
    fhBed.close();
    fai_destroy(fhFai);
    for (auto handle : handlesBam)
        delete handle;


    return 0;
}
