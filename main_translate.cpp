#include <iostream>
#include <fstream>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bedrecord.h"
#include "aminoacidtable.h"

int main_translate(const int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;

    // parse command line arguments
    ParserArgv parser(argc, argv);

    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::translate::error, provide a BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::translate::error, provide a FASTA file." << std::endl;
        return 1;
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::translate::error, failed to load fasta reference " << fileFasta << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::translate::error, failed to open BED reference " << fileBed << std::endl;
        return 1;
    }

    // loop over BED record
    AminoAcidTable aainfo;
    std::string line;
    while (std::getline(fhBed, line)) {

        // parse bed line
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // read fasta sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        // translate codons
        std::cout << ">" << bed.name << std::endl;
        int counter = 0;
        for (int c = bed.cdsStart; c < bed.cdsEnd; c += 3) {
            char codon_seq[4];
            std::strncpy(codon_seq, &sequence[c], 3);
            codon_seq[3] = '\0';
            std::string codon(codon_seq);
            char AA = aainfo.translate(codon);
            std::cout << AA;
            counter++;
            if (counter >= 80) {
                std::cout << std::endl;
                counter = 0;
            }
        }
        std::cout << std::endl;
    }


    // destructors
    fai_destroy(fhFai);
    fhBed.close();


    return 0;
}
