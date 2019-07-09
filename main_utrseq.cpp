#include <iostream>
#include <fstream>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bedrecord.h"

void exportSequence(const char *sequence, const char *name, int start, int end);

int main_utrseq (int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    bool use5pUTR = true;

    // parse command line arguments
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::codonrate::error, provide a BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::codonrate::error, provide a FASTA file." << std::endl;
        return 1;
    }

    if (parser.find("--5putr"))
        use5pUTR = true;

    if (parser.find("--3putr"))
        use5pUTR = false;

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::codonrate::error, failed to load fasta reference " << fileFasta << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::codonrate::error, failed to open BED reference " << fileBed << std::endl;
        return 1;
    }

    // loop over BED record
    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        int span5pUTR = bed.cdsStart;
        int span3pUTR = bed.span - bed.cdsEnd;

        if ((0 < span5pUTR) && use5pUTR) {
            exportSequence(sequence, bed.name.c_str(), std::max(0, bed.cdsStart - 500), bed.cdsStart);
        }

        if ((0 < span3pUTR) && (use5pUTR==false)) {
            exportSequence(sequence, bed.name.c_str(), bed.cdsEnd, std::min(bed.span, bed.cdsEnd + 500));
        }

        if (sequence)
            free(sequence);

    }

    // destructors
    fai_destroy(fhFai);
    fhBed.close();


    return 0;
}


void exportSequence(const char *sequence, const char *name, int start, int end)
{
    int span = end - start;
    char region[span + 1];
    std::strncpy(region, &sequence[start], span);
    region[span] = '\0';

    std::cout << name << "\t" << region << "\n";
    return;
}
