#include <iostream>
#include <fstream>
#include <iomanip>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bedrecord.h"

double getGCcontent(char *sequence, int start, int end);

int main_features(const int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;

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
    std::cout << "#transcript\tgene\tlength.Gene\tlength.5pUTR\tlength.CDS\tlength.3pUTR\tgc.5pUTR\tgc.CDS\tgc.3pUTR" << std::endl;

    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        double gc_5pUTR = getGCcontent(sequence, 0, bed.cdsStart);
        double gc_CDS = getGCcontent(sequence, bed.cdsStart, bed.cdsEnd);
        double gc_3pUTR = getGCcontent(sequence, bed.cdsEnd, bed.span);

        std::cout << std::fixed;
        std::cout << std::setprecision(4);
        std::cout << bed.transcript << "\t" 
                  << bed.gene << "\t"
                  << bed.span << "\t"
                  << bed.cdsStart << "\t"
                  << bed.cdsSpan << "\t"
                  << (bed.span - bed.cdsEnd) << "\t"
                  << gc_5pUTR << "\t"
                  << gc_CDS << "\t"
                  << gc_3pUTR << std::endl;

        if (sequence)
            free(sequence);

    }

    // destructors
    fai_destroy(fhFai);
    fhBed.close();
    
    return 0;
}


double getGCcontent(char *sequence, int start, int end)
{
    double gc = 0.0;
    int span = end - start;

    // no sequence
    if (span <= 0)
        return gc;

    char *seqStart = sequence + start;
    char *seqEnd = sequence + end;

    int count_GC = 0;
    while (seqStart < seqEnd) {

        if ((*seqStart == 'G') || (*seqStart == 'C'))
            count_GC++;

        seqStart++;
    }

    return 100.0 * count_GC / span; 
}
