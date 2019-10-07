#include <iostream>
#include <fstream>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bedrecord.h"

int main_gcref(const int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;

    // parse command line arguments
    ParserArgv parser(argc, argv);

    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::gcref::error, provide a BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::gcref::error, provide a FASTA file." << std::endl;
        return 1;
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::gcref::error, failed to load fasta reference " << fileFasta << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::gcref::error, failed to open BED reference " << fileBed << std::endl;
        return 1;
    }

    // loop over BED record
    std::string line;
    while (std::getline(fhBed, line)) {

        // parse bed line
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // read fasta sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        int gc_all = 0;
        int gc_fputr = 0;
        int gc_ini = 0;
        int gc_cds = 0;
        int gc_ter = 0;
        int gc_tputr = 0;

        int iniStart = bed.cdsStart;
        int iniEnd = std::min(bed.cdsStart + 30, bed.cdsEnd);

        int terStart = std::max(bed.cdsEnd - 30, bed.cdsStart);
        int terEnd = bed.cdsEnd;

        for (int i = 0; i < bed.span; ++i) {
            char base = sequence[i];

            if (!((base == 'C') || (base == 'G')))
                continue; // no GC match

            gc_all++;

            if (i < bed.cdsStart)
                gc_fputr++;

            if (bed.cdsEnd < i )
                gc_tputr++;

            if ((bed.cdsStart <= i) && (i <= bed.cdsEnd))
                gc_cds++;

            if ((iniStart <= i) && (i <= iniEnd))
                gc_ini++;

            if ((terStart <= i) && (i <= terEnd))
                gc_ter++;
        }

        // calculate gc
        std::cout << bed.name << "\t" <<
                     static_cast<double>(gc_all)/bed.span << "\t" <<
                     static_cast<double>(gc_fputr)/(std::max(bed.cdsStart,1)) << "\t" <<
                     static_cast<double>(gc_ini)/(iniEnd - iniStart) << "\t" <<
                     static_cast<double>(gc_cds)/bed.cdsSpan << "\t" <<
                     static_cast<double>(gc_ter)/(terEnd - terStart) << "\t" <<
                     static_cast<double>(gc_tputr)/(std::max(bed.span - bed.cdsEnd,1)) << std::endl;
    }


    // destructors
    fai_destroy(fhFai);
    fhBed.close();


    return 0;
}
