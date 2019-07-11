#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"
#include "aminoacids.h"


int main_pausing(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::vector <BamHandle*> handlesBam;

    const int backgroundWindow_basic = 150;
    const int backgroundWindow_flank = 5;
    const int skipCodons = 10;

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

    // output header
    std::cout << "# name\tcodon\tn.codons\tn.asites\tsum.zscore" << std::endl;

    // loop over BED record
    std::string line;
    int line_counter = 0;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        int lengthCodons = bed.cdsSpan / 3;

        // check if codon length is too small
        if (lengthCodons <= 2*skipCodons + 2*backgroundWindow_flank + 1) continue;

        // calculate A-site codon coverage
        std::vector<int> codons(static_cast<size_t>(lengthCodons), 0);
        for (auto handle : handlesBam)
            handle->calculateSiteCoverage(codons, bed.transcript, 0, bed.span, bed.cdsStart, true);

        // count empty codons
        //int emptyCodons = static_cast<int>(std::count(codons.begin(), codons.end(), 0));
        //double emptyRatio = static_cast<double>(emptyCodons) / lengthCodons;

        // get average coverage
        double background_average = static_cast<double>(std::accumulate(codons.begin(), codons.end(), 0)) / lengthCodons;

        // fiter based on coverage and empty codons
        //if (background_average < 0.1 || emptyRatio > 0.75) continue;
        if (background_average < 0.1) continue;


        // calculate background
        int backgrounWindow_offset = std::min(backgroundWindow_basic, lengthCodons - 2*skipCodons);
        double background_basic = static_cast<double>(std::accumulate(codons.begin() + skipCodons, codons.begin() + backgrounWindow_offset, 0)) / backgrounWindow_offset;

        // retrieve sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        std::vector<int>::const_iterator it;
        std::vector<int>::const_iterator it_prev;
        std::vector<int>::const_iterator it_next;
        auto aa = AminoAcids();

        for (it = codons.begin() + skipCodons; it != (codons.end() - skipCodons); ++it)
        {
            // skip empty codons
            if (*it == 0) continue;

            // calculate previous background
            double background_prev = 0.0;
            it_prev = it - backgroundWindow_flank;
            if (codons.begin() <= it_prev)
                background_prev = static_cast<double>(std::accumulate(it_prev, it, background_prev)) / backgroundWindow_flank;

            // calculate next background
            double background_next = 0.0;
            it_next = it + backgroundWindow_flank + 1;
            if ((it + 1) < codons.end() && (it_next <= codons.end()))
                background_next = static_cast<double>(std::accumulate(it+1, it_next, background_next)) / backgroundWindow_flank;

            // background value
            double background = std::max(background_next, background_prev);
            background = std::max(background, background_basic);

            // skip if background is empty
            if (background == 0.0) continue;

            // zscore
            double zscore = (*it - background) / std::sqrt(background);

            // current codon code
            int idx_codon = static_cast<int>(it - codons.begin());
            int idx_nucleotide = idx_codon * 3 + bed.cdsStart;
            char codon[4];
            std::strncpy(codon, &sequence[idx_nucleotide], 3);
            codon[3] = '\0';

            aa.addTime(std::string(codon), zscore, static_cast<double>(*it));

            // debug output
            //std::cout << idx_codon << "\t" << codon << "\t" << *it << "\t" << background << "\t" << zscore << std::endl;
        }

        // output of AminoAcid map
        aa.log(bed.name);

        if (sequence)
            free(sequence);

        line_counter++;
    }

    std::cerr << "used genes: " << line_counter << std::endl;

    // destructors
    fhBed.close();
    fai_destroy(fhFai);
    for (auto handle : handlesBam)
        delete handle;


    return 0;
}
