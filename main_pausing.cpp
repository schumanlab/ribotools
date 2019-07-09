#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bamhandle.h"
#include "bedrecord.h"


void slidingAverageForward(std::vector<double> &sa, const std::vector<int> &data, int window_size)
{
    int window_elements = 0;
    int window_sum = 0;
    std::vector<double>::iterator it_sa = sa.begin();
    for (std::vector<int>::const_iterator it = data.begin(); it != data.end(); ++it) {
        window_sum += *it;
        window_elements++;
        if (window_elements > window_size) {
            window_elements--;
            window_sum -= *(it - window_size);
        }
        *it_sa = static_cast<double>(window_sum) / window_elements;
        ++it_sa;
    }
}


void slidingAverageReverse(std::vector<double> &sa, const std::vector<int> &data, int window_size)
{
    int window_elements = 0;
    int window_sum = 0;
    std::vector<double>::reverse_iterator it_sa = sa.rbegin();
    for (std::vector<int>::const_reverse_iterator it = data.rbegin(); it != data.rend(); ++it) {
        window_sum += *it;
        window_elements++;
        if (window_elements > window_size) {
            window_elements--;
            window_sum -= *(it - window_size);
        }
        *it_sa = static_cast<double>(window_sum) / window_elements;
        ++it_sa;
    }

}



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

        int countCodons = bed.cdsSpan / 3;

        // calculate A-site codon coverage
        std::vector<int> codons(static_cast<size_t>(countCodons), 0);
        for (auto handle : handlesBam)
            handle->calculateSiteCoverage(codons, bed.transcript, 0, bed.span, bed.cdsStart, true);

        //int i = 0;
        //for(int val : codons)
        //    std::cout << i++ << "\t" << val << std::endl;

        // calculate background
        int spanWindow = std::min(150, countCodons);
        double background_basic = static_cast<double>(std::accumulate(codons.begin(), codons.begin() + spanWindow, 0)) / spanWindow;
        std::vector<double> codons_next(static_cast<size_t>(countCodons), 0.0);
        std::vector<double> codons_prev(static_cast<size_t>(countCodons), 0.0);
        slidingAverageForward(codons_prev, codons, 5);
        slidingAverageReverse(codons_next, codons, 5);


        // retrieve sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        for (int idx_codon = 10; idx_codon < (countCodons - 10); ++idx_codon) {

            // calculate background next
            int idx_next = idx_codon + 1;
            double background_next = (idx_next < countCodons) ? codons_next.at(idx_next) : 0.0;

            // calculate background previous
            int idx_prev = idx_codon - 1;
            double background_prev = (0 < idx_prev) ? codons_prev.at(idx_prev) : 0.0;

            // background value
            double background = std::max(background_next, background_prev);
            background = std::max(background, background_basic);

            // zscore
            double zscore = (codons.at(idx_codon) - background) / std::sqrt(background);


            int idx_nucleotide = idx_codon * 3 + bed.cdsStart;
            char codon[4];
            std::strncpy(codon, &sequence[idx_nucleotide], 3);
            codon[3] = '\0';
            std::cout << idx_codon << "\t" << codon << "\t" << codons.at(idx_codon) << "\t" << background << "\t" << zscore << std::endl;

        }


        if (sequence)
            free(sequence);
    }

    // destructors
    fhBed.close();
    fai_destroy(fhFai);
    for (auto handle : handlesBam)
        delete handle;


    return 0;
}
