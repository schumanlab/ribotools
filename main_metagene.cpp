#include <iostream>
#include <fstream>
#include <numeric>
#include <map>

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
            BamHandle *handle = new BamHandle(fileNameNext, 255, 0);
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
    double edge_min = -0.1;
    double edge_max = 1.1;
    int bin_counts = 100;
    double bin_width = (edge_max - edge_min) / (bin_counts - 1);
    std::vector<double> metageneSum(bin_counts, 0.0);
    std::vector<int> metageneNorm(bin_counts, 0);

    // loop over BED records
    int count_genes = 0;
    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // calculate footprint coverage
        std::vector<int> fc(bed.span, 0);
        for (auto handle : handlesBam)
            handle->calculateFootprintCoverage(fc, bed.transcript, 0, bed.span);

        // average coverage in ORF
        int countReadsORF = 0;
        for (int k = bed.cdsStart; k < bed.cdsEnd; ++k)
            countReadsORF += fc[k];
        double afc = static_cast<double>(countReadsORF) / bed.cdsSpan;
        if (afc < 0.1) continue;

        // metagene coverage
        std::vector<double> histSum(bin_counts, 0.0);
        std::vector<int> histNorm(bin_counts, 0);
        for (int k = 0; k < bed.span; ++k) {
            double idx = static_cast<double>(k - bed.cdsStart) / bed.cdsSpan;
            idx = (idx - edge_min) / bin_width;
            int xbin = static_cast<int>(idx);
            if ((0 <= xbin) && (xbin < bin_counts)) {
                histSum[xbin] += (fc[k] / afc);
                histNorm[xbin]++;
            }
        }

        // add to main 
        for (int k = 0; k < bin_counts; ++k) {

            if (histNorm[k] > 0) {
                metageneSum[k] += (histSum[k] / histNorm[k]);
                metageneNorm[k]++;
            }

        }

        count_genes++;
    }

    for (int k = 0; k < bin_counts; ++k)
        std::cout << k << "\t" << metageneSum[k] << "\t" << metageneNorm[k] << std::endl;

    std::cerr << "genes: " << count_genes << std::endl;

    // destructors
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}
