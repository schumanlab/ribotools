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
    double edge_min = -0.1;
    double edge_max = 1.1;
    int bin_counts = 100;
    double bin_width = (edge_max - edge_min) / (bin_counts - 1);
    int buffer_size = bin_counts * handlesBam.size();
    std::vector<double> metageneSum(buffer_size, 0.0);
    std::vector<int> metageneNorm(buffer_size, 0);

    // loop over BED records
    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;
        bed.parseExons();

        int codonsCDS = bed.cdsSpan / 3;

        // calculate coverage per file
        int f = 0;
        for (auto handle : handlesBam) {
            
            std::map<int, int> depth;
            std::vector<double> histSum(bin_counts, 0.0);
            std::vector<int> histNorm(bin_counts, 0);
            
            // retrieve codon depth
            handle->codonDepth(depth, bed.transcript, bed.span, bed.cdsStart);
            
            // calculate average coverage in ORF
            double averageDepth = 0.0;
            int averageNorm = 0;
            for (auto value : depth) {
                if ((0 <= value.first) && (value.first <= codonsCDS)) {
                    averageDepth += static_cast<double>(value.second);
                    averageNorm++;
                }
            }
            averageDepth /= averageNorm;
            
            if (averageDepth > 1.0) {

                for (auto value : depth) {
                    double x = static_cast<double>(value.first) / (bed.cdsSpan/3);
                    int xbin = static_cast<int>((x - edge_min) / bin_width);
                    if ((0 <= xbin) && (xbin < bin_counts)) {
                        histSum[xbin] += value.second / averageDepth;
                        histNorm[xbin]++;
                    }
                    //std::cout << xbin << "\t" << x << "\t" << value.first << "\t" << value.second << "\t" << value.second / averageDepth << std::endl;
                }
                
                for (int k = 0; k < bin_counts; ++k) {

                    double histAverage = 0.0;
                    int idx = f + k;
                    if (histNorm[k] != 0) {
                        histAverage = histSum[k] / histNorm[k];
                        metageneNorm[idx]++;
                    }

                    metageneSum[idx] += histAverage;
                }

            }

            f += bin_counts;
        }
        
    }

    int k = 0;
    for (int f = 0; f < metageneSum.size(); ++f) {
        
        std::cout << k << "\t" << metageneSum[f] << "\t" << metageneNorm[f] << std::endl;
        k++;
        if (k >= bin_counts)
            k = 0;
    }




    // destructors
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}