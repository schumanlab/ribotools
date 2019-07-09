#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bedrecord.h"

struct AminoAcid {
    char letter;
    double weight;
    std::string codon;
    std::string code;
    std::string name;
};

void readCodonWeights(std::unordered_map<std::string, double> &codonWeights, const std::string &fileName);

int main_codonrate(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::string fileRates;

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

    if (!(parser.find("--rates") && parser.next(fileRates))) {
        std::cerr << "ribotools::codonrate::error, provide a RATE file." << std::endl;
        return 1;
    }


    std::unordered_map<std::string, double> codonWeights;
    readCodonWeights(codonWeights, fileRates);

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

        double logSum = 0.0;
        int norm = 0;

        for (int c = bed.cdsStart; c < bed.cdsEnd; c += 3) {
            char codon[4];
            std::strncpy(codon, &sequence[c], 3);
            codon[3] = '\0';

            if (std::strchr(codon, 'N')) continue;
            
            auto next = codonWeights.find(std::string(codon));
            if (next != codonWeights.end()) {
                logSum += std::log(next->second);
                norm++;
            }            
        }
        
        // write results
        if (norm > 0) {
            std::cout << bed.transcript << "\t" << bed.gene << "\t" << std::exp(logSum / norm) << std::endl;
        }
            
        
        if (sequence)
            free(sequence);
    }

    // destructors
    fai_destroy(fhFai);
    fhBed.close();
    return 0;
}


void readCodonWeights(std::unordered_map<std::string, double> &codonWeights, const std::string &fileName)
{
    std::ifstream fhs;
    
    fhs.open(fileName);
    if (!fhs.is_open()) {
        std::cerr << "ribotools::codonrate::error, failed to open rates table " << fileName << std::endl;
        return;
    }

    // parse file
    std::string line;
    while (std::getline(fhs, line)) {
        
        // skip header
        if (line.at(0) == '#') continue; 

        // split line
        std::istringstream iss(line);
        AminoAcid next;

        iss >> next.codon >> next.letter >> next.code >> next.name >> next.weight;

        codonWeights[next.codon] = next.weight;
        
    }

    fhs.close();
}
