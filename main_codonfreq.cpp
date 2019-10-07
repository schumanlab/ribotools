#include <iostream>
#include <fstream>
#include <unordered_map>

#include <htslib/faidx.h>

#include "parserargv.h"
#include "bedrecord.h"

struct AminoAcid {
    char letter;
    int count;
    std::string codon;
    std::string code;
    std::string name;
};

void readCodonTable(std::unordered_map<std::string, AminoAcid> &codonMap, const std::string &fileName);
void printCodonTable(std::unordered_map<std::string, AminoAcid> &codonMap, int totalCodons);

int main_codonfreq(const int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::string fileTable;
    std::unordered_map<std::string, AminoAcid> codonMap;

    // parse command line arguments
    ParserArgv parser(argc, argv);
    if (!(parser.find("--table") && parser.next(fileTable))) {
        std::cerr << "ribotools::codonfreq::error, codon table is required." << std::endl;
        return 1;
    }

    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::codonfreq::error, provide a BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::codonfreq::error, provide a FASTA file." << std::endl;
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

    readCodonTable(codonMap, fileTable);
    int totalCodons = 0;

    // loop over BED record
    std::string line;
    while (std::getline(fhBed, line)) {
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;
        
        //std::cout << bed.gene << std::endl;

        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        for (int c = bed.cdsStart; c < bed.cdsEnd; c += 3) {
            char codon[4];
            std::strncpy(codon, &sequence[c], 3);
            codon[3] = '\0';

            if (std::strchr(codon, 'N')) continue;
            
            auto next = codonMap.find(std::string(codon));
            if (next != codonMap.end()) {
                next->second.count++;
                totalCodons++;
            }
            
        }

        if (sequence)
            free(sequence);
    }

    // destructors
    fai_destroy(fhFai);
    fhBed.close();

    printCodonTable(codonMap, totalCodons);
    
    return 0;
}

void printCodonTable(std::unordered_map<std::string, AminoAcid> &codonMap, int totalCodons)
{
    // print codon map
    for (std::pair<std::string, AminoAcid> next : codonMap) {
        std::cout << next.first << "\t" 
                  << next.second.letter << "\t"
                  << next.second.code << "\t"
                  << next.second.name << "\t"
                  << static_cast<double>(next.second.count) / totalCodons << std::endl;
    }
}


void readCodonTable(std::unordered_map<std::string, AminoAcid> &codonMap, const std::string &fileName)
{
    std::ifstream fhs;
    fhs.open(fileName);
    if (!fhs.is_open()) {
        std::cerr << "ribotools::codonfreq::error, failed to read codon table " << fileName << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fhs, line)) {
        std::istringstream iss(line);
        AminoAcid next;
        iss >> next.codon >> next.letter >> next.code >> next.name;
        next.count = 0;
        codonMap[next.codon] = next;
    }

    fhs.close();
}
