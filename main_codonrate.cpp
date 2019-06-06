#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include <htslib/faidx.h>

#include "parserargv.hpp"
#include "bedrecord.hpp"

int parseRates(int &elementSize, std::unordered_map<std::string, std::vector<double>> &rateTable, const std::string &fileName);

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
        std::cerr << "ribotools::codonrate::error, provide a BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("--rates") && parser.next(fileRates))) {
        std::cerr << "ribotools::codonrate::error, provide a BED file." << std::endl;
        return 1;
    }


    std::unordered_map<std::string, std::vector<double>> rateTable;
    int elementSize = 0;
    if (parseRates(elementSize, rateTable, fileRates) > 0) {
        std::cerr << "ribotools::codonrate::error, failed to parse rates." << std::endl;
        return 1;
    }

    /*
    std::string testKey = "AAA";
    std::unordered_map<std::string, std::vector<double>>::const_iterator iot = rateTable.find(testKey);
    if (iot != rateTable.end()) {
        std::cout << iot->first;
        
        for (std::vector<double>::const_iterator it = iot->second.begin(); it != iot->second.end(); ++it) {
            std::cout << "\t" << *it;
        }
        std::cout << std::endl;
    }
    */

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
        bed.parseExons();

        int queryStart = bed.cdsStart;
        int queryEnd = bed.cdsEnd;
        int querySpan = queryEnd - queryStart;

        std::vector<double> resultSum(elementSize);
        std::vector<int> resultCount(elementSize);
        int i = 0;
        char *rnaSeq = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        for (int n = 0; n < querySpan; n += 3) {
            char codonSeq[4];
            memset(codonSeq, '\0', 4);
            std::strncpy(codonSeq, &rnaSeq[queryStart + n], 3);
            codonSeq[3] = '\0';

            if (std::strchr(codonSeq, 'N'))
                continue;
            
            std::unordered_map<std::string, std::vector<double>>::const_iterator iot = rateTable.find(std::string(codonSeq));
            if (iot != rateTable.end()) {
                i = 0;
                for (std::vector<double>::const_iterator it = iot->second.begin(); it != iot->second.end(); ++it) {
                    double value = *it;
                    
                    if (value > 0) {
                        resultSum[i] += log(1/value);
                        resultCount[i]++;
                    }
            
                    i++;        
                }
            }
        }

        // write results
        std::cout << bed.name << "\t" << querySpan;
        for (int j = 0; j < elementSize; ++j) {
            std::cout << "\t" << exp(resultSum[j]/resultCount[j]);
        }
        std::cout << std::endl;

        if (rnaSeq)
            free(rnaSeq);
    }

    // destructors
    fai_destroy(fhFai);
    fhBed.close();
    return 0;
}


int parseRates(int &elementSize, std::unordered_map<std::string, std::vector<double>> &rateTable, const std::string &fileName)
{
    std::ifstream fhs;
    std::unordered_map<std::string, std::vector<double>> rate_table;
    fhs.open(fileName);
    if (!fhs.is_open()) {
        std::cerr << "ribotools::codonrate::error, failed to open rates table " << fileName << std::endl;
        return 1;
    }

    std::string line;

    // read header
    std::getline(fhs, line);
    std::istringstream iss(line);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>());
    elementSize = tokens.size() - 4;

    // read each line
    while (std::getline(fhs, line)) {
        
        // split line
        std::istringstream iss(line);
        std::vector<std::string> tokens(std::istream_iterator<std::string>{iss},
                                         std::istream_iterator<std::string>());

        // map key
        std::string key = tokens[0];
        std::vector<double> rates;

        for (std::vector<std::string>::size_type i = 4; i != tokens.size(); i++) {
            double value = std::stod(tokens[i]);
            rates.push_back(value);
        }

        rateTable.insert({key, rates});
    }

    fhs.close();
    return 0;
}