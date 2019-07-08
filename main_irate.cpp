#include <iostream>
#include <fstream>
#include <unordered_map>

#include "parserargv.h"
#include "bedrecord.h"
#include "bamhandle.h"

void readTranslationEfficiency(std::unordered_map<std::string, double> &te_map, const std::string &fileName);

int main_irate(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileTE;
    std::vector<BamHandle*> handlesBam;
    std::unordered_map<std::string, double> te_map;

    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::irate::error, provide BED file." << std::endl;
        return 1;
    }

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            auto handle = new BamHandle(fileNameNext, 255, 0);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::irate::error, provide BAM file." << std::endl;
        return 1;
    }

    if (!(parser.find("--te") && parser.next(fileTE))) {
        std::cerr << "ribotools::irate::error, provide TE file." << std::endl;
        return 1;
    }
    readTranslationEfficiency(te_map, fileTE);

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "codonNFC::Error, failed to open BED reference" << std::endl;
        return 1;
    }

    // loop over BED records
    int countRecords = 0;
    const double Xi = 0.015; //0.2742;
    std::string line;
    while (std::getline(fhBed, line)) {

        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // translation efficiency
        auto record = te_map.find(bed.gene);
        if (record == te_map.end()) continue;
        double te = record->second;
        
        // calculate footprint coverage
        std::vector<int> fc(bed.cdsSpan, 0);
        for (auto handle : handlesBam)
            handle->calculateFootprintCoverage(fc, bed.transcript, bed.cdsStart, bed.cdsEnd);

        // calculate Nc
        int Nc = bed.cdsSpan / 3;

        // average ribosomes per codon
        int cRiboReadsTotal = 0;
        int cRiboFound = 0;
        for (int x : fc) {
            cRiboReadsTotal += x;
            if (x > 0) cRiboFound++;
        }
        if (cRiboReadsTotal == 0) continue;
        
        // average ribosomes in initiation
        int cRiboReadsInitiation = 0;
        for (int k = 3; k < 33; ++k) {
            cRiboReadsInitiation += fc[k];
        }
        if (cRiboReadsInitiation == 0) continue;

        
        // calculate <T>
        double T_exp = (Nc - 1) * (1/4.0); // 4 AA / sec as rate

        // calculate <Ro>
        double Ro_exp = te * Xi;

        // calculate R_o
        double Ro = (static_cast<double>(cRiboReadsInitiation) / cRiboReadsTotal) * Ro_exp * (Nc - 1);




        // calculate initiation rate
        double ir = (Ro_exp * (Nc - 1)) / (T_exp * (1 - Ro));
        std::cout << bed.name << "\t" 
                  << Nc << "\t"
                  << (cRiboFound/3) << "\t" 
                  << cRiboReadsTotal << "\t" 
                  << cRiboReadsInitiation << "\t" 
                  << Ro_exp << "\t"
                  << T_exp << "\t"
                  << Ro << "\t"
                  << ir << std::endl;
        
        countRecords++;
    }

    std::cerr << "Bed Records: " << countRecords << std::endl;

    // destructors
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}


void readTranslationEfficiency(std::unordered_map<std::string, double> &te_map, const std::string &fileName)
{
    std::ifstream fhs;
    fhs.open(fileName);
    if (!fhs.is_open()) {
        std::cerr << "ribotools::irate::error, failed to read translation efficiency " << fileName << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fhs, line)) {
        std::istringstream iss(line);
        std::string gene;
        double rep1, rep2, rep3, repAvg;

        iss >> gene >> rep1 >> rep2 >> rep3 >> repAvg;
    
        te_map[gene] = repAvg;
    }

    fhs.close();
}