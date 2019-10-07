#include <iostream>
#include <fstream>
#include <unordered_map>

#include "parserargv.h"
#include "bedrecord.h"
#include "bamhandle.h"

int main_irate(const int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<BamHandle*> handlesRFP; // ribosome footprints
    std::vector<BamHandle*> handlesRNA; // rna seq

    // parse command line parameters
    ParserArgv parser(argc, argv);

    // parse BED file
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::irate::error, provide BED file." << std::endl;
        return 1;
    }

    // parse BAM files from Ribosome Footprints
    if (parser.find("--rfp")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            auto handle = new BamHandle(fileNameNext, 255, 0);
            handlesRFP.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::irate::error, provide BAM file from ribosome footprints." << std::endl;
        return 1;
    }

    // parse BAM files from RNA seq
    if (parser.find("--rna")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            auto handle = new BamHandle(fileNameNext, 255, 0);
            handlesRNA.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::irate::error, provide BAM file from RNASeq." << std::endl;
        return 1;
    }

    // check if handles size is matching
    if (handlesRFP.size() != handlesRNA.size()) {
        std::cerr << "ribotools::irate::error, provide equal number of BAM files for Ribosome Footprint(--rfp) and RNASeq(--rna)" << std::endl;
        return 1;
    }

    // count reads per file
    std::cerr << "# counting unique reads in ribosome footprint BAMs" << std::endl;
    for (auto handle : handlesRFP) {
        handle->countUniqueReads();
        std::cerr << handle->name() << "\t" << handle->reads() << std::endl;
    }

    std::cerr << "# counting unique reads in RNASeq BAMs" << std::endl;
    for (auto handle : handlesRNA) {
        handle->countUniqueReads();
        std::cerr << handle->name() << "\t" << handle->reads() << std::endl;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "codonNFC::Error, failed to open BED reference" << std::endl;
        return 1;
    }

    // loop over bed records
    std::string line;
    while (std::getline(fhBed, line)) {

        // parse bed line
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // stats per file
        std::vector<BamHandle*>::iterator iteratorRFP;
        std::vector<BamHandle*>::iterator iteratorRNA;

        std::cout << bed.gene << "\t" << bed.span << "\t" << bed.cdsSpan;

        for (iteratorRFP = handlesRFP.begin(), iteratorRNA = handlesRNA.begin();
             (iteratorRFP != handlesRFP.end()) && (iteratorRNA != handlesRNA.end());
             ++iteratorRFP, ++iteratorRNA) {

            std::vector<int> rfpc(static_cast<size_t>(bed.cdsSpan), 0);
            int reads_rna = (*iteratorRNA)->readsPerRegion(bed.transcript, 0, bed.span);

            (*iteratorRFP)->calculateSiteCoverage(rfpc, bed.transcript, bed.cdsStart, bed.cdsEnd, bed.cdsStart, true);
            int reads_rfp = 0;
            int reads_ini = 0;
            for (size_t k = 3; k < rfpc.size(); ++k) {

                reads_rfp += rfpc.at(k);

                if (k < 33)
                    reads_ini += rfpc.at(k);

            }


            //if ((reads_rfp > 0) && (reads_ini > 0) && (reads_rna > 0))
                std::cout << "\t" << (*iteratorRFP)->reads() << "\t" << reads_rfp << "\t" << reads_ini << "\t" <<
                (*iteratorRNA)->reads() << "\t" << reads_rna;
        }
        std::cout << std::endl;


    }


    // close bed file
    fhBed.close();

    // destroy ribosome footprint handles
    for (auto handle : handlesRFP)
        delete handle;

    // destroy RNASeq handles
    for (auto handle : handlesRNA)
        delete handle;

    return 0;
}

/*
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
*/
