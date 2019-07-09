#include <iostream>
#include <fstream>
#include <unordered_map>
#include <numeric>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parserargv.h"
#include "bedrecord.h"
#include "bamhandle.h"
#include "aminoacids.h"

void readCodonParameters(AminoAcids &aainfo, const std::string &fileName);
void normalizedFootprintCoveragePerCodon(const std::vector<int> &fc, double afc, const char *sequence, int qStart, int qEnd);
void calculateMTDR(const std::string &name, const AminoAcids &aainfo, const std::vector<int> &fc, double afc, const char *sequence, int qStart, int qEnd);

int main_mtdr(int argc, const char *argv[])
{
    bool calculateNFC = true;
    std::string fileBed;
    std::string fileFasta;
    std::string fileParams;
    std::vector<BamHandle*> handlesBam;
    AminoAcids aainfo;
    
    
    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::mfdr::error, provide BED file." << std::endl;
        return 1;
    }
        
    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::mfdr::error, provide FASTA file." << std::endl;
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
        std::cerr << "ribotools::mfdr::error, provide BAM file." << std::endl;
        return 1;
    }

    if ((parser.find("--params") && parser.next(fileParams))) {
        calculateNFC = false;
        readCodonParameters(aainfo, fileParams);
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "codonNFC::Error, failed to open BED reference" << std::endl;
        return 1;
    }


    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "codonNFC::Error, failed to load fasta reference" << std::endl;
        return 1;
    }

    // loop over BED records
    int countRecords = 0;
    std::string line;
    while (std::getline(fhBed, line)) {

        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // calculate footprint coverage
        std::vector<int> fc(bed.span, 0);
        for (auto handle : handlesBam)
            handle->calculateFootprintCoverage(fc, bed.transcript, 0, bed.span);

        // calculate average in CDS
        int qStart = bed.cdsStart;
        int qEnd = bed.cdsEnd;
        int qOffset = 0;
        if (calculateNFC) {
            qOffset = 60; // 20 codons to avoid initiaiton / termination
            qStart += 60; // 20 codons to avoid initiaiton
            qEnd -= 60; // 20 codons to avoid termination

        }
        double afc = std::accumulate(fc.begin()+qOffset, fc.end()-qOffset, 0) / (qEnd - qStart);
        
        // filter based on average
        if (afc == 0.0) continue;
        if (calculateNFC && (afc < 1.0)) continue;
        

        // estimate NFC per codon
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        
        if (calculateNFC) {
            normalizedFootprintCoveragePerCodon(fc, afc, sequence, qStart, qEnd);
        }
        else {
            calculateMTDR(bed.name, aainfo, fc, afc, sequence, qStart, qEnd);
        }
        
        countRecords++;

        if (sequence)
            free(sequence);
        
    }

    std::cerr << "Bed Records: " << countRecords << std::endl;

    // destructors
    fai_destroy(fhFai);
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}

void normalizedFootprintCoveragePerCodon(const std::vector<int> &fc, double afc, const char *sequence, int qStart, int qEnd)
{
    // print codons
    for (int c = qStart; c < qEnd; c += 3) {
        double nfc = (fc[c] + fc[c + 1] + fc[c + 2]) / (3 * afc);
        char codonSeq[4];
        std::strncpy(codonSeq, &sequence[c], 3);
        codonSeq[3] = '\0';
        if (std::strchr(codonSeq, 'N')) continue;
        
        if ((0.0 < nfc) && (nfc <= 10.0))
            std::cout << codonSeq << "\t" << nfc << std::endl;
    }
}



void calculateMTDR(const std::string &name, const AminoAcids &aainfo, const std::vector<int> &fc, double afc, const char *sequence, int qStart, int qEnd)
{
    int counter = 0;
    int fast = 0;
    int slow = 0;
    
    double logSum = 0.0;
    
    for (int c = qStart; c < qEnd; c += 3) {
        //double codonAFC = (fc[c] + fc[c + 1] + fc[c + 2]) / 3.0;
        
        char codonSeq[4];
        std::strncpy(codonSeq, &sequence[c], 3);
        codonSeq[3] = '\0';
        std::string codon(codonSeq);

        if (codon.length() != 3) continue;
        if (std::strchr(codonSeq, 'N')) continue;

        double timeValue = aainfo.timeDecoding(codon);
        fast++;

        /*
        if ((codonAFC / afc) > 3.0) {
            timeValue = aainfo.timePausing(codon);
            slow++;
            fast--;
        }
        */
        
        // geometric mean accumulation
        logSum += std::log(timeValue);

        counter++;
    }

    // geometric mean
    if (counter > 0)
        std::cout << name << "\t" << counter << "\t" << fast << "\t" << slow << "\t" << 1/std::exp(logSum / counter) << std::endl;
}



void readCodonParameters(AminoAcids &aainfo, const std::string &fileName)
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

        std::string codon;
        char letter;
        std::string code;
        std::string name;
        int count;
        double timeDecoding;
        double timePausing;

        iss >> codon >> letter >> code >> name >> count >> timeDecoding >> timePausing;
        aainfo.add(codon, count);
        aainfo.setTimeDecoding(codon, timeDecoding);
        aainfo.setTimePausing(codon, timePausing);
    }

    fhs.close();
}
