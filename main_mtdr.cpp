#include <iostream>
#include <memory>
#include <fstream>
#include <numeric>

#include <htslib/faidx.h>


#include "argumentparser.h"
#include "bedrecord.h"
#include "bamhandle.h"
#include "aminoacidtable.h"
#include "version.h"

//void readCodonParameters(AminoAcidTable &aainfo, const std::string &fileName);
//void normalizedFootprintCoveragePerCodon(const std::vector<int> &fc, double afc, const char *sequence, int qStart, int qEnd);
//void calculateMTDR(const std::string &name, const AminoAcidTable &aainfo, const std::vector<int> &fc, double afc, const char *sequence, int qStart, int qEnd);




int main_mtdr(const int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::string fileTable;
    std::vector<std::string> filesBam;
    std::vector<std::shared_ptr<BamHandle>> handlesBam;
    int skipCodons;

    auto p = ArgumentParser("mtdr", std::string(VERSION), "calculates mean translation decoding rate");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentRequired("sequence").setKeyShort("-f").setKeyLong("--fasta").setHelp("FASTA file containing transcript sequence");
    p.addArgumentRequired("footprints").setKeyShort("-b").setKeyLong("--bam").setHelp("BAM file containing footprint coverage").setCount(-1);
    p.addArgumentOptional("aatable").setKeyShort("-t").setKeyLong("--table").setHelp("Amino Acid table that contains average decoding rate").setDefaultValue<std::string>("");
    p.addArgumentOptional("skip").setKeyShort("-s").setKeyLong("--skip").setDefaultValue<int>(10).setHelp("number of codons to skip after/before start/stop codon");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        fileFasta = p.get<std::string>("sequence");
        filesBam = p.get<std::vector<std::string>>("footprints");
        fileTable = p.get<std::string>("aatable");
        skipCodons = p.get<int>("skip");
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // open BAM handles
    for (const auto &fileName : filesBam) {
        std::shared_ptr<BamHandle> handle = std::make_shared<BamHandle>(fileName, 255, 0);
        handlesBam.emplace_back(handle);
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::pausing::error, failed to load FASTA file " << fileFasta << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::pausing::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // loop over BED record
    auto aatable = AminoAcidTable();
    if (!fileTable.empty()) {
        aatable.load(fileTable);
        //aatable.write();
    }

    std::string line;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::stringstream iss(line);
        iss >> bed;

        // filter by ORF length
        int lengthORF = bed.cdsSpan / 3;
        if (lengthORF <= (2 * skipCodons + 1)) continue;

        // accumulate coverage from all files
        std::vector<int> codons(static_cast<std::size_t>(lengthORF), 0);
        for (auto handle : handlesBam)
            handle->calculateSiteCoverage(codons, bed.transcript, 0, bed.span, bed.cdsStart, true);

        // calculate background
        double background = static_cast<double>(std::accumulate(codons.begin() + skipCodons, codons.end() - skipCodons, 0)) / (lengthORF - 2*skipCodons);
        if (background < 0.1) continue;

        // codon NFC or MTDR
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        std::vector<int>::const_iterator it;
        double score_observed = 0.0;
        double score_expected = 0.0;
        double score_paused = 0.0;
        int score_count = 0;
        for (it = codons.begin() + skipCodons; it != (codons.end() - skipCodons); ++it) {

            // skip empty codons
            int fc = (*it);
            if (fc <= 0) continue;

            // retrieve current codon
            int idx_codon = static_cast<int>(it - codons.begin());
            int idx_nucleotide = idx_codon * 3 + bed.cdsStart;
            char tag[4];
            std::strncpy(tag, &sequence[idx_nucleotide], 3);
            tag[3] = '\0';

            auto codon = std::string(tag);
            if (codon.find('N') != std::string::npos)
                continue;

            if (fileTable.empty()) { // NFC
                //std::cout << "\t" << codon << "\t" << (fc / background) << std::endl;
                std::cout << idx_codon << "\t" << codon << "\t" << fc << std::endl;
            }
            else { // MTDR
                score_observed += std::log(fc / background);
                score_expected += std::log(aatable.value(codon));
                score_paused += std::log(aatable.score(codon));
                score_count++;
            }

        }

        if (!fileTable.empty()) {
            score_observed = std::exp(score_observed / score_count);
            score_expected = std::exp(score_expected / score_count);
            score_paused = std::exp(score_paused / score_count);
            std::cout << bed.gene << "\t" << bed.cdsSpan / 3 << "\t" << score_observed << "\t" << score_expected << "\t" << score_paused << std::endl;
        }

        if (sequence)
            free(sequence);

        break;

    }

    // destructors
    fhBed.close();
    fai_destroy(fhFai);

    return 0;
}
    /*
    bool calculateNFC = true;
    std::string fileBed;
    std::string fileFasta;
    std::string fileParams;
    std::vector<BamHandle*> handlesBam;
    AminoAcidTable aainfo;
    
    
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



void calculateMTDR(const std::string &name, const AminoAcidTable &aainfo, const std::vector<int> &fc, double afc, const char *sequence, int qStart, int qEnd)
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

        double timeValue = aainfo.value(codon);
        fast++;


        if ((codonAFC / afc) > 3.0) {
            timeValue = aainfo.timePausing(codon);
            slow++;
            fast--;
        }

        
        // geometric mean accumulation
        logSum += std::log(timeValue);

        counter++;
    }

    // geometric mean
    if (counter > 0)
        std::cout << name << "\t" << counter << "\t" << fast << "\t" << slow << "\t" << 1/std::exp(logSum / counter) << std::endl;
}



void readCodonParameters(AminoAcidTable &aainfo, const std::string &fileName)
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
        aainfo.addParameters(codon, count);
        aainfo.setValue(codon, timeDecoding);
        aainfo.setScore(codon, timePausing);
    }

    fhs.close();
}
*/
