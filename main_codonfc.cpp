#include <iostream>
#include <fstream>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parserargv.h"
#include "bedrecord.h"
#include "bamhandle.h"

void calculateFootprintCoverage(std::vector<int> &fc, BamHandle *handle, const std::string &qName, int qStart, int qEnd);
double calculateAverageFootprintCoverage(const std::vector<int> &fc, int qStart, int qEnd);
void normalizedFootprintCoveragePerCodon(const std::vector<int> &fc, double fcAverage, char *sequence, int qStart, int qEnd);

int main_codonfc(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::vector<BamHandle*> handlesBam;
    
    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::codonfc::error, provide BED file." << std::endl;
        return 1;
    }
        
    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::codonfc::error, provide FASTA file." << std::endl;
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
        std::cerr << "ribotools::metagene::error, provide BAM file." << std::endl;
        return 1;
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
    std::string line;
    while (std::getline(fhBed, line)) {

        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;
        bed.parseExons();

        // calculate footprint coverage
        std::vector<int> fc(bed.span, 0);
        for (auto handle : handlesBam)
            calculateFootprintCoverage(fc, handle, bed.transcript, 0, bed.span);

        // calculate average in CDS
        double fcAverage = calculateAverageFootprintCoverage(fc, bed.cdsStart, bed.cdsEnd);
        if (fcAverage < 1.0) continue;

        // estimate NFC per codon
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        
        normalizedFootprintCoveragePerCodon(fc, fcAverage, sequence, bed.cdsStart, bed.cdsEnd);

        if (sequence)
            free(sequence);
    }


    // destructors
    fai_destroy(fhFai);
    fhBed.close();
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}

void normalizedFootprintCoveragePerCodon(const std::vector<int> &fc, double fcAverage, char *sequence, int qStart, int qEnd)
{
    
    // print codons
    for (int c = qStart; c < qEnd; c += 3) {
        double nfc = (fc[c] + fc[c + 1] + fc[c + 2]) / (3 * fcAverage);
        char codonSeq[4];
        std::strncpy(codonSeq, &sequence[c], 3);
        codonSeq[3] = '\0';
        if (std::strchr(codonSeq, 'N')) continue;
        
        //if ((0.0 < nfc) && (nfc <= 10.0))
            std::cout << codonSeq << "\t" << nfc << std::endl;
    }

}


double calculateAverageFootprintCoverage(const std::vector<int> &fc, int qStart, int qEnd)
{
    qStart = std::max(0, qStart);
    qEnd = std::min(static_cast<int>(fc.size()), qEnd);
    if (qEnd - qStart < 1) return 0.0;

    int sum = 0;
    for (int k = qStart; k < qEnd; k++)
        sum += fc[k];
    
    return static_cast<double>(sum) / (qEnd - qStart);
}


void calculateFootprintCoverage(std::vector<int> &fc, BamHandle *handle, const std::string &qName, int qStart, int qEnd)
{
    bam1_t *alignment = bam_init1();
    handle->query(qName, qStart, qEnd);

    while (handle->readBam(alignment) > 0) {

        int readStart = alignment->core.pos;
        int readLength = bam_cigar2qlen(alignment->core.n_cigar, bam_get_cigar(alignment));

        // accumulate P-site per read
        int readPsite = readStart + readLength - 17 - qStart;
        if ((0<= readPsite) && (readPsite < fc.size()))
            fc[readPsite]++;
    }
    
    if (alignment)
        bam_destroy1(alignment);

}
