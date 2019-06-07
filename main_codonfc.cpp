#include <iostream>
#include <fstream>
#include <numeric>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parserargv.h"
#include "bedrecord.h"

void calculateFootprintCoverage(std::vector<double> &fc, BedRecord &bed, bam_hdr_t *hdrBam, hts_idx_t *fhBai, htsFile *fhBam);
void normalizedFootprintCoveragePerCodon(std::vector<double> &fc, double fcAverage, int offset, BedRecord &bed, faidx_t *fhFai);

int main_codonfc(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::string fileBam;
    
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
    
    if (!(parser.find("--bam") && parser.next(fileBam))) {
        std::cerr << "ribotools::codonfc::error, provide BAM file." << std::endl;
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


    // open BAM file
    htsFile *fhBam = hts_open(fileBam.c_str(), "r");
    if (!fhBam) {
        std::cerr << "codonNFC::Error, failed to open BAM file " << fileBam << std::endl;
        return 1;
    }


    // load BAI file
    hts_idx_t *fhBai = hts_idx_load(fileBam.c_str(), HTS_FMT_BAI);
    if (!fhBai) {
        std::cerr << "codonNFC::Error, failed to load BAM index " << fileBam << ".bai" << std::endl;
        return 1;
    }


    // read BAM header
    bam_hdr_t *hdrBam = sam_hdr_read(fhBam);
    if (!hdrBam) {
        std::cerr << "codonNFC::Error, failed to load BAM header " << std::endl;
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
        int offset = 20;
        std::vector<double> fc;
        calculateFootprintCoverage(fc, bed, hdrBam, fhBai, fhBam);
        double fcAverage = std::accumulate(fc.begin() + offset, fc.end() - offset, 0.0) / fc.size();

        // estimate NFC per codon
        normalizedFootprintCoveragePerCodon(fc, fcAverage, offset, bed, fhFai);
    }


    // destructors
    bam_hdr_destroy(hdrBam);
    hts_idx_destroy(fhBai);
    hts_close(fhBam);
    fai_destroy(fhFai);
    fhBed.close();

    return 0;
}

void normalizedFootprintCoveragePerCodon(std::vector<double> &fc, double fcAverage, int offset, BedRecord &bed, faidx_t *fhFai)
{
    int codonSpan = fc.size() - offset;
    char *rnaSeq = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

    // print codons
    for (int c = offset; c < codonSpan; ++c) {
        double nfc = fc[c] / fcAverage;
        char codonSeq[4];
        std::strncpy(codonSeq, &rnaSeq[bed.cdsStart + c * 3], 3);
        codonSeq[3] = '\0';
        if (std::strchr(codonSeq, 'N'))
            continue;
        
        if ((0.0 < nfc) && (nfc <= 10.0))
            std::cout << codonSeq << "\t" << nfc << std::endl;
    }

    if (rnaSeq)
        free(rnaSeq);
}



void calculateFootprintCoverage(std::vector<double> &fc, BedRecord &bed, bam_hdr_t *hdrBam, hts_idx_t *fhBai, htsFile *fhBam)
{
    int queryTid = bam_name2id(hdrBam, bed.transcript().c_str());
    int queryStart = bed.cdsStart;
    int queryEnd = bed.cdsEnd;
    int querySpan = (queryEnd - queryStart) / 3;
    
    hts_itr_t *itrBam = bam_itr_queryi(fhBai, queryTid, queryStart, queryEnd);
    bam1_t *algBam = bam_init1();
    int ret = 0;
    fc.resize(querySpan, 0.0);
    while ((ret = sam_itr_next(fhBam, itrBam, algBam)) >= 0) {
        int readStart = algBam->core.pos;
        int readLength = bam_cigar2qlen(algBam->core.n_cigar, bam_get_cigar(algBam));

        // calculate A-site coverage
        int index = readStart + readLength/2 - queryStart;
        index = index - (index % 3);
        index = index / 3;
        if ((0 <= index) && (index < querySpan))
                fc[index] += 1.0;

        // calculate full coverage
        /*
        for (int k = 0; k < readLength; k += 3) {
            int index = (readStart + k - queryStart);
            index = index - (index % 3);
            index = index / 3;
            if ((0 <= index) && (index < querySpan))
                fc[index] += 1.0;
        }
        */
    }

    if (algBam)
        bam_destroy1(algBam);

    if (itrBam)
        bam_itr_destroy(itrBam);
}