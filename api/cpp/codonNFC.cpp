#include <iostream>
#include <fstream>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parsercommands.hpp"
#include "bedrecord.hpp"

int main(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::string fileBam;
    
    // parse command line parameters
    ParserCommands input(argc, argv);

    if (!(input.findOption("-bed") && input.nextArgument(fileBed))) {
        std::cerr << "codonNFC::Error, provide BED file." << std::endl;
        return 1;
    }
        
    if (!(input.findOption("-fasta") && input.nextArgument(fileFasta))) {
        std::cerr << "codonNFC::Error, provide FASTA file." << std::endl;
        return 1;
    }
    
    if (!(input.findOption("-bams") && input.nextArgument(fileBam))) {
        std::cerr << "codonNFC::Error, provide BAM file." << std::endl;
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
    std::string bedLine;
    while (std::getline(fhBed, bedLine)) {

        auto bed = BedRecord();
        std::istringstream iss(bedLine);
        iss >> bed;
        bed.parseExons();
        char *rnaSeq = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        int queryStart = bed.cdsStart;
        int queryEnd = bed.cdsEnd;
        int querySpan = queryEnd - queryStart;
        int queryTid = bam_name2id(hdrBam, bed.transcript().c_str());
        std::cout << bed.transcript() << ":" << queryStart << "-" << queryEnd << std::endl;

        
        hts_itr_t *itrBam = bam_itr_queryi(fhBai, queryTid, queryStart, queryEnd);
        bam1_t *algBam = bam_init1();
        int test = 0;
        int ret = 0;
        while ((ret = sam_itr_next(fhBam, itrBam, algBam)) >= 0) {
            
            test++;
        }
        std::cout << bed.name << "\t" << test << std::endl;

        if (algBam)
            bam_destroy1(algBam);

        if (itrBam)
            bam_itr_destroy(itrBam);

        if (rnaSeq)
            free(rnaSeq);
    }

    // destructors
    bam_hdr_destroy(hdrBam);
    hts_idx_destroy(fhBai);
    hts_close(fhBam);
    fai_destroy(fhFai);
    fhBed.close();
    
    return 0;
}

