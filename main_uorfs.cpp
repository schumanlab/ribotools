#include <iostream>
#include <fstream>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "parserargv.h"
#include "bedrecord.h"

bool uORFs_checkStart(const std::string &codonSeq);
bool uORFs_checkStop(const std::string &codonSeq);

int countReadsOverRegion(const std::string &chrom, int chromStart, int chromEnd, bam_hdr_t *hdrBam, hts_idx_t *fhBai, htsFile *fhBam);

int main_uorfs(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::string fileBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (!(parser.find("--bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::uorfs::error, provide BED file." << std::endl;
        return 1;
    }
        
    if (!(parser.find("--fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::uorfs::error, provide FASTA file." << std::endl;
        return 1;
    }
    
    if (!(parser.find("--bam") && parser.next(fileBam))) {
        std::cerr << "ribotools::uorfs::error, provide BAM file." << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::uorfs::error, failed to open BED reference" << std::endl;
        return 1;
    }


    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::uorfs::error, failed to load fasta reference" << std::endl;
        return 1;
    }


    // open BAM file
    htsFile *fhBam = hts_open(fileBam.c_str(), "r");
    if (!fhBam) {
        std::cerr << "ribotools::uorfs::error, failed to open BAM file " << fileBam << std::endl;
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

    std::cout << "name\tstart\tstop\tcodons\tframe\treads\tORF.start\tORF.reads\tseq" << std::endl;

    // loop over BED records
    std::string line;
    while (std::getline(fhBed, line)) {

        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;
        bed.parseExons();
        int countORF = countReadsOverRegion(bed.transcript, bed.cdsStart, bed.cdsEnd, hdrBam, fhBai, fhBam);
        if (countORF == 0)
            continue;

        char *rnaSeq = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        for (int frame = 0; frame < 3; ++frame) {
            bool flagStart = false;
            int uorfStart = 0;
            for (int b = frame; b < bed.cdsStart; b += 3) {
                char codonSeq[4];
                std::strncpy(codonSeq, &rnaSeq[b], 3);
                codonSeq[3] = '\0';

                if (uORFs_checkStart(std::string(codonSeq))) {
                    flagStart = true;
                    uorfStart = b;
                }

                if (flagStart && uORFs_checkStop(std::string(codonSeq))) {
                    flagStart = false;
                    int uorfEnd = b + 3;
                    int uorfSpan = uorfEnd - uorfStart;
                    char *uorfSeq = (char *)malloc((uorfSpan + 1)*sizeof(char));
                    std::strncpy(uorfSeq, &rnaSeq[uorfStart], uorfSpan);
                    uorfSeq[uorfSpan] = '\0';

                    int countUORF = countReadsOverRegion(bed.transcript, uorfStart, uorfEnd - 3, hdrBam, fhBai, fhBam);

                    if ((uorfSpan > 6) && (countUORF > 0))
                        std::cout << bed.name << "\t" << uorfStart << "\t" << uorfEnd << "\t" << uorfSpan / 3 << "\t" << frame << "\t" << countUORF << "\t" << bed.cdsStart << "\t" << countORF << "\t" << uorfSeq << std::endl;

                    if (uorfSeq)
                        free(uorfSeq);
                }
            }
        }

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


bool uORFs_checkStart(const std::string &codonSeq)
{
    if (codonSeq == "ATG")
        return true;
    return false;
}


bool uORFs_checkStop(const std::string &codonSeq)
{
    if (codonSeq == "TAA")
        return true;
    
    if (codonSeq == "TAG")
        return true;

    if (codonSeq == "TGA")
        return true;

    return false;
}


int countReadsOverRegion(const std::string &chrom, int chromStart, int chromEnd, bam_hdr_t *hdrBam, hts_idx_t *fhBai, htsFile *fhBam)
{
    int chromTid = bam_name2id(hdrBam, chrom.c_str());
    hts_itr_t *itrBam = bam_itr_queryi(fhBai, chromTid, chromStart, chromEnd);
    bam1_t *algBam = bam_init1();
    int ret = 0;
    int count = 0;
    while((ret = sam_itr_next(fhBam, itrBam, algBam)) >= 0)
        count++;
    
    if (algBam)
        bam_destroy1(algBam);

    if (itrBam)
        bam_itr_destroy(itrBam);
    
    return count;
}