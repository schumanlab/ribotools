#include <iostream>
#include <fstream>
#include <sstream>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>


#include "bedrecord.hpp"




int main(int argc, char *argv[])
{
    auto fileBed = std::string(argv[1]);
    auto fileBam = std::string(argv[2]);

    std::ifstream fh;
    std::string line;
    
    // constructor
    samFile *fhBam = sam_open(fileBam.c_str(), "r");
    bam_hdr_t *hdrBam = sam_hdr_read(fhBam);
    hts_idx_t *idxBam = sam_index_load(fhBam, fileBam.c_str());
    bam1_t *aln = bam_init1();
    hts_itr_t *iter;
    fh.open(fileBed);

    char *readSeq = (char*)malloc(sizeof(char) * 64);


    while (std::getline(fh, line)) {
        auto bed = BedRecord();
        std::istringstream isline(line);
        isline >> bed;
        bed.parseExons();

        std::cout << bed.transcript() << std::endl;
        int queryStart = static_cast<int>(bed.cdsStart);// + 60;
        int queryEnd = static_cast<int>(bed.cdsEnd);// - 60;


        int tid = bam_name2id(hdrBam, bed.transcript().c_str());
        iter = sam_itr_queryi(idxBam, tid, queryStart, queryEnd);
        int ret = 0;
        while (sam_itr_next(fhBam, iter, aln) >= 0) {
            ret++;
            int readStart = aln->core.pos;
            int readOffset = readStart - bed.cdsStart;
            int readLength = bam_cigar2qlen(aln->core.n_cigar, bam_get_cigar(aln));
            
            for (int b = 0; b < aln->core.l_qseq; b++) {
                readSeq[b] = seq_nt16_str[bam_seqi(bam_get_seq(aln), b)];
            }
            readSeq[aln->core.l_qseq] = '\0';
            std::cout << readStart << "\t" << readOffset << "\t" << readLength << "\t" << readSeq << std::endl;
            
        }
        std::cout << ret << std::endl;


    }

    // destructor
    free(readSeq);
    fh.close();
    bam_destroy1(aln);
    bam_hdr_destroy(hdrBam);
    hts_itr_destroy(iter);
    hts_idx_destroy(idxBam);
    sam_close(fhBam);

    return 0;
}
