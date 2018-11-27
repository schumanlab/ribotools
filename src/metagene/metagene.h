#ifndef METAGENE_H
#define METAGENE_H

#include <iostream>
#include <sstream>
#include <map>

#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/khash.h>

#include "BedRecord.h"

class MetaGene
{
public:
    MetaGene();
    ~MetaGene();

    void open(const std::string &fileBed, const std::string &fileBam);
    void pileup();

private:
    BGZF *m_fhBed;
    samFile *m_fhBam;
    hts_idx_t *m_fhBai;
    bam_hdr_t *m_header;
    bam1_t *m_bam;

    KHASH_MAP_INIT_INT(32, uint32_t);

    void error(const std::string &errorMessage);
    static inline uint32_t bam_calqlen(const bam1_core_t *core, const uint32_t *cigar);
};

#endif /* METAGENE_H */