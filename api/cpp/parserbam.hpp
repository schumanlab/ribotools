#ifndef PARSERBAM_H
#define PARSERBAM_H

#include <iostream>
#include <htslib/hts.h>
#include <htslib/sam.h>

class ParserBam
{
public:
    explicit ParserBam(const std::string &fileName = "", int mapq = 255);
    ~ParserBam();

    void query(const std::string &chrom, int chromStart, int chromEnd);
    int next();

private:
    htsFile *m_bam;
    bam_hdr_t *m_header;
    hts_idx_t *m_bai;
    hts_itr_t *m_iter;
    bam1_t *m_alignment;
    int m_min_mapQ;
};

#endif /* PARSERBAM_H */