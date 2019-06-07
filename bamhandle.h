#ifndef BAMHANDLE_H
#define BAMHANDLE_H

#include <iostream>
#include <vector>
#include <string>
#include <htslib/hts.h>
#include <htslib/sam.h>

class BamHandle
{
public:
    explicit BamHandle(const std::string &fileName);
    ~BamHandle();

    int orfDepth(std::vector<double> &cov, const std::string &chrom, int chromStart, int chromEnd, int orfStart);

private:
    htsFile *bam;
    hts_idx_t *bai;
    bam_hdr_t *header;
};

#endif