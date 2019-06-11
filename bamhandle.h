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

    void codonDepth(std::vector<double> &depth, const std::string &name, int geneSpan, int cdsStart, int offset);

private:
    htsFile *bam;
    hts_idx_t *bai;
    bam_hdr_t *header;
};

#endif