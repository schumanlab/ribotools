#ifndef BAMHANDLE_H
#define BAMHANDLE_H

#include <iostream>
#include <map>

#include <htslib/hts.h>
#include <htslib/sam.h>

class BamHandle
{
public:
    explicit BamHandle(const std::string &fileName);
    ~BamHandle();

    void codonDepth(std::map<int, int> &depth, const std::string &name, int geneSpan, int cdsStart);

    std::string name;

private:
    htsFile *bam;
    hts_idx_t *bai;
    bam_hdr_t *header;
};

#endif