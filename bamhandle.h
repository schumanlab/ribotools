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

private:
    htsFile *bam;
    hts_idx_t *bai;
    bam_hdr_t *header;
};

#endif