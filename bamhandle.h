#ifndef BAMHANDLE_H
#define BAMHANDLE_H

#include <iostream>
#include <vector>
#include <algorithm>

#include <htslib/hts.h>
#include <htslib/sam.h>

class BamHandle
{
public:
    explicit BamHandle(const std::string &fileName, int mapq, int length);
    ~BamHandle();

    std::string name();
    void query(const std::string &queryChrom, int queryStart, int queryEnd);
    int readBam(bam1_t *b);
    void calculateFootprintCoverage(std::vector<int> &fc, const std::string &qName, int qStart, int qEnd);
    void calculateSiteCoverage(std::vector<int> &fc, const std::string &qName, int qStart, int qEnd, int qRef, bool useAsite);

private:
    int m_mapq;
    int m_length;
    std::string m_file;
    htsFile *m_bam;
    hts_idx_t *m_bai;
    bam_hdr_t *m_header;
    hts_itr_t *m_iterator;
};

#endif
