#ifndef BAMHANDLE_H
#define BAMHANDLE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

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
    void countUniqueReads();
    int readsPerRegion(const std::string &qName, int qStart, int qEnd);
    int calculateGCcontent(double &gc_mean, double &gc_std);
    void calculateGCperORF(double &gc_mean, double &gc_M2, double &gc_var, int &readCount);
    void calculateBaseContent(int bufferLength, std::vector<int> &base_A, std::vector<int> &base_C, std::vector<int> &base_G, std::vector<int> &base_T, std::vector<int> &base_N);

    int reads() const {return m_reads;}

private:
    int m_mapq;
    int m_length;
    int m_reads;
    std::string m_file;
    htsFile *m_bam;
    hts_idx_t *m_bai;
    bam_hdr_t *m_header;
    hts_itr_t *m_iterator;
};

#endif
