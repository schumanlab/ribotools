#ifndef BAMIO_H
#define BAMIO_H

#include <string>
#include <vector>
#include <memory>

#include <htslib/hts.h>
#include <htslib/sam.h>

class Bam
{
public:
    explicit Bam();
    ~Bam();

    bool isOpen() const;
    bool open(const std::string &fileName);
    uint8_t mapq() const {return m_mapq;}
    void setMapQ(uint8_t value) {m_mapq = value;}
    uint16_t readLength() const {return m_readLength;}
    void setReadLength(uint16_t value) {m_readLength = value;}
    void rewind();
    uint32_t count();
    bool query(const std::string queryChrom, int queryStart, int queryEnd);

private:
    uint8_t m_mapq;
    uint16_t m_readLength;
    std::string m_name;
    htsFile *m_handleBam;
    hts_idx_t *m_handleBai;
    bam_hdr_t *m_handleHeader;
    hts_itr_t *m_handleIterator;
    bam1_t *m_alignment;

    int readBam(bam1_t *b);
};


class BamIO
{
public:
    explicit BamIO() {}
    explicit BamIO(const std::string fileName) {}
    explicit BamIO(const std::vector<std::string> fileNames) {}

    void open(const std::string fileName) {}
    void open(const std::vector<std::string> fileNames) {}

private:
    std::vector<std::shared_ptr<Bam>> m_bams;


};

#endif // BAMIO_H
