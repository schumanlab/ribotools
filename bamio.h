#ifndef BAMIO_H
#define BAMIO_H

#include <string>
#include <vector>
#include <memory>

#include <htslib/hts.h>
#include <htslib/sam.h>

class BamAuxiliary
{
public:
    explicit BamAuxiliary(const std::string &fileBam, int mapq = 0, int minlen = 0);
    ~BamAuxiliary();

    bool isOpen() const;
    const std::string name() const;
    const std::string what() const;
    bool next();
    bool query(const std::string &queryChrom, int queryStart, int queryEnd);
    int count();
    int readStart() const;
    int readEnd() const;
    int readLength() const;
    std::vector<int> depth(const std::string &queryChrom, int queryStart, int queryEnd);

private:

    struct bam_aux_t {
        htsFile *fp;
        bam_hdr_t *hdr;
        hts_idx_t *idx;
        hts_itr_t *iter;
        int min_mapQ, min_len;
    };

    bam_aux_t m_aux;
    bam1_t *m_bam;
    std::string m_name;
    std::string m_error;

    static int read_bam(void *data, bam1_t *b);
};


class BamIO
{
public:
    explicit BamIO(const std::string &fileBam, int mapq = 0, int minlen = 0) :
        m_isOpen(true),
        m_error("")
    {
        std::vector<std::string> files;
        files.push_back(fileBam);
        BamIO(files, mapq, minlen);
    }

    explicit BamIO(const std::vector<std::string> &filesBam, int mapq = 0, int minlen = 0) :
        m_isOpen(true),
        m_error("")
    {
        for (auto fileName : filesBam) {
            auto handle = std::make_shared<BamAuxiliary>(fileName, mapq, minlen);
            if (!handle->isOpen()) {
                m_isOpen = false;
                m_error = handle->what();
                return;
            }
            aux.emplace_back(handle);
        }
    }


    bool isOpen() const {return m_isOpen;}
    const std::string what() const {return m_error;}

    std::vector<std::shared_ptr<BamAuxiliary>> aux;

private:
    bool m_isOpen;
    std::string m_error;
};



/*

class BamIO
{
public:
    explicit BamIO();
    explicit BamIO(const std::string &fileName, uint8_t _mapq = 0, uint16_t _readLength = 0);
    ~BamIO();

    bool isOpen() const;
    bool open(const std::string &fileName, uint8_t _mapq = 0, uint16_t _readLength = 0);
    const std::string name() const;

    uint8_t mapq() const {return m_mapq;}
    void setMapQ(uint8_t value) {m_mapq = value;}

    int32_t readLength() const {
        const uint8_t *p_cigar = m_alignment->data + m_alignment->core.l_qname;
        return bam_cigar2qlen(static_cast<int>(m_alignment->core.n_cigar), reinterpret_cast<const uint32_t *>(p_cigar));
    }
    int32_t readStart() const {return m_alignment->core.pos;}
    int32_t readEnd() const {return m_alignment->core.pos + (m_alignment->core.n_cigar ? bam_cigar2rlen(m_alignment->core.n_cigar, bam_get_cigar(m_alignment)) : 1);}
    void setReadLength(uint16_t value) {m_readLength = value;}
    int32_t count();
    bool query(const std::string &queryChrom, int queryStart, int queryEnd);
    bool next();
    const std::string error() const {return m_error;}
    std::vector<int> depth(const std::string &queryChrom, int queryStart, int queryEnd);

private:
    uint8_t m_mapq;
    uint16_t m_readLength;
    std::string m_name;
    std::string m_error;
    htsFile *m_handleBam;
    hts_idx_t *m_handleBai;
    bam_hdr_t *m_handleHeader;
    hts_itr_t *m_handleIterator;
    bam1_t *m_alignment;







};



*/

#endif // BAMIO_H
