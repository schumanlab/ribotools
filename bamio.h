#ifndef BAMIO_H
#define BAMIO_H

#include <iostream>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

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
    void setReadLength(uint16_t value) {m_readLength = value;}
    uint32_t count();
    bool query(const std::string &queryChrom, int queryStart, int queryEnd);
    bool next();
    const std::string error() const {return m_error;}
    void depth(std::vector<int> &data, const std::string &queryChrom, int queryStart, int queryEnd);

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


    struct mplp_aux_t {
        htsFile *fp;
        bam_hdr_t *hdr;
        hts_itr_t *iter;
        int min_mapQ, min_len;
    };


    static int read_bam(void *data, bam1_t *b) {
        mplp_aux_t *aux = static_cast<mplp_aux_t*>(data);
        int ret;
        while (1) {
            ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
            if ( ret<0 ) break;
            if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
            if ( static_cast<int>(b->core.qual) < aux->min_mapQ ) continue;
            //if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue; // filter based on read length
            break;
        }
        return ret;
    }

};





#endif // BAMIO_H
