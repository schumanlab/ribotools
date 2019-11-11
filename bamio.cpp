#include "bamio.h"

Bam::Bam() :
    m_mapq(0),
    m_readLength(0),
    m_handleBam(nullptr),
    m_handleBai(nullptr),
    m_handleHeader(nullptr),
    m_handleIterator(nullptr),
    m_alignment(nullptr)
{}

Bam::~Bam() {
    if (m_alignment)
        bam_destroy1(m_alignment);

    if (m_handleIterator)
        bam_itr_destroy(m_handleIterator);

    if (m_handleHeader)
        bam_hdr_destroy(m_handleHeader);

    if (m_handleBai)
        hts_idx_destroy(m_handleBai);

    if (m_handleBam)
        hts_close(m_handleBam);
}

bool Bam::isOpen() const {
    return (m_handleBam && m_handleBai && m_handleHeader && m_alignment);
}

bool Bam::open(const std::string &fileName) {
    m_handleBam = hts_open(fileName.c_str(), "r");
    if (!m_handleBam)
        return false;

    m_handleBai = hts_idx_load(fileName.c_str(), HTS_FMT_BAI);
    if (!m_handleBai)
        return false;

    m_handleHeader = sam_hdr_read(m_handleBam);
    if (!m_handleHeader)
        return false;

    m_alignment = bam_init1();
    if (!m_alignment)
        return false;

    return true;
}

void Bam::rewind() {
    // check if iterator is set and reset
    if (m_handleIterator) {
        bam_itr_destroy(m_handleIterator);
        m_handleIterator = nullptr;
    }
}

uint32_t Bam::count() {
    uint32_t counter = 0;
    // read bam file record by record
    while (readBam(m_alignment) > 0)
        counter++;
    return counter;
}

bool Bam::query(const std::string queryChrom, int queryStart, int queryEnd) {
    int queryTid = bam_name2id(m_handleHeader, queryChrom.c_str());
    rewind();
    m_handleIterator = bam_itr_queryi(m_handleBai, queryTid, queryStart, queryEnd);
    return (m_handleIterator == nullptr);
}

int Bam::readBam(bam1_t *b) {
    int ret;
    while (1) {

        // read next iterator or next sam line
        ret = m_handleIterator ?
                    sam_itr_next(m_handleBam, m_handleIterator, b) :
                    sam_read1(m_handleBam, m_handleHeader, b);

        // invalid record
        if ( ret < 0) break;

        // skip special
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;

        // filter by maping quality
        if (m_mapq > 0) {
            if ( static_cast<uint8_t>(b->core.qual) < m_mapq ) continue;
        }

        // filter by minimum read length
        if (m_readLength > 0) {
            const uint8_t *p_cigar = b->data + b->core.l_qname;
            if (bam_cigar2qlen(static_cast<int>(b->core.n_cigar), reinterpret_cast<const uint32_t *>(p_cigar)) < m_readLength) continue;
        }

        // successful read
        break;
    }
    return ret;
}
