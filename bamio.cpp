#include "bamio.h"

BamIO::BamIO() :
    m_mapq(0),
    m_readLength(0),
    m_handleBam(nullptr),
    m_handleBai(nullptr),
    m_handleHeader(nullptr),
    m_handleIterator(nullptr),
    m_alignment(nullptr)
{}

BamIO::BamIO(const std::string &fileName, uint8_t _mapq, uint16_t _readLength) :
    m_mapq(0),
    m_readLength(0),
    m_handleBam(nullptr),
    m_handleBai(nullptr),
    m_handleHeader(nullptr),
    m_handleIterator(nullptr),
    m_alignment(nullptr)
{
    open(fileName, _mapq, _readLength);
}

BamIO::~BamIO() {
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

bool BamIO::isOpen() const {
    return (m_handleBam && m_handleBai && m_handleHeader && m_alignment);
}

bool BamIO::open(const std::string &fileName, uint8_t _mapq, uint16_t _readLength) {

    m_handleBam = hts_open(fileName.c_str(), "r");
    if (!m_handleBam) {
        m_error = "BamIO::error, failed to open BAM file " + fileName;
        return false;
    }

    m_handleBai = hts_idx_load(fileName.c_str(), HTS_FMT_BAI);
    if (!m_handleBai) {
        m_error = "BamIO::error, failed to laod BAI index " + fileName + ".bai";
        return false;
    }

    m_handleHeader = sam_hdr_read(m_handleBam);
    if (!m_handleHeader) {
        m_error = "BamIO::error, failed to parse BAM header";
        return false;
    }

    m_alignment = bam_init1();
    if (!m_alignment) {
        m_error = "BamIO::error, failed to allocate alignment memory";
        return false;
    }

    setMapQ(_mapq);
    setReadLength(_readLength);

    return true;
}


uint32_t BamIO::count() {
    uint32_t counter = 0;
    // read bam file record by record
    while (next())
        counter++;
    return counter;
}

bool BamIO::query(const std::string queryChrom, int queryStart, int queryEnd) {
    int queryTid = bam_name2id(m_handleHeader, queryChrom.c_str());
    m_handleIterator = bam_itr_queryi(m_handleBai, queryTid, queryStart, queryEnd);
    return (m_handleIterator == nullptr);
}

bool BamIO::next() {
    int ret = -1;
    while (1) {

        // read next iterator or next sam line
        ret = m_handleIterator ?
                    sam_itr_next(m_handleBam, m_handleIterator, m_alignment) :
                    sam_read1(m_handleBam, m_handleHeader, m_alignment);

        // invalid record
        if ( ret < 0) break;

        // skip special
        if ( m_alignment->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;

        // filter by maping quality
        if (m_mapq > 0) {
            if ( static_cast<uint8_t>(m_alignment->core.qual) < m_mapq ) continue;
        }

        // filter by minimum read length
        if (m_readLength > 0) {
            const uint8_t *p_cigar = m_alignment->data + m_alignment->core.l_qname;
            if (bam_cigar2qlen(static_cast<int>(m_alignment->core.n_cigar), reinterpret_cast<const uint32_t *>(p_cigar)) < m_readLength) continue;
        }

        // successful read
        break;
    }
    return ret > 0;
}


