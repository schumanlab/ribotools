#include "bamio.h"

BamIO::BamIO() :
    m_mapq(0),
    m_readLength(0),
    m_name(""),
    m_handleBam(nullptr),
    m_handleBai(nullptr),
    m_handleHeader(nullptr),
    m_handleIterator(nullptr),
    m_alignment(nullptr)
{}

BamIO::BamIO(const std::string &fileName, uint8_t _mapq, uint16_t _readLength) :
    m_mapq(0),
    m_readLength(0),
    m_name(""),
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
    m_name = fileName;

    return true;
}

const std::string BamIO::name() const {
    std::string tag = m_name;

    // remove path
    const size_t idx_path = tag.find_last_of("\\/");
    if (std::string::npos != idx_path)
        tag.erase(0, idx_path + 1);

    // remove extension
    const size_t idx_extension = tag.rfind('.');
    if (std::string::npos != idx_extension)
        tag.erase(idx_extension);

    return tag;
}


uint32_t BamIO::count() {
    uint32_t counter = 0;
    // read bam file record by record
    while (next())
        counter++;
    return counter;
}

bool BamIO::query(const std::string &queryChrom, int queryStart, int queryEnd) {
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

void BamIO::depth(std::vector<int> &coverage, const std::string &queryChrom, int queryStart, int queryEnd) {

    query(queryChrom, queryStart, queryEnd);

    mplp_aux_t data;
    mplp_aux_t *p_data = &data;
    void *v_data = p_data;
    p_data->fp = m_handleBam;
    p_data->hdr = m_handleHeader;
    p_data->iter = m_handleIterator;
    p_data->min_mapQ = m_mapq;
    p_data->min_len = m_readLength;


    bam_mplp_t mplp = bam_mplp_init(1, read_bam,  &v_data);
    bam_mplp_set_maxcnt(mplp, INT_MAX);

    bam_pileup1_t plp;
    const bam_pileup1_t *p_plp = &plp;
    const bam_pileup1_t **pp_plp = &p_plp;

    int ret, tid, pos, n_plp;
    std::vector<int>::iterator it;

    while ((ret=bam_mplp_auto(mplp, &tid, &pos, &n_plp, pp_plp)) > 0) {

        if ((pos < queryStart) || (queryEnd <= pos)) continue;

        int m = 0;
        for (int j = 0; j < n_plp; ++j) {
            const bam_pileup1_t *p = p_plp + j;
            if (p->is_del || p->is_refskip) ++m;
            else if (p->qpos < p->b->core.l_qseq && bam_get_qual(p->b)[p->qpos] < data.min_mapQ) ++m;
        }


        it = coverage.begin() + pos - queryStart;
        if (it < coverage.end())
            *it += n_plp - m;
    }

    bam_mplp_destroy(mplp);
}


