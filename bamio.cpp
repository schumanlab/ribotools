#include "bamio.h"

BamAuxiliary::BamAuxiliary(const std::string &fileBam, int mapq, int minlen) :
    m_name(fileBam),
    m_error("")
{
    // open HTS file
    m_aux.fp = hts_open(fileBam.c_str(), "r");
    if (!m_aux.fp) {
        m_error = "BamAuxiliary::error, failed to open BAM file " + fileBam;
        return;
    }

    // load HTS index
    m_aux.idx = hts_idx_load(fileBam.c_str(), HTS_FMT_BAI);
    if (!m_aux.idx) {
        m_error = "BamAuxiliary::error, failed to laod BAI index " + fileBam + ".bai";
        return;
    }

    // load SAM header
    m_aux.hdr = sam_hdr_read(m_aux.fp);
    if (!m_aux.hdr) {
        m_error = "BamAuxiliary::error, failed to parse BAM header";
        return;
    }

    // set empty iterator
    m_aux.iter = nullptr;

    // alocate BAM record
    m_bam = bam_init1();
    if (!m_bam) {
        m_error = "BamAuxiliary::error, failed to allocate alignment memory";
        return;
    }

    // set defaults
    m_aux.min_mapQ = mapq;
    m_aux.min_len = minlen;
}

BamAuxiliary::~BamAuxiliary() {
    if (m_bam)
        bam_destroy1(m_bam);

    if (m_aux.iter)
        bam_itr_destroy(m_aux.iter);

    if (m_aux.hdr)
        bam_hdr_destroy(m_aux.hdr);

    if (m_aux.idx)
        hts_idx_destroy(m_aux.idx);

    if (m_aux.fp)
        hts_close(m_aux.fp);
}

bool BamAuxiliary::isOpen() const { return m_error.empty();}

const std::string BamAuxiliary::name() const {
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

const std::string BamAuxiliary::what() const {return m_error;}

bool BamAuxiliary::next() {return read_bam(&m_aux, m_bam) > 0;}

bool BamAuxiliary::query(const std::string &queryChrom, int queryStart, int queryEnd) {
    int queryTid = bam_name2id(m_aux.hdr, queryChrom.c_str());
    m_aux.iter = bam_itr_queryi(m_aux.idx, queryTid, queryStart, queryEnd);
    return (m_aux.iter != nullptr);
}

int BamAuxiliary::count() {
    int counter = 0;
    while (next())
        counter++;
    return counter;
}

int BamAuxiliary::readStart() const {return m_bam ? m_bam->core.pos : 0;}

int BamAuxiliary::readEnd() const {return m_bam ?
                m_bam->core.pos +
                (m_bam->core.n_cigar ?
                     bam_cigar2rlen(m_bam->core.n_cigar, bam_get_cigar(m_bam)) : 1) : 0;}

int BamAuxiliary::readLength() const {
    if (!m_bam)
        return 0;
    const uint8_t *p_cigar = m_bam->data + m_bam->core.l_qname;
    return bam_cigar2qlen(static_cast<int>(m_bam->core.n_cigar), reinterpret_cast<const uint32_t *>(p_cigar));
}

std::vector<int> BamAuxiliary::depth(const std::string &queryChrom, int queryStart, int queryEnd) {
    std::vector<int> cov;

    // check query range
    if (queryEnd <= queryStart)
        return cov;

    // check if region exists
    if (!query(queryChrom, queryStart, queryEnd))
        return cov;

    // resize vector
    cov.resize(static_cast<std::size_t>(queryEnd - queryStart));

    // pileup algorithm htslib
    bam_aux_t *p_aux = &m_aux;
    void *v_aux = p_aux;
    bam_pileup1_t plp;
    const bam_pileup1_t *p_plp = &plp;
    int ret, tid, pos, n_plp = 0;
    std::vector<int>::iterator ov = cov.begin();

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, &v_aux);
    bam_mplp_set_maxcnt(mplp, INT_MAX);

    while ((ret = bam_mplp_auto(mplp, &tid, &pos, &n_plp, &p_plp)) > 0) {

        // filter position by range
        if ((pos < queryStart) || (queryEnd <= pos)) continue;

        // filter piled reads by alignment quality and read length
        int m = 0;
        for (int j = 0; j < n_plp; ++j) {
            const bam_pileup1_t *p = p_plp + j;
            if (p->is_del || p->is_refskip) ++m;
            else if (p->qpos < p->b->core.l_qseq && bam_get_qual(p->b)[p->qpos] < m_aux.min_mapQ) ++m;
        }

        // update coverage
        ov = cov.begin() + pos - queryStart;
        if (ov < cov.end())
            *ov = n_plp - m;
    }

    bam_mplp_destroy(mplp);

    return cov;
}

int BamAuxiliary::read_bam(void *data, bam1_t *b) {
    bam_aux_t *aux = static_cast<bam_aux_t*>(data);
    int ret;
    while (1) {
        ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( static_cast<int>(b->core.qual) < aux->min_mapQ ) continue;
        if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue; // filter based on read length
        break;
    }
    return ret;
}
