#include "parserbam.hpp"

ParserBam::ParserBam(const std::string &fileName, int mapq) :
    m_bam(nullptr),
    m_header(nullptr),
    m_bai(nullptr),
    m_iter(nullptr),
    m_alignment(nullptr),
    m_min_mapQ(mapq)
{
    m_bam = hts_open(fileName.c_str(), "r");
    if (!m_bam) {
        std::cerr << "ParserBam::Error, failed to open BAM file " << fileName << std::endl;
        return;
    }

    m_header = sam_hdr_read(m_bam);
    if (!m_header) {
        std::cerr << "ParserBam::Error, failed to read BAM header " << fileName << std::endl;
        return;
    }

    m_bai = hts_idx_load(fileName.c_str(), HTS_FMT_BAI);
    if (!m_bai) {
        std::cerr << "ParserBam::Error, failed to read BAM index " << fileName << ".bai" << std::endl;
        return;
    }

    m_alignment = bam_init1();
    if (!m_alignment) {
        std::cerr << "ParserBam::Error, failed to initialize alignemnt record" << std::endl;
        return;
    }
}

ParserBam::~ParserBam()
{
    if (m_alignment != nullptr)
        bam_destroy1(m_alignment);

    if (m_iter != nullptr)
        hts_itr_destroy(m_iter);

    if (m_bai != nullptr)
        hts_idx_destroy(m_bai);

    if (m_header != nullptr)
        bam_hdr_destroy(m_header);

    if (m_bam != nullptr)
        hts_close(m_bam);
}


void ParserBam::query(const std::string &chrom, int chromStart, int chromEnd)
{
    int chromTid = bam_name2id(m_header, chrom.c_str());

    if (m_iter)
        hts_itr_destroy(m_iter);
    
    m_iter = bam_itr_queryi(m_bai, chromTid, chromStart, chromEnd);
}


int ParserBam::next()
{
    int ret;
    while (1) {
        ret = m_iter ? sam_itr_next(m_bam, m_iter, m_alignment) : sam_read1(m_bam, m_header, m_alignment);
        if (ret < 0) break;
        if (m_alignment->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        if ((int)m_alignment->core.qual < m_min_mapQ) continue;
        break;
    }
    return ret;
}