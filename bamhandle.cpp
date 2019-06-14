#include "bamhandle.h"

BamHandle::BamHandle(const std::string &fileName, int mapq, int length) :
    m_mapq(mapq),
    m_length(length),
    m_file(fileName),
    m_bam(nullptr),
    m_bai(nullptr),
    m_header(nullptr),
    m_iterator(nullptr)
{
    // open BAM file
    m_bam = hts_open(fileName.c_str(), "r");
    if (!m_bam) {
        std::cerr << "ribotools::bamhandle::error, failed to open BAM file " << m_file << std::endl;
        return;
    }

    // load BAI index
    m_bai = hts_idx_load(fileName.c_str(), HTS_FMT_BAI);
    if (!m_bai) {
        std::cerr << "ribotools::bamhandle::error, failed to load BAI index " << m_file << std::endl;
        return;
    }

    // read BAM header
    m_header = sam_hdr_read(m_bam);
    if (!m_header) {
        std::cerr << "ribotools::bamhandle::error, failed to read BAM header " << m_file << std::endl;
        return;
    }

}


BamHandle::~BamHandle()
{
    // destroy iterator
    if (m_iterator)
        bam_itr_destroy(m_iterator);

    // destroy header
    if (m_header)
        bam_hdr_destroy(m_header);
    
    // destory index
    if (m_bai)
        hts_idx_destroy(m_bai);

    // close BAM file
    if (m_bam)
        hts_close(m_bam);
}


std::string BamHandle::name()
{
    std::string tag = m_file;

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


void BamHandle::query(const std::string &queryChrom, int queryStart, int queryEnd)
{
    int queryTid = bam_name2id(m_header, queryChrom.c_str());
    m_iterator = bam_itr_queryi(m_bai, queryTid, queryStart, queryEnd);
}


int BamHandle::readBam(bam1_t *b)
{
    int ret;
    while (1) {
        ret = m_iterator ? sam_itr_next(m_bam, m_iterator, b) : sam_read1(m_bam, m_header, b);
        if ( ret < 0) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( (int)b->core.qual < m_mapq ) continue;
        if ( m_length && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < m_length ) continue;
        break;
    }
    return ret;
}


/*
void BamHandle::codonDepth(std::map<int, int> &depth, const std::string &name, int geneSpan, int cdsStart)
{
    int chromTid = bam_name2id(header, name.c_str());

    hts_itr_t *iterator = bam_itr_queryi(bai, chromTid, 0, geneSpan);
    bam1_t *alignment = bam_init1();
    int ret = 0;
    while ((ret = sam_itr_next(bam, iterator, alignment)) >= 0) {

        int readStart = alignment->core.pos;
        int readLength = bam_cigar2qlen(alignment->core.n_cigar, bam_get_cigar(alignment));

        // calculate P-site per read
        int readPsite = readStart + readLength - 17 - cdsStart;
        
        //std::cout << readStart << "\t" << readEnd << "\t" << readLength << "\t" << readStart - cdsStart << "\t" << readEnd - cdsStart << std::endl;
        if (readPsite < 0)
            readPsite -= 2;
        
        readPsite = (readPsite - (readPsite % 3)) / 3;
        depth[readPsite]++;
         
        
    }

    if (alignment)
        bam_destroy1(alignment);
    
    if (iterator)
        bam_itr_destroy(iterator);

    //return count;
}
 */