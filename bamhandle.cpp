#include "bamhandle.h"

BamHandle::BamHandle(const std::string &fileName, int mapq, int length) :
    m_mapq(mapq),
    m_length(length),
    m_reads(0),
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
        if ( static_cast<int>(b->core.qual) < m_mapq ) continue;
        if ( m_length && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < m_length ) continue;
        break;
    }
    return ret;
}


void BamHandle::calculateFootprintCoverage(std::vector<int> &fc, const std::string &qName, int qStart, int qEnd)
{
    bam1_t *alignment = bam_init1();
    query(qName, qStart, qEnd);

    auto pileup = [](int &n){n++;};

    while (readBam(alignment) > 0) {

        // read linear coordinates
        int readStart = alignment->core.pos;
        int readLength = bam_cigar2qlen(static_cast<int>(alignment->core.n_cigar), bam_get_cigar(alignment));

        int offsetStart = (qStart < readStart) ? (readStart - qStart) : 0;
        int offsetEnd = ((offsetStart + readLength) < qEnd) ? (offsetStart + readLength) : (qEnd - qStart);
        std::for_each(fc.begin() + offsetStart, fc.begin() + offsetEnd, pileup);

    }
    
    if (alignment)
        bam_destroy1(alignment);
}


void BamHandle::calculateSiteCoverage(std::vector<int> &fc, const std::string &qName, int qStart, int qEnd, int qRef, bool useAsite)
{
    bam1_t *alignment = bam_init1();
    query(qName, qStart, qEnd);

    while (readBam(alignment) > 0) {

        // read linear coordinates
        int readStart = alignment->core.pos;
        int readPsite = (readStart + 12 - qRef);
        int codonIndex = (readPsite - (readPsite % 3)) / 3;

        if (useAsite)
            codonIndex += 1; // increment P-site codon with 1 to get A-site codon

        if ((0 <= codonIndex) && (static_cast<size_t>(codonIndex) < fc.size()))
            fc.at(codonIndex)++;

    }

    if (alignment)
        bam_destroy1(alignment);
}


void BamHandle::countUniqueReads()
{
    bam1_t *alignment = bam_init1();
    m_reads = 0; // reset counter

    // check if iterator is set and reset
    if (m_iterator) {
        bam_itr_destroy(m_iterator);
        m_iterator = nullptr;
    }

    // read bam file record by record
    while (readBam(alignment) > 0)
        m_reads++;

    if (alignment)
        bam_destroy1(alignment);
}


int BamHandle::readsPerRegion(const std::string &qName, int qStart, int qEnd)
{
    bam1_t *alignment = bam_init1();
    int reads = 0;
    query(qName, qStart, qEnd);

    // read alignments from iterator
    while (readBam(alignment) > 0)
        reads++;


    if (alignment)
        bam_destroy1(alignment);

    return reads;
}


void BamHandle::calculateGCcontent(double &gc_mean, double &gc_M2, double &gc_var, int &readCount)
{
    bam1_t *alignment = bam_init1();

    // read alignments from iterator
    while (readBam(alignment) > 0) {

        int gc = 0;
        int readLength = alignment->core.l_qseq;
        uint8_t *sequence = bam_get_seq(alignment);

        for (int b = 0; b < readLength; ++b) {
            char base = "=ACMGRSVTWYHKDBN"[bam_seqi(sequence, b)];
            if ((base == 'C') || (base == 'G'))
                    gc++;
        }

        double gc_ratio = static_cast<double>(gc) / readLength;
        double gc_delta = gc_ratio - gc_mean;
        gc_mean += gc_delta / (readCount + 1);
        gc_M2 += gc_delta * (gc_ratio - gc_mean);
        gc_var = gc_M2 / (readCount + 1);
        readCount++;
        std::cout << readLength << "\t" << gc_ratio << std::endl;
    }

    if (alignment)
        bam_destroy1(alignment);
}


void BamHandle::calculateGCbases(uint64_t &baseCount, uint64_t &gcCount)
{
    bam1_t *alignment = bam_init1();

    // read alignments from iterator
    while (readBam(alignment) > 0) {

        int gc = 0;
        int readLength = alignment->core.l_qseq;
        uint8_t *sequence = bam_get_seq(alignment);

        for (int b = 0; b < readLength; ++b) {
            char base = "=ACMGRSVTWYHKDBN"[bam_seqi(sequence, b)];
            if ((base == 'C') || (base == 'G'))
                    gc++;
        }

        baseCount += static_cast<uint64_t>(readLength);
        gcCount += static_cast<uint64_t>(gc);

    }

    if (alignment)
        bam_destroy1(alignment);
}


void BamHandle::calculateBaseContent(int bufferLength, std::vector<int> &base_A, std::vector<int> &base_C, std::vector<int> &base_G, std::vector<int> &base_T, std::vector<int> &base_N)
{
    bam1_t *alignment = bam_init1();

    // check if iterator is set and reset
    if (m_iterator) {
        bam_itr_destroy(m_iterator);
        m_iterator = nullptr;
    }


    // read alignments from iterator
    while (readBam(alignment) > 0) {

        int readLength = alignment->core.l_qseq;
        uint8_t *sequence = bam_get_seq(alignment);
        int loopLength = std::min(bufferLength, readLength);

        for (int b = 0; b < loopLength; ++b) {
            char base = "=ACMGRSVTWYHKDBN"[bam_seqi(sequence, b)];
            if (base == 'A')
                base_A.at(b)++;
            else if (base == 'C')
                base_C.at(b)++;
            else if (base == 'G')
                base_G.at(b)++;
            else if (base == 'T')
                base_T.at(b)++;
            else
                base_N.at(b)++;
        }


    }

    if (alignment)
        bam_destroy1(alignment);
}
