#include "bamhandle.h"

BamHandle::BamHandle(const std::string &fileName) :
    bam(nullptr),
    bai(nullptr),
    header(nullptr)
{
    // open BAM file
    bam = hts_open(fileName.c_str(), "r");
    if (!bam) {
        std::cerr << "ribotools::bamhandle::error, failed to open BAM file " << fileName << std::endl;
        return;
    }

    // load BAI index
    bai = hts_idx_load(fileName.c_str(), HTS_FMT_BAI);
    if (!bai) {
        std::cerr << "ribotools::bamhandle::error, failed to load BAI index " << fileName << std::endl;
        return;
    }

    // read BAM header
    header = sam_hdr_read(bam);
    if (!header) {
        std::cerr << "ribotools::bamhandle::error, failed to read BAM header " << fileName << std::endl;
        return;
    }
}


BamHandle::~BamHandle()
{
    // destroy header
    if (header)
        bam_hdr_destroy(header);
    
    // destory index
    if (bai)
        hts_idx_destroy(bai);

    // close BAM file
    if (bam)
        hts_close(bam);
}


int BamHandle::orfDepth(std::vector<double> &cov, const std::string &chrom, int chromStart, int chromEnd, int orfStart)
{
    int count = 0;
    int chromTid = bam_name2id(header, chrom.c_str());
    int chromSpan = (chromEnd - chromStart) / 3;
    hts_itr_t *iterator = bam_itr_queryi(bai, chromTid, chromStart, chromEnd);
    bam1_t *alignment = bam_init1();
    cov.resize(chromSpan, 0.0);
    int ret = 0;
    while ((ret = sam_itr_next(bam, iterator, alignment)) >= 0) {

        int readStart = alignment->core.pos;
        int readLength = bam_cigar2qlen(alignment->core.n_cigar, bam_get_cigar(alignment));

        // calculate A-site coverage
        int index = readStart + readLength/2 - orfStart;
        index = (index + 3) - (index % 3);
        index = index / 3;
        if ((0<= index) && (index < chromSpan)) {
            cov[index] += 1.0;
            count++;
        }

    }

    if (alignment)
        bam_destroy1(alignment);
    
    if (iterator)
        bam_itr_destroy(iterator);

    return count;
}