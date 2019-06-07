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