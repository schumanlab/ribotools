#include "bamhandle.h"

BamHandle::BamHandle(const std::string &fileName) :
    name(fileName),
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
