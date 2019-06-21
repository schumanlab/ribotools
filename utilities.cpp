#include "utilities.h"

double calculateAverageFootprintCoverage(const std::vector<int> &fc, int qStart, int qEnd)
{
    qStart = std::max(0, qStart);
    qEnd = std::min(static_cast<int>(fc.size()), qEnd);
    if (qEnd - qStart < 1) return 0.0;

    int sum = 0;
    int norm = 0;
    for (int k = qStart; k < qEnd; k++) {
        sum += fc[k];
        if (fc[k] > 0) norm++;
    }
        
    
    return static_cast<double>(sum) / norm;
}


void calculateFootprintCoverage(std::vector<int> &fc, BamHandle *handle, const std::string &qName, int qStart, int qEnd)
{
    bam1_t *alignment = bam_init1();
    handle->query(qName, qStart, qEnd);

    while (handle->readBam(alignment) > 0) {

        int readStart = alignment->core.pos;
        int readLength = bam_cigar2qlen(alignment->core.n_cigar, bam_get_cigar(alignment));

        // accumulate P-site per read
        int readPsite = readStart + readLength - 17 - qStart;
        if ((0<= readPsite) && (readPsite < fc.size()))
            fc[readPsite]++;
    }
    
    if (alignment)
        bam_destroy1(alignment);

}