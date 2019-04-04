#ifndef COVERAGE_H
#define COVERAGE_H

#include "mex.h"
#include "class_handle.h"
#include <iostream>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

class Coverage
{
public:
    explicit Coverage(char *fileName);
    ~Coverage();

    void query(char *chrom, int32_t chromStart, int32_t chromEnd);

private:
    kstring_t m_buffer;
    hts_itr_t *m_iterator;
    htsFile *m_handleFile;
    tbx_t *m_handleIndex;

};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* COVERAGE_H */