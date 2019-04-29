#ifndef COVERAGE_H
#define COVERAGE_H

#include "mex.h"
#include "class_handle.h"
#include "bedreduced.h"
#include <iostream>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

class Coverage
{
public:
    explicit Coverage(char *fileName);
    ~Coverage();

    void query(int32_t *track, int32_t *counts, char *chrom, int32_t chromStart, int32_t chromEnd);

private:
    kstring_t m_buffer;
    hts_itr_t *m_iterator;
    htsFile *m_handleFile;
    tbx_t *m_handleIndex;

};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* COVERAGE_H */