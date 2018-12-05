#ifndef METAGENE_H
#define METAGENE_H

#include <iostream>
#include <vector>

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

#include "bedrecord.hpp"

class MetaGene
{
public:
    MetaGene(const std::string &file);
    ~MetaGene();

    void parse(const BedRecord &bed);

private:
    hts_itr_t *m_iterator;
    htsFile * m_handleFile;
    tbx_t * m_handleIndex;

};

#endif /* METAGENE_H */