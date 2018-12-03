#ifndef METAGENE_H
#define METAGENE_H

#include <iostream>
#include <sstream>

#include <htslib/bgzf.h>
#include <htslib/regidx.h>
#include <htslib/kstring.h>

#include "argparse.hpp"
#include "bedrecord.hpp"
#include "version.hpp"

class MetaGene
{
public:
    MetaGene();
    ~MetaGene();

    bool parse(int argc, char const *argv[]);
    void pileup();
    
private:
    int m_nbins;
    float m_depth;
    float m_boundLower;
    float m_boundUpper;
    std::string m_fileBed;
    std::vector<std::string> m_filesGbed;

    BGZF *m_fhBed;
    std::vector<BGZF *> m_fhGbed;
    std::vector<regidx_t *> m_fhTabix;

    static int tabixParse(const char *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *usr);
    static void tabixFree(void *payload);


    void open();
    void readBed();
    void queryBed(const BedRecord &bed);

    void help();
    void error(const std::string &errorMessage);
    
};

#endif /* METAGENE_H */