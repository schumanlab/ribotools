#ifndef METAGENE_H
#define METAGENE_H

#include <iostream>
#include <sstream>
#include <map>

#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/khash.h>

#include "argparse.hpp"
#include "bedrecord.hpp"
#include "version.hpp"

class MetaGene
{
public:
    MetaGene();
    ~MetaGene();

    bool parse(int argc, char const *argv[]);
    
private:
    int m_nbins;
    float m_depth;
    float m_boundLower;
    float m_boundUpper;
    std::string m_fileBed;
    std::vector<std::string> m_filesGbed;

    void help();
    void error(const std::string &errorMessage);
    
};

#endif /* METAGENE_H */