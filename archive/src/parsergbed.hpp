#ifndef PARSERGBED_H
#define PARSERGBED_H

#include <iostream>

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

class ParserGbed
{
public:
    ParserGbed();
    ~ParserGbed();

    void open(const std::string &file);
    int grab(const std::string &query);
    int next();

    kstring_t buffer;

private:
    hts_itr_t *m_iterator;
    htsFile *m_handleFile;
    tbx_t *m_handleIndex;
};

#endif /* PARSERGBED_H */