#ifndef PARSERBED_H
#define PARSERBED_H

#include <iostream>
#include <htslib/hts.h>
#include <htslib/kstring.h>

class ParserBed
{
public:
    ParserBed(const std::string &fileName = "");
    ~ParserBed();

    int next();
    kstring_t buffer;

private:
    htsFile *m_fh;
};

#endif /* PARSERBED_H */