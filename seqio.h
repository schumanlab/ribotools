#ifndef SEQIO_H
#define SEQIO_H

#include <iostream>
#include <htslib/faidx.h>

class SeqIO
{
public:
    explicit SeqIO();
    explicit SeqIO(const std::string &fileName);
    ~SeqIO();
    bool open(const std::string &fileName);
    bool isOpen() const {return m_handleFai != nullptr;}
    const std::string error() const {return m_error;}
    const char *sequence() const {return m_sequence;}
    bool fetch(const std::string &queryName, int queryStart, int queryEnd);

private:
    faidx_t *m_handleFai;
    char *m_sequence;
    std::string m_error;

    void load(const std::string &fileName);
};

#endif // SEQIO_H
