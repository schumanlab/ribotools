#include "seqio.h"

SeqIO::SeqIO() :
    m_handleFai(nullptr),
    m_sequence(nullptr),
    m_error("")
{}

SeqIO::SeqIO(const std::string &fileName) :
    m_handleFai(nullptr),
    m_sequence(nullptr),
    m_error("")
{
    load(fileName);
}

SeqIO::~SeqIO() {

    if (m_sequence)
        free(m_sequence);

    if (m_handleFai)
        fai_destroy(m_handleFai);
}

bool SeqIO::open(const std::string &fileName) {
    load(fileName);
    if (isOpen())
        return true;
    return false;
}

bool SeqIO::fetch(const std::string &queryName, int queryStart, int queryEnd) {
    if (m_sequence)
        free(m_sequence);
    int querySpan = queryEnd - queryStart;
    m_sequence = faidx_fetch_seq(m_handleFai, queryName.c_str(), queryStart, queryEnd, &querySpan);
    if (!m_sequence)
        return false;
    return true;
}

void SeqIO::load(const std::string &fileName) {
    m_handleFai = fai_load(fileName.c_str());
    if (!m_handleFai) {
        m_error = "seqio::error, failed to load FASTA file " + fileName;
    }
}
