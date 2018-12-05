#include "parserbed.hpp"

ParserBed::ParserBed(const std::string &fileName) :
    m_fh(nullptr)
{
    if (!fileName.empty())
        m_fh = hts_open(fileName.c_str(), "r");

}


ParserBed::~ParserBed()
{
    if (m_fh != nullptr)
    {
        free(ks_release(&buffer));
        hts_close(m_fh);
    }
}


int ParserBed::next()
{
    return hts_getline(m_fh, '\n', &buffer);
}

