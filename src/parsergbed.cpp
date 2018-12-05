#include "parsergbed.hpp"

ParserGbed::ParserGbed(const std::string &file) :
    m_iterator(nullptr),
    m_handleFile(nullptr),
    m_handleIndex(nullptr)
{
    m_handleFile = hts_open(file.c_str(), "r");
    if (!m_handleFile)
    {
        std::cerr << "MetaGene: could not open coverage file " << file << std::endl;
        return;
    }
        
    m_handleIndex = tbx_index_load(file.c_str());
    if (!m_handleIndex)
    {
        std::cerr << "MetaGene: could not open index file " << file << ".bai" << std::endl;
        return;
    }
}


ParserGbed::~ParserGbed()
{
    if (m_iterator != nullptr)
    {
        tbx_itr_destroy(m_iterator);
    }

    if (m_handleFile != nullptr)
    {
        hts_close(m_handleFile);
    }
    
    if (m_handleIndex != nullptr)
    {
        tbx_destroy(m_handleIndex);
    }
}


int ParserGbed::grab(const std::string &query)
{
    m_iterator = tbx_itr_querys(m_handleIndex, query.c_str());
    if (!m_iterator)
    {
        std::cerr << "ParserGbed::Error: failed to query index " << query << std::endl;
        return 1;
    }

    return 0;
}


int ParserGbed::next()
{
    int code = tbx_itr_next(m_handleFile, m_handleIndex, m_iterator, &buffer);
    m_iterator->i++;
    return code;
}