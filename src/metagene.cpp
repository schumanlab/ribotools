#include "metagene.hpp"

MetaGene::MetaGene(const std::string &file) :
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


MetaGene::~MetaGene()
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


void MetaGene::parse(const BedRecord &bed)
{
    return;
}

