#include "bedio.h"

bool Bed12::parse(const std::string &line)
{
    std::istringstream iss(line);
    std::string blockSizesList;
    std::string blockStartsList;

    // parse stream
    iss >>
            m_chrom >>
            m_chromStart >>
            m_chromEnd >>
            m_name >>
            m_score >>
            m_strand >>
            m_thickStart >>
            m_thickEnd >>
            m_itemRgb >>
            m_blocks >>
            blockSizesList >>
            blockStartsList;

    if (!iss.eof() || iss.fail())
        return false;

    // parse blocks
    parseBlocks(m_blockSizes, blockSizesList);
    parseBlocks(m_blockStarts, blockStartsList);

    // parse exons
    parseExons();

    return true;
}

bool Bed12::overlap(const std::string &queryChrom, int queryStart, int queryEnd)
{
    return (m_chrom.compare(queryChrom) == 0) && (m_chromStart <= queryEnd) && (queryStart <= m_chromEnd);
}

const std::string Bed12::name(std::size_t index, const char delimiter) const
{
    if (index == 0)
        return m_name;

    // split label by delimiter
    std::string _name = m_name;
    std::istringstream iss(m_name);
    std::string token;
    std::vector<std::string> listOfTokens;
    while (std::getline(iss, token, delimiter))
        listOfTokens.push_back(token);

    index--;
    if (index < listOfTokens.size())
        _name = listOfTokens.at(index);

    return _name;
}

void Bed12::parseBlocks(std::vector<int> &array, const std::string &list)
{
    std::stringstream iss(list);
    array = {};
    int value = 0;

    while (iss >> value) {
        array.push_back(value);
        if (iss.peek() == ',' || iss.peek() == ' ')
            iss.ignore();
    }
}

void Bed12::parseExons()
{
    m_span = 0;
    std::vector<int>::const_iterator itStarts;
    std::vector<int>::const_iterator itSizes;
    for (itStarts = m_blockStarts.begin(), itSizes = m_blockSizes.begin();
         itStarts != m_blockStarts.end() && itSizes != m_blockSizes.end();
         ++itStarts, ++itSizes) {

        int exonStart = m_chromStart + *itStarts;
        int exonEnd = exonStart + *itSizes;

        if ((exonStart <= m_thickStart) && (m_thickStart <= exonEnd))
            m_orfStart = (m_thickStart - exonStart) + m_span;

        if ((exonStart <= m_thickEnd) && (m_thickEnd <= exonEnd))
            m_orfEnd = (m_thickEnd - exonStart) + m_span;

        m_span += *itSizes;
    }

    // swap orfStart and orfEnd if negative strand
    if (m_strand == '-') {
        m_orfStart = m_span - m_orfStart;
        m_orfEnd = m_span - m_orfEnd;
        std::swap(m_orfStart, m_orfEnd);
    }

    // correct for unsynced ORF
    m_orfSpan = m_orfEnd - m_orfStart;
    m_orfSpan = m_orfSpan - (m_orfSpan % 3);
    m_orfEnd = m_orfStart + m_orfSpan;
}

BedIO::BedIO() {}

BedIO::BedIO(const std::string &fileName) {
    m_fileStream.open(fileName);
    if (!m_fileStream.is_open())
        m_error = "bedio::error, failed to open " + fileName;
}

BedIO::~BedIO() {
    if (m_fileStream.is_open())
        m_fileStream.close();
}

bool BedIO::open(const std::string &fileName) {
    m_fileStream.open(fileName);
    if (!m_fileStream.is_open()) {
        m_error = "bedio::error, failed to open " + fileName;
        return false;
    }
    return true;
}

bool BedIO::next() {

    bool isNext = false;

    if (std::getline(m_fileStream, m_buffer)) {

        if (!m_bed.parse(m_buffer))
            m_error = "bedio::error, failed to parse line " + m_buffer;
        else
            isNext = true;
    }

    if (m_fileStream.bad()) {
        m_error = "bedio::error, failed while reading line";
        isNext = false;
    }

    return isNext;
}

void BedIO::rewind() {
    m_fileStream.clear();
    m_fileStream.seekg(0);
}
