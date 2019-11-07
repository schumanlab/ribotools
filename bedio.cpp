#include "bedio.h"

BedIO::BedIO(const std::string &fileName)
{
    m_fileStream.open(fileName);
    if (!m_fileStream.is_open())
        throw std::logic_error("ribotools::bedio::error, failed to open BED file " + fileName);
}


BedIO::~BedIO()
{
    if (m_fileStream.is_open())
        m_fileStream.close();
}


bool BedIO::next()
{
    bool isNext = false;
    if (std::getline(m_fileStream, m_buffer)) {
        m_bed.parse(m_buffer);
        isNext = true;
    }

    return isNext;
}


void Bed12::parse(const std::string &line)
{
    std::istringstream iss(line);
    std::string blockSizesList;
    std::string blockStartsList;

    // parse stream
    iss >>
            chrom >>
            chromStart >>
            chromEnd >>
            label >>
            score >>
            strand >>
            thickStart >>
            thickEnd >>
            itemRgb >>
            blocks >>
            blockSizesList >>
            blockStartsList;

    if (!iss.eof() || iss.fail())
        throw std::logic_error("ribotools::bedio::error, failed to convert value\n" + line);

    // parse blocks
    parseBlocks(blockSizes, blockSizesList);
    parseBlocks(blockStarts, blockStartsList);

    if (blockSizes.size() != blockStarts.size())
        throw std::logic_error("ribotools::bedio::error, blocks mismatch\n" + line);

    // parse exons
    parseExons();

}

bool Bed12::overlap(const std::string &queryChrom, int queryStart, int queryEnd)
{
    return (chrom.compare(queryChrom) == 0) && (chromStart <= queryEnd) && (queryStart <= chromEnd);
}

std::string Bed12::name(std::size_t index, const char delimiter) const
{
    if (index == 0)
        return label;

    // split label by delimiter
    std::string _name = label;
    std::istringstream iss(label);
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
    span = 0;
    std::vector<int>::const_iterator itStarts;
    std::vector<int>::const_iterator itSizes;
    for (itStarts = blockStarts.begin(), itSizes = blockSizes.begin();
         itStarts != blockStarts.end() && itSizes != blockSizes.end();
         ++itStarts, ++itSizes) {

        int exonStart = chromStart + *itStarts;
        int exonEnd = exonStart + *itSizes;

        if ((exonStart <= thickStart) && (thickStart <= exonEnd))
            orfStart = (thickStart - exonStart) + span;

        if ((exonStart <= thickEnd) && (thickEnd <= exonEnd))
            orfEnd = (thickEnd - exonStart) + span;

        span += *itSizes;
    }

    // swap orfStart and orfEnd if negative strand
    if (strand == '-') {
        orfStart = span - orfStart;
        orfEnd = span - orfEnd;
        std::swap(orfStart, orfEnd);
    }

    // correct for unsynced ORF
    orfSpan = orfEnd - orfStart;
    orfSpan = orfSpan - (orfSpan % 3);
    orfEnd = orfStart + orfSpan;
}
