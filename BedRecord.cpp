#include "bedrecord.h"

BedRecord::BedRecord() :
    strand(0),
    chromStart(0),
    chromEnd(0),
    thickStart(0),
    thickEnd(0),
    score(0),
    blocks(0),
    cdsStart(0),
    cdsEnd(0),
    span(0),
    chrom(""),
    name(""),
    itemRgb("")
{

}


void BedRecord::swap(BedRecord &other)
{
    std::swap(chrom, other.chrom);
    std::swap(chromStart, other.chromStart);
    std::swap(chromEnd, other.chromEnd);
    std::swap(name, other.name);
    std::swap(score, other.score);
    std::swap(strand, other.strand);
    std::swap(thickStart, other.thickStart);
    std::swap(thickEnd, other.thickEnd);
    std::swap(itemRgb, other.itemRgb);
    std::swap(blocks, other.blocks);
    std::swap(blockSizes, other.blockSizes);
    std::swap(blockStarts, other.blockStarts);
}

void BedRecord::listToArray(std::vector<int32_t> &array, const std::string &list)
{
    std::stringstream ss(list);
    int32_t value;
    while (ss >> value)
    {
        array.push_back(value);
        if (ss.peek() == ',' || ss.peek() == ' ')
            ss.ignore();
    }
}


std::istream& operator>> (std::istream& in, BedRecord &data)
{
    BedRecord temp;
    std::string tempBlockSizes;
    std::string tempBlockStarts;
    if ((in >> temp.chrom) &&
        (in >> temp.chromStart) &&
        (in >> temp.chromEnd) &&
        (in >> temp.name) &&
        (in >> temp.score) &&
        (in >> temp.strand) &&
        (in >> temp.thickStart) &&
        (in >> temp.thickEnd) &&
        (in >> temp.itemRgb) &&
        (in >> temp.blocks) &&
        (in >> tempBlockSizes) &&
        (in >> tempBlockStarts))
    {
        BedRecord::listToArray(temp.blockSizes, tempBlockSizes);
        BedRecord::listToArray(temp.blockStarts, tempBlockStarts);
        data.swap(temp);
    }
    return in;
}


void BedRecord::parseExons()
{
    for (int32_t k = 0; k < blocks; k++)
    {
        int32_t exonStart = chromStart + blockStarts[k];
        int32_t exonEnd = exonStart + blockSizes[k];
        ExonNode exon = {exonStart, exonEnd, span};
        m_exonTree.insert(exon);

        if ((exonStart <= thickStart) && (thickStart <= exonEnd))
        {
            cdsStart = (thickStart - exonStart) + span;
        }
        
        if ((exonStart <= thickEnd) && (thickEnd <= exonEnd))
        {
            cdsEnd = (thickEnd - exonStart) + span;
        }

        span += blockSizes[k];
    }

    cdsSpan = cdsEnd - cdsStart;
}


bool BedRecord::toLinear(int32_t &readOffset, int32_t readStart)
{
    bool found = false;
    ExonNode query = {readStart, readStart, 0};
    auto search = m_exonTree.find(query);

    if (search != m_exonTree.end())
    {
        readOffset = (readStart - search->exonStart) + search->offset;
        found = true;
    }

    return found;
}