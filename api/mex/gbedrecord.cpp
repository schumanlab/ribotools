#include "gbedrecord.h"

GbedRecord::GbedRecord() :
    chromStart(0),
    chromEnd(0),
    depth(0),
    chrom("")
{

}


void GbedRecord::swap(GbedRecord &other)
{
    std::swap(chrom, other.chrom);
    std::swap(chromStart, other.chromStart);
    std::swap(chromEnd, other.chromEnd);
    std::swap(depth, other.depth);
}


std::istream& operator>> (std::istream& in, GbedRecord &data)
{
    GbedRecord temp;
    if ((in >> temp.chrom) &&
        (in >> temp.chromStart) &&
        (in >> temp.chromEnd) &&
        (in >> temp.depth))
    {
        data.swap(temp);
    }
    return in;
}
