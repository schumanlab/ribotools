#ifndef BEDREDUCED_H
#define BEDREDUCED_H

#include <sstream>
#include <string>

class BedReduced
{
public:
    BedReduced() : chromStart(0), chromEnd(0), chrom("") {}

    uint32_t chromStart;
    uint32_t chromEnd;
    std::string chrom;

    friend std::istream& operator>> (std::istream& in, BedReduced &data)
    {
        BedReduced temp;
        if ((in >> temp.chrom) &&
            (in >> temp.chromStart) &&
            (in >> temp.chromEnd))
        {
            data.swap(temp);
        }
        return in;
    }

private:
    void swap(BedReduced &other)
    {
        std::swap(chrom, other.chrom);
        std::swap(chromStart, other.chromStart);
        std::swap(chromEnd, other.chromEnd);
    }
};

#endif /* BEDREDUCED_H */
