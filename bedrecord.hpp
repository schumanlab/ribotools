#ifndef BEDRECORD_H
#define BEDRECORD_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>

class BedRecord
{
public:
    BedRecord();

    uint8_t strand;
    uint32_t chromStart;
    uint32_t chromEnd;
    uint32_t thickStart;
    uint32_t thickEnd;
    int32_t score;
    int32_t blocks;
    int32_t cdsStart;
    int32_t cdsEnd;
    int32_t cdsSpan;
    int32_t span;
    std::string chrom;
    std::string name;
    std::string itemRgb;
    std::vector<uint32_t> blockSizes;
    std::vector<uint32_t> blockStarts;

    static void listToArray(std::vector<uint32_t> &array, const std::string &list);
    friend std::istream& operator>> (std::istream& in, BedRecord &data);
    void parseExons();
    bool toLinear(int32_t &readOffset, int32_t readStart);

private:
    struct ExonNode
    {
        int32_t exonStart;
        int32_t exonEnd;
        int32_t offset;
    };
    
    struct ExonCompare
    {
        // overlapping ranges are considered equivalent
        bool operator()(const ExonNode& lhv, const ExonNode& rhv) const
        {
            return lhv.exonStart < rhv.exonStart && lhv.exonEnd < rhv.exonStart;
        }
    };

    std::set<ExonNode, ExonCompare> m_exonTree;
    
    void swap(BedRecord &other);
};

#endif /* BEDLINE_H */