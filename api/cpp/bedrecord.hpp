#ifndef BEDRECORD_H
#define BEDRECORD_H

#include <sstream>
#include <string>
#include <vector>

class BedRecord
{
public:
    BedRecord();

    uint8_t strand;
    uint32_t chromStart;
    uint32_t chromEnd;
    uint32_t thickStart;
    uint32_t thickEnd;
    uint32_t score;
    uint32_t blocks;
    uint32_t cdsStart;
    uint32_t cdsEnd;
    uint32_t cdsSpan;
    uint32_t span;
    std::string chrom;
    std::string name;
    std::string itemRgb;
    std::vector<uint32_t> blockSizes;
    std::vector<uint32_t> blockStarts;

    static void listToArray(std::vector<uint32_t> &array, const std::string &list);
    friend std::istream& operator>> (std::istream& in, BedRecord &data);
    void parseExons();
    bool overlap(const std::string &readChrom, uint32_t readStart, uint32_t readEnd);
    std::string transcript();


private:
    void swap(BedRecord &other);
};

#endif /* BEDLINE_H */