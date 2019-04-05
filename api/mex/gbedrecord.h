#ifndef GBEDRECORD_H
#define GBEDRECORD_H

#include <sstream>
#include <string>


class GbedRecord
{
public:
    GbedRecord();

    uint32_t chromStart;
    uint32_t chromEnd;
    int32_t depth;
    std::string chrom;

    friend std::istream& operator>> (std::istream& in, GbedRecord &data);

private:
    void swap(GbedRecord &other);
};

#endif /* GBEDRECORD_H */
