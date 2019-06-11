#ifndef BEDRECORD_H
#define BEDRECORD_H

#include <sstream>
#include <string>
#include <vector>

class BedRecord
{
public:
    BedRecord();

    char strand;
    int chromStart;
    int chromEnd;
    int thickStart;
    int thickEnd;
    int score;
    int blocks;
    int cdsStart;
    int cdsEnd;
    int cdsSpan;
    int span;
    std::string chrom;
    std::string name;
    std::string itemRgb;
    std::string transcript;
    std::string gene;
    std::vector<int> blockSizes;
    std::vector<int> blockStarts;

    static void listToArray(std::vector<int> &array, const std::string &list);
    friend std::istream& operator>> (std::istream& in, BedRecord &data);
    void parseExons();
    bool overlap(const std::string &readChrom, int readStart, int readEnd);

    int nextCodon(int position);
    int prevCodon(int position);

private:
    void swap(BedRecord &other);
    static int closestNumber(int n, int m);
};

#endif /* BEDLINE_H */