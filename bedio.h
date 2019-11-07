#ifndef BEDIO_H
#define BEDIO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

class Bed12
{
public:
    explicit Bed12(){}
    void parse(const std::string &line);
    bool overlap(const std::string &queryChrom, int queryStart, int queryEnd);
    std::string name(std::size_t index = 0, const char delimiter = ';') const;

    char strand;
    int chromStart;
    int chromEnd;
    int thickStart;
    int thickEnd;
    int orfStart;
    int orfEnd;
    int orfSpan;
    int span;
    int score;
    int blocks;

    std::string chrom;
    std::string label;
    std::string itemRgb;

    std::vector<int> blockSizes;
    std::vector<int> blockStarts;

private:
    void parseBlocks(std::vector<int> &array, const std::string &list);
    void parseExons();

};


class BedIO
{
public:
    explicit BedIO(const std::string &fileName);
    ~BedIO();

    bool next();
    std::string line() const {return m_buffer;}
    Bed12 bed() const {return m_bed;}

private:
    std::ifstream m_fileStream;
    std::string m_buffer;
    Bed12 m_bed;
};

#endif // BEDIO_H
