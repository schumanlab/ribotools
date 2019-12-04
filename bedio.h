#ifndef BEDIO_H
#define BEDIO_H

#include <fstream>
#include <sstream>
#include <vector>

class Bed12
{
public:
    explicit Bed12() :
        m_strand('+'),
        m_chromStart(0),
        m_chromEnd(0),
        m_thickStart(0),
        m_thickEnd(0),
        m_orfStart(0),
        m_orfEnd(0),
        m_orfSpan(0),
        m_span(0),
        m_score(0),
        m_blocks(0),
        m_chrom(""),
        m_name(""),
        m_itemRgb("")
    {}
    bool parse(const std::string &line);
    bool overlap(const std::string &queryChrom, int queryStart, int queryEnd);

    char strand() const {return m_strand;}
    int chromStart() const {return m_chromStart;}
    int chromEnd() const {return m_chromEnd;}
    int score() const {return m_score;}
    int thickStart() const {return m_thickStart;}
    int thickEnd() const {return m_thickEnd;}
    int blocks() const {return m_blocks;}
    int span() const {return m_span;}
    int orfSpan() const {return m_orfSpan;}
    int orfStart() const {return m_orfStart;}
    int orfEnd() const {return m_orfEnd;}

    const std::string chrom() const {return m_chrom;}
    const std::string name(std::size_t index = 0, const char delimiter = ';') const;
    const std::string itemRgb() const {return m_itemRgb;}

private:   
    char m_strand;
    int m_chromStart;
    int m_chromEnd;
    int m_thickStart;
    int m_thickEnd;
    int m_orfStart;
    int m_orfEnd;
    int m_orfSpan;
    int m_span;
    int m_score;
    int m_blocks;

    std::string m_chrom;
    std::string m_name;
    std::string m_itemRgb;

    std::vector<int> m_blockSizes;
    std::vector<int> m_blockStarts;

    void parseBlocks(std::vector<int> &array, const std::string &list);
    void parseExons();
};


class BedIO
{
public:
    explicit BedIO();
    explicit BedIO(const std::string &fileName);
    ~BedIO();

    bool open(const std::string &fileName);

    bool isOpen() const {return m_fileStream.is_open();}
    const std::string line() const {return m_buffer;}
    const std::string error() const {return m_error;}
    Bed12 bed() const {return m_bed;}
    bool next();
    void rewind();

private:
    std::ifstream m_fileStream;
    std::string m_buffer;
    std::string m_error;
    Bed12 m_bed;
};

#endif // BEDIO_H
