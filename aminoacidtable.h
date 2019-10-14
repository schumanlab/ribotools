#ifndef AMINOACIDTABLE_H
#define AMINOACIDTABLE_H

#include <iostream>
#include <map>
#include <cmath>
#include <memory>

class AminoAcid {
public:
    explicit AminoAcid(char letter = 'X',
                       std::string codon = "NNN",
                       std::string code = "Xaa",
                       std::string name = "unknown",
                       int count = 0,
                       double value = 0.0,
                       double score = 0.0) :
        letter(letter),
        codon(codon),
        code(code),
        name(name),
        count(count),
        value(value),
        score(score) {}

    char letter;
    std::string codon;
    std::string code;
    std::string name;
    int count;
    double value;
    double score;
};


class AminoAcidTable
{
public:
    explicit AminoAcidTable();

    void reset();

    void setCount(const std::string &codon, int count);
    int count(const std::string &codon) const;

    void setValue(const std::string &codon, double value);
    double value(const std::string &codon) const;

    void setScore(const std::string &codon, double score);
    double score(const std::string &codon) const;

    void addParameters(const std::string &codon, int count = 0, double value = 0, double score = 0);
    void addPausingScore(const std::string &codon, double zscore);
    
    char translate(const std::string &codon) const;
    
    void write();
    void log(const std::string &label);

    void test();

    void calculateRSCU();

private:

    std::map<std::string, std::shared_ptr<AminoAcid>> m_table;
    std::multimap<char, std::shared_ptr<AminoAcid>> m_synonyms;
};

#endif
