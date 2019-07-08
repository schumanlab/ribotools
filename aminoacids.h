#ifndef AMINOACIDS_H
#define AMINOACIDS_H

#include <iostream>
#include <unordered_map>

class AminoAcids
{
public:
    explicit AminoAcids();

    void add(const std::string &codon, int value);
    void reset();

    void setTimeDecoding(const std::string &codon, double timeValue);
    void setTimePausing(const std::string &codon, double timeValue);
    
    double timeDecoding(const std::string &codon) const;
    double timePausing(const std::string &codon) const;
    
    void write();

private:
    struct Node
    {
        char letter;
        int count;
        double timeDecoding;
        double timePausing;
        std::string codon;
        std::string code;
        std::string name;
    };

    std::unordered_map<std::string, Node> m_map;
};

#endif