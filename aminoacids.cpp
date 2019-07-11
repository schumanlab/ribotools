#include "aminoacids.h"

AminoAcids::AminoAcids()
{
    Node aa;

    aa = {'*', 0, 0.0, 0.0, "TAA", "Ter", "Stop_Codon"};
    m_map["TAA"] = aa;
    aa = {'*', 0, 0.0, 0.0, "TAG", "Ter", "Stop_Codon"};
    m_map["TAG"] = aa;
    aa = {'*', 0, 0.0, 0.0, "TGA", "Ter", "Stop_Codon"};
    m_map["TGA"] = aa;
    aa = {'A', 0, 0.0, 0.0, "GCA", "Ala", "Alanine"};
    m_map["GCA"] = aa;
    aa = {'A', 0, 0.0, 0.0, "GCC", "Ala", "Alanine"};
    m_map["GCC"] = aa;
    aa = {'A', 0, 0.0, 0.0, "GCG", "Ala", "Alanine"};
    m_map["GCG"] = aa;
    aa = {'A', 0, 0.0, 0.0, "GCT", "Ala", "Alanine"};
    m_map["GCT"] = aa;
    aa = {'C', 0, 0.0, 0.0, "TGC", "Cys", "Cysteine"};
    m_map["TGC"] = aa;
    aa = {'C', 0, 0.0, 0.0, "TGT", "Cys", "Cysteine"};
    m_map["TGT"] = aa;
    aa = {'D', 0, 0.0, 0.0, "GAC", "Asp", "Aspartic_Acid"};
    m_map["GAC"] = aa;
    aa = {'D', 0, 0.0, 0.0, "GAT", "Asp", "Aspartic_Acid"};
    m_map["GAT"] = aa;
    aa = {'E', 0, 0.0, 0.0, "GAA", "Glu", "Glutamic_Acid"};
    m_map["GAA"] = aa;
    aa = {'E', 0, 0.0, 0.0, "GAG", "Glu", "Glutamic_Acid"};
    m_map["GAG"] = aa;
    aa = {'F', 0, 0.0, 0.0, "TTC", "Phe", "Phenylalanine"};
    m_map["TTC"] = aa;
    aa = {'F', 0, 0.0, 0.0, "TTT", "Phe", "Phenylalanine"};
    m_map["TTT"] = aa;
    aa = {'G', 0, 0.0, 0.0, "GGA", "Gly", "Glycine"};
    m_map["GGA"] = aa;
    aa = {'G', 0, 0.0, 0.0, "GGC", "Gly", "Glycine"};
    m_map["GGC"] = aa;
    aa = {'G', 0, 0.0, 0.0, "GGG", "Gly", "Glycine"};
    m_map["GGG"] = aa;
    aa = {'G', 0, 0.0, 0.0, "GGT", "Gly", "Glycine"};
    m_map["GGT"] = aa;
    aa = {'H', 0, 0.0, 0.0, "CAC", "His", "Histidine"};
    m_map["CAC"] = aa;
    aa = {'H', 0, 0.0, 0.0, "CAT", "His", "Histidine"};
    m_map["CAT"] = aa;
    aa = {'I', 0, 0.0, 0.0, "ATA", "Ile", "Isoleucine"};
    m_map["ATA"] = aa;
    aa = {'I', 0, 0.0, 0.0, "ATC", "Ile", "Isoleucine"};
    m_map["ATC"] = aa;
    aa = {'I', 0, 0.0, 0.0, "ATT", "Ile", "Isoleucine"};
    m_map["ATT"] = aa;
    aa = {'K', 0, 0.0, 0.0, "AAA", "Lys", "Lysine"};
    m_map["AAA"] = aa;
    aa = {'K', 0, 0.0, 0.0, "AAG", "Lys", "Lysine"};
    m_map["AAG"] = aa;
    aa = {'L', 0, 0.0, 0.0, "CTA", "Leu", "Leucine"};
    m_map["CTA"] = aa;
    aa = {'L', 0, 0.0, 0.0, "CTC", "Leu", "Leucine"};
    m_map["CTC"] = aa;
    aa = {'L', 0, 0.0, 0.0, "CTG", "Leu", "Leucine"};
    m_map["CTG"] = aa;
    aa = {'L', 0, 0.0, 0.0, "CTT", "Leu", "Leucine"};
    m_map["CTT"] = aa;
    aa = {'L', 0, 0.0, 0.0, "TTA", "Leu", "Leucine"};
    m_map["TTA"] = aa;
    aa = {'L', 0, 0.0, 0.0, "TTG", "Leu", "Leucine"};
    m_map["TTG"] = aa;
    aa = {'M', 0, 0.0, 0.0, "ATG", "Met", "Methionine"};
    m_map["ATG"] = aa;
    aa = {'N', 0, 0.0, 0.0, "AAC", "Asn", "Asparagine"};
    m_map["AAC"] = aa;
    aa = {'N', 0, 0.0, 0.0, "AAT", "Asn", "Asparagine"};
    m_map["AAT"] = aa;
    aa = {'P', 0, 0.0, 0.0, "CCA", "Pro", "Proline"};
    m_map["CCA"] = aa;
    aa = {'P', 0, 0.0, 0.0, "CCC", "Pro", "Proline"};
    m_map["CCC"] = aa;
    aa = {'P', 0, 0.0, 0.0, "CCG", "Pro", "Proline"};
    m_map["CCG"] = aa;
    aa = {'P', 0, 0.0, 0.0, "CCT", "Pro", "Proline"};
    m_map["CCT"] = aa;
    aa = {'Q', 0, 0.0, 0.0, "CAA", "Gln", "Glutamine"};
    m_map["CAA"] = aa;
    aa = {'Q', 0, 0.0, 0.0, "CAG", "Gln", "Glutamine"};
    m_map["CAG"] = aa;
    aa = {'R', 0, 0.0, 0.0, "AGA", "Arg", "Arginine"};
    m_map["AGA"] = aa;
    aa = {'R', 0, 0.0, 0.0, "AGG", "Arg", "Arginine"};
    m_map["AGG"] = aa;
    aa = {'R', 0, 0.0, 0.0, "CGA", "Arg", "Arginine"};
    m_map["CGA"] = aa;
    aa = {'R', 0, 0.0, 0.0, "CGC", "Arg", "Arginine"};
    m_map["CGC"] = aa;
    aa = {'R', 0, 0.0, 0.0, "CGG", "Arg", "Arginine"};
    m_map["CGG"] = aa;
    aa = {'R', 0, 0.0, 0.0, "CGT", "Arg", "Arginine"};
    m_map["CGT"] = aa;
    aa = {'S', 0, 0.0, 0.0, "AGC", "Ser", "Serine"};
    m_map["AGC"] = aa;
    aa = {'S', 0, 0.0, 0.0, "AGT", "Ser", "Serine"};
    m_map["AGT"] = aa;
    aa = {'S', 0, 0.0, 0.0, "TCA", "Ser", "Serine"};
    m_map["TCA"] = aa;
    aa = {'S', 0, 0.0, 0.0, "TCC", "Ser", "Serine"};
    m_map["TCC"] = aa;
    aa = {'S', 0, 0.0, 0.0, "TCG", "Ser", "Serine"};
    m_map["TCG"] = aa;
    aa = {'S', 0, 0.0, 0.0, "TCT", "Ser", "Serine"};
    m_map["TCT"] = aa;
    aa = {'T', 0, 0.0, 0.0, "ACA", "Thr", "Threonine"};
    m_map["ACA"] = aa;
    aa = {'T', 0, 0.0, 0.0, "ACC", "Thr", "Threonine"};
    m_map["ACC"] = aa;
    aa = {'T', 0, 0.0, 0.0, "ACG", "Thr", "Threonine"};
    m_map["ACG"] = aa;
    aa = {'T', 0, 0.0, 0.0, "ACT", "Thr", "Threonine"};
    m_map["ACT"] = aa;
    aa = {'V', 0, 0.0, 0.0, "GTA", "Val", "Valine"};
    m_map["GTA"] = aa;
    aa = {'V', 0, 0.0, 0.0, "GTC", "Val", "Valine"};
    m_map["GTC"] = aa;
    aa = {'V', 0, 0.0, 0.0, "GTG", "Val", "Valine"};
    m_map["GTG"] = aa;
    aa = {'V', 0, 0.0, 0.0, "GTT", "Val", "Valine"};
    m_map["GTT"] = aa;
    aa = {'W', 0, 0.0, 0.0, "TGG", "Trp", "Tryptophan"};
    m_map["TGG"] = aa;
    aa = {'Y', 0, 0.0, 0.0, "TAC", "Tyr", "Tyrosine"};
    m_map["TAC"] = aa;
    aa = {'Y', 0, 0.0, 0.0, "TAT", "Tyr", "Tyrosine"};
    m_map["TAT"] = aa;

}


void AminoAcids::add(const std::string &codon, int value)
{
    auto node = m_map.find(codon);
    if (node == m_map.end()) return;
    node->second.count += value;
}


void AminoAcids::reset()
{
    for (auto node : m_map)
        node.second.count = 0;
}


void AminoAcids::setTimeDecoding(const std::string &codon, double valueTime)
{
    auto node = m_map.find(codon);
    if (node == m_map.end()) return;
    node->second.timeDecoding = valueTime;
}


void AminoAcids::setTimePausing(const std::string &codon, double valueTime)
{
    auto node = m_map.find(codon);
    if (node == m_map.end()) return;
    node->second.timePausing = valueTime;
}


void AminoAcids::addTime(const std::string &codon, double timePausing, double timeDecoding)
{
    auto node = m_map.find(codon);
    if (node == m_map.end()) return;
    node->second.timePausing += timePausing;
    node->second.timeDecoding += timeDecoding;
    node->second.count++;
}


double AminoAcids::timeDecoding(const std::string &codon) const
{
    auto node = m_map.find(codon);
    if (node == m_map.end()) return 0.0;
    return node->second.timeDecoding;
}


double AminoAcids::timePausing(const std::string &codon) const
{
    auto node = m_map.find(codon);
    if (node == m_map.end()) return 0.0;
    return node->second.timePausing;
}


void AminoAcids::write()
{
    std::cout << "#codon\tletter\tcode\tname\tcount\ttime.decoding\ttime.pausing" << std::endl;
    for (auto node : m_map)
        std::cout << node.second.codon << "\t"
                  << node.second.letter << "\t"
                  << node.second.code << "\t"
                  << node.second.name << "\t"
                  << node.second.count << "\t"
                  << node.second.timeDecoding << "\t"
                  << node.second.timePausing << std::endl;
}


void AminoAcids::log(const std::string &label)
{
    for (auto node : m_map) {
        std::cout << label << "\t"
                  << node.second.codon << "\t"
                  << node.second.count << "\t"
                  << node.second.timeDecoding << "\t"
                  << node.second.timePausing << std::endl;
    }
}
