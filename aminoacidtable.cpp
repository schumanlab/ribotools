#include "aminoacidtable.h"

AminoAcidTable::AminoAcidTable()
{
    m_table["TAA"] = std::make_unique<AminoAcid>('*', "TAA", "Ter", "Stop_Codon");
    m_table["TAG"] = std::make_shared<AminoAcid>('*', "TAG", "Ter", "Stop_Codon");
    m_table["TGA"] = std::make_shared<AminoAcid>('*', "TGA", "Ter", "Stop_Codon");
    m_table["GCA"] = std::make_shared<AminoAcid>('A', "GCA", "Ala", "Alanine");
    m_table["GCC"] = std::make_shared<AminoAcid>('A', "GCC", "Ala", "Alanine");
    m_table["GCG"] = std::make_shared<AminoAcid>('A', "GCG", "Ala", "Alanine");
    m_table["GCT"] = std::make_shared<AminoAcid>('A', "GCT", "Ala", "Alanine");
    m_table["TGC"] = std::make_shared<AminoAcid>('C', "TGC", "Cys", "Cysteine");
    m_table["TGT"] = std::make_shared<AminoAcid>('C', "TGT", "Cys", "Cysteine");
    m_table["GAC"] = std::make_shared<AminoAcid>('D', "GAC", "Asp", "Aspartic_Acid");
    m_table["GAT"] = std::make_shared<AminoAcid>('D', "GAT", "Asp", "Aspartic_Acid");
    m_table["GAA"] = std::make_shared<AminoAcid>('E', "GAA", "Glu", "Glutamic_Acid");
    m_table["GAG"] = std::make_shared<AminoAcid>('E', "GAG", "Glu", "Glutamic_Acid");
    m_table["TTC"] = std::make_shared<AminoAcid>('F', "TTC", "Phe", "Phenylalanine");
    m_table["TTT"] = std::make_shared<AminoAcid>('F', "TTT", "Phe", "Phenylalanine");
    m_table["GGA"] = std::make_shared<AminoAcid>('G', "GGA", "Gly", "Glycine");
    m_table["GGC"] = std::make_shared<AminoAcid>('G', "GGC", "Gly", "Glycine");
    m_table["GGG"] = std::make_shared<AminoAcid>('G', "GGG", "Gly", "Glycine");
    m_table["GGT"] = std::make_shared<AminoAcid>('G', "GGT", "Gly", "Glycine");
    m_table["CAC"] = std::make_shared<AminoAcid>('H', "CAC", "His", "Histidine");
    m_table["CAT"] = std::make_shared<AminoAcid>('H', "CAT", "His", "Histidine");
    m_table["ATA"] = std::make_shared<AminoAcid>('I', "ATA", "Ile", "Isoleucine");
    m_table["ATC"] = std::make_shared<AminoAcid>('I', "ATC", "Ile", "Isoleucine");
    m_table["ATT"] = std::make_shared<AminoAcid>('I', "ATT", "Ile", "Isoleucine");
    m_table["AAA"] = std::make_shared<AminoAcid>('K', "AAA", "Lys", "Lysine");
    m_table["AAG"] = std::make_shared<AminoAcid>('K', "AAG", "Lys", "Lysine");
    m_table["CTA"] = std::make_shared<AminoAcid>('L', "CTA", "Leu", "Leucine");
    m_table["CTC"] = std::make_shared<AminoAcid>('L', "CTC", "Leu", "Leucine");
    m_table["CTG"] = std::make_shared<AminoAcid>('L', "CTG", "Leu", "Leucine");
    m_table["CTT"] = std::make_shared<AminoAcid>('L', "CTT", "Leu", "Leucine");
    m_table["TTA"] = std::make_shared<AminoAcid>('L', "TTA", "Leu", "Leucine");
    m_table["TTG"] = std::make_shared<AminoAcid>('L', "TTG", "Leu", "Leucine");
    m_table["ATG"] = std::make_shared<AminoAcid>('M', "ATG", "Met", "Methionine");
    m_table["AAC"] = std::make_shared<AminoAcid>('N', "AAC", "Asn", "Asparagine");
    m_table["AAT"] = std::make_shared<AminoAcid>('N', "AAT", "Asn", "Asparagine");
    m_table["CCA"] = std::make_shared<AminoAcid>('P', "CCA", "Pro", "Proline");
    m_table["CCC"] = std::make_shared<AminoAcid>('P', "CCC", "Pro", "Proline");
    m_table["CCG"] = std::make_shared<AminoAcid>('P', "CCG", "Pro", "Proline");
    m_table["CCT"] = std::make_shared<AminoAcid>('P', "CCT", "Pro", "Proline");
    m_table["CAA"] = std::make_shared<AminoAcid>('Q', "CAA", "Gln", "Glutamine");
    m_table["CAG"] = std::make_shared<AminoAcid>('Q', "CAG", "Gln", "Glutamine");
    m_table["AGA"] = std::make_shared<AminoAcid>('R', "AGA", "Arg", "Arginine");
    m_table["AGG"] = std::make_shared<AminoAcid>('R', "AGG", "Arg", "Arginine");
    m_table["CGA"] = std::make_shared<AminoAcid>('R', "CGA", "Arg", "Arginine");
    m_table["CGC"] = std::make_shared<AminoAcid>('R', "CGC", "Arg", "Arginine");
    m_table["CGG"] = std::make_shared<AminoAcid>('R', "CGG", "Arg", "Arginine");
    m_table["CGT"] = std::make_shared<AminoAcid>('R', "CGT", "Arg", "Arginine");
    m_table["AGC"] = std::make_shared<AminoAcid>('S', "AGC", "Ser", "Serine");
    m_table["AGT"] = std::make_shared<AminoAcid>('S', "AGT", "Ser", "Serine");
    m_table["TCA"] = std::make_shared<AminoAcid>('S', "TCA", "Ser", "Serine");
    m_table["TCC"] = std::make_shared<AminoAcid>('S', "TCC", "Ser", "Serine");
    m_table["TCG"] = std::make_shared<AminoAcid>('S', "TCG", "Ser", "Serine");
    m_table["TCT"] = std::make_shared<AminoAcid>('S', "TCT", "Ser", "Serine");
    m_table["ACA"] = std::make_shared<AminoAcid>('T', "ACA", "Thr", "Threonine");
    m_table["ACC"] = std::make_shared<AminoAcid>('T', "ACC", "Thr", "Threonine");
    m_table["ACG"] = std::make_shared<AminoAcid>('T', "ACG", "Thr", "Threonine");
    m_table["ACT"] = std::make_shared<AminoAcid>('T', "ACT", "Thr", "Threonine");
    m_table["GTA"] = std::make_shared<AminoAcid>('V', "GTA", "Val", "Valine");
    m_table["GTC"] = std::make_shared<AminoAcid>('V', "GTC", "Val", "Valine");
    m_table["GTG"] = std::make_shared<AminoAcid>('V', "GTG", "Val", "Valine");
    m_table["GTT"] = std::make_shared<AminoAcid>('V', "GTT", "Val", "Valine");
    m_table["TGG"] = std::make_shared<AminoAcid>('W', "TGG", "Trp", "Tryptophan");
    m_table["TAC"] = std::make_shared<AminoAcid>('Y', "TAC", "Tyr", "Tyrosine");
    m_table["TAT"] = std::make_shared<AminoAcid>('Y', "TAT", "Tyr", "Tyrosine");
}


void AminoAcidTable::reset()
{
    for (auto aa : m_table) {
        aa.second->count = 0;
        aa.second->value = 0.0;
        aa.second->score = 0.0;
    }
}


void AminoAcidTable::setCount(const std::string &codon, int count)
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return;
    aa->second->count = count;
}

int AminoAcidTable::count(const std::string &codon) const
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return 0;
    return aa->second->count;
}


void AminoAcidTable::setValue(const std::string &codon, double value)
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return;
    aa->second->value = value;
}


double AminoAcidTable::value(const std::string &codon) const
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return 0.0;
    return aa->second->value;
}


void AminoAcidTable::setScore(const std::string &codon, double score)
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return;
    aa->second->score = score;
}


double AminoAcidTable::score(const std::string &codon) const
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return 0.0;
    return aa->second->score;
}


void AminoAcidTable::addParameters(const std::string &codon, int count, double value, double score)
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return;
    aa->second->count += count;
    aa->second->value += value;
    aa->second->score += score;
}


void AminoAcidTable::addPausingScore(const std::string &codon, double zscore)
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return;

    /*
    if (zscore > 0.0) {
        aa->second->score += zscore;
        aa->second->count++;
    }
    else {
        aa->second->value += zscore;
    }
    */
    aa->second->count++;
    aa->second->score = std::max(zscore, aa->second->score);
    //aa->second->score += zscore;
    aa->second->value = std::min(zscore, aa->second->value);
}



char AminoAcidTable::translate(const std::string &codon) const
{
    auto aa = m_table.find(codon);
    if (aa == m_table.end()) return 'X';
    return aa->second->letter;
}


void AminoAcidTable::write()
{
    std::cout << "#codon\tletter\tcode\tname\tcount\tvalue\tscore" << std::endl;
    for (auto aa : m_table)
        std::cout << aa.second->codon << "\t"
                  << aa.second->letter << "\t"
                  << aa.second->code << "\t"
                  << aa.second->name << "\t"
                  << aa.second->count << "\t"
                  << aa.second->value << "\t"
                  << aa.second->score << std::endl;
}


void AminoAcidTable::log(const std::string &label)
{
    for (auto aa : m_table) {
        std::cout << label << "\t"
                  << aa.second->codon << "\t"
                  << aa.second->count << "\t"
                  << aa.second->value << "\t"
                  << aa.second->score << std::endl;
    }
}


void AminoAcidTable::calculateRSCU()
{
    std::multimap<char, std::shared_ptr<AminoAcid>> table_synonyms;
    std::unordered_map<char, double> table_maxscore;

    // fill synonyms table
    double weight_total = 0.0;
    for (auto aa : m_table) {
        char letter = aa.second->letter;
        table_synonyms.insert(std::pair<char, std::shared_ptr<AminoAcid>>(letter, aa.second));
        weight_total += aa.second->value;
    }

    // relative frequency
    for (auto aa : m_table)
        aa.second->value = static_cast<double>(aa.second->value) / weight_total;

    // find maximum value per synonym
    std::multimap<char, std::shared_ptr<AminoAcid>>::iterator it_synonym;
    std::unordered_map<char, double>::iterator it_maxscore;

    for(it_synonym = table_synonyms.begin(); it_synonym != table_synonyms.end(); ++it_synonym) {

        it_maxscore = table_maxscore.find(it_synonym->first);
        if (it_maxscore == table_maxscore.end()) {
            table_maxscore[it_synonym->first] = it_synonym->second->value;
        }
        else {
            table_maxscore[it_synonym->first] = std::max(it_synonym->second->value, it_maxscore->second);
        }

    }

    // calculate score based on max table
    for (auto aa : m_table) {
        it_maxscore = table_maxscore.find(aa.second->letter);
        if (it_maxscore == table_maxscore.end()) continue;
        aa.second->score = aa.second->value / it_maxscore->second;
    }

}


void AminoAcidTable::load(const std::string &fileName)
{
    std::ifstream fh;
    fh.open(fileName);
    if (!fh.is_open()) {
        std::cerr << "AminoAcidTable::error, failed to open amino acid table " << fileName << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fh, line)) {
        std::istringstream iss(line);
        std::string codon;
        int count;
        double value;
        double score;
        iss >> codon >> count >> value >> score;
        if (iss.fail()) continue;

        setCount(codon, count);
        setValue(codon, value);
        setScore(codon, score);
    }

    fh.close();
}
