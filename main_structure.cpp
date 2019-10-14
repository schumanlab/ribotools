#include <iostream>
#include <memory>
#include <fstream>
#include <numeric>
#include <unordered_map>

#include <htslib/faidx.h>

#include "argumentparser.h"
#include "bamhandle.h"
#include "bedrecord.h"
#include "version.h"

struct sinfo {
    int count;
    double score_sum;
    double score_max;
};

int main_structure(const int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta_sequence;
    std::string fileFasta_structure;
    std::vector<std::string> filesBam;
    std::vector<std::shared_ptr<BamHandle>> handlesBam;

    // prepare argument parser
    auto p = ArgumentParser("pausing", std::string(VERSION), "calculates z-score pausing score per codon");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentRequired("sequence").setKeyShort("-f").setKeyLong("--fasta").setHelp("FASTA file containing protein sequence");
    p.addArgumentRequired("structure").setKeyShort("-s").setKeyLong("--structure").setHelp("FASTA file containing secondary protein structure");
    p.addArgumentPositional("alignment").setCount(-1).setHelp("list of RiboSeq BAM files");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        fileFasta_sequence = p.get<std::string>("sequence");
        fileFasta_structure = p.get<std::string>("structure");
        filesBam = p.get<std::vector<std::string>>("alignment");
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // open BAM handles
    for (const auto &fileName : filesBam) {
        std::shared_ptr<BamHandle> handle = std::make_shared<BamHandle>(fileName, 255, 0);
        handlesBam.emplace_back(handle);
    }

    // open FASTA files
    faidx_t *fhFai_sequence = fai_load(fileFasta_sequence.c_str());
    if (!fhFai_sequence) {
        std::cerr << "ribotools::structure::error, failed to load FASTA file " << fileFasta_sequence << std::endl;
        return 1;
    }

    faidx_t *fhFai_structure = fai_load(fileFasta_structure.c_str());
    if (!fhFai_structure) {
        std::cerr << "ribotools::structure::error, failed to load FASTA file " << fileFasta_structure << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::structure::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // loop over BED record
    std::string line;
    int line_counter = 0;
    std::vector<char> keys = {'C','E','H'};
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::stringstream iss(line);
        iss >> bed;

        // ORF length
        int lengthCodons = bed.cdsSpan / 3;

        // accumulate coverage from all files
        std::vector<int> codons(static_cast<std::size_t>(lengthCodons), 0);
        for (auto handle : handlesBam)
            handle->calculateSiteCoverage(codons, bed.transcript, 0, bed.span, bed.cdsStart, true);

        // filter based on coverage
        double background_average = static_cast<double>(std::accumulate(codons.begin(), codons.end(), 0)) / lengthCodons;
        if (background_average < 0.1) continue;

        // retrive sequence / structure
        char *seq_aa = faidx_fetch_seq(fhFai_sequence, bed.name.c_str(), 0, lengthCodons, &lengthCodons);
        char *seq_st = faidx_fetch_seq(fhFai_structure, bed.name.c_str(), 0, lengthCodons, &lengthCodons);

        // loop over codons
        int idx_aa = 0;
        int idx_st = 0;
        std::unordered_map<char,sinfo> map_structure;

        for (std::vector<int>::const_iterator it = codons.begin(); it != codons.end(); ++it) {

            char aa = 'X';
            char st = 'X';

            if (idx_aa < lengthCodons)
                aa = seq_aa[idx_aa++];

            if (aa == 'X' || aa == '*')
                continue;

            if (idx_st < lengthCodons)
                st = seq_st[idx_st++];

            if (st == 'X')
                continue;

            double score = static_cast<double>(*it) / background_average;

            // add record to map
            std::unordered_map<char, sinfo>::iterator mit = map_structure.find(st);
            if (mit == map_structure.end()) {
                map_structure[st] = {1, score, score};
            }
            else {
                mit->second.count++;
                mit->second.score_sum += score;
                mit->second.score_max = std::max(score, mit->second.score_max);
            }

            //std::cout << idx_aa << "\t" << idx_st << "\t" << aa << "\t" << st << "\t" << score << std::endl;

        }

        // write out result
        std::cout << bed.name << "\t" << lengthCodons;
        for (auto key : keys) {
            auto it = map_structure.find(key);
            if (it == map_structure.end()) {
                std::cout << "\t" << key << "0\t0.0\t0.0";
            }
            else {
                std::cout << "\t"
                          << it->first << "\t"
                          << it->second.count << "\t"
                          << it->second.score_sum << "\t"
                          << it->second.score_max;

            }
        }
        std::cout << std::endl;




        if (seq_aa)
            free(seq_aa);

        if (seq_st)
            free(seq_st);

        line_counter++;

        //break;
    }

    std::cerr << "Line.Counter= " << line_counter << std::endl;

    // destructors
    fhBed.close();
    fai_destroy(fhFai_sequence);
    fai_destroy(fhFai_structure);

    return 0;
}
