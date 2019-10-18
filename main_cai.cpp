#include <iostream>
#include <memory>
#include <fstream>
#include <numeric>

#include <htslib/faidx.h>

#include "argumentparser.h"
#include "bamhandle.h"
#include "bedrecord.h"
#include "aminoacidtable.h"
#include "version.h"


int calculateBackgroundFrequency(AminoAcidTable &aatable, const std::string &fileBed, const std::string &fileFasta, const std::vector<std::string> &filesBam);
int calculateCodonAdaptationIndex(AminoAcidTable &aatable, const std::string &fileBed, const std::string &fileFasta, const std::string &fileStructure);

int main_cai(const int argc, const char *argv[])
{
    std::string fileBed_query;
    std::string fileBed_reference;
    std::string fileFasta;
    std::string fileStructure;
    std::vector<std::string> filesBam;

    // arguments
    auto p = ArgumentParser("cai", std::string(VERSION), "calculates codon adaptation index (CAI)");
    p.addArgumentRequired("query").setKeyShort("-q").setKeyLong("--query").setHelp("BED file containing query transcript annotation");
    p.addArgumentRequired("reference").setKeyShort("-r").setKeyLong("--reference").setHelp("BED file containing reference transcript annotation");
    p.addArgumentRequired("sequence").setKeyShort("-f").setKeyLong("--fasta").setHelp("FASTA file containing transcript sequence");
    p.addArgumentOptional("structure").setKeyShort("-s").setKeyLong("--structure").setHelp("FASTA file containing protein structure").setDefaultValue<std::string>("");
    p.addArgumentOptional("weight").setCount(-1).setKeyShort("-b").setKeyLong("--bam").setHelp("list of BAM files to weight reference codons").setDefaultValue<std::string>("");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed_query = p.get<std::string>("query");
        fileBed_reference = p.get<std::string>("reference");
        fileFasta = p.get<std::string>("sequence");
        fileStructure = p.get<std::string>("structure");
        filesBam = p.get<std::vector<std::string>>("weight");
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // calculate background frequency
    AminoAcidTable aatable;
    if (calculateBackgroundFrequency(aatable, fileBed_reference, fileFasta, filesBam) != 0) {
        std::cerr << "ribotools::cai::error, failed to calculate codons background frequency" << std::endl;
        return 1;
    }

    // calculate relative synonymous codon usage (RSCU)
    aatable.calculateRSCU();
    //aatable.write();

    // calculate codon adaptation index CAI
    if (calculateCodonAdaptationIndex(aatable, fileBed_query, fileFasta, fileStructure) != 0) {
        std::cerr << "ribotools::cai::error, failed to calculate codon adaptation index" << std::endl;
        return 1;
    }

    return 0;
}


/**
 * @brief calculateCodonAdaptationIndex
 * @param aatable amino acid table class
 * @param fileBed query bed file for CAI score
 * @param fileFasta DNA sequence
 * @return 0 for success
 */
int calculateCodonAdaptationIndex(AminoAcidTable &aatable, const std::string &fileBed, const std::string &fileFasta, const std::string &fileStructure)
{
    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::cai::error, failed to load FASTA file " << fileFasta << std::endl;
        return 1;
    }

    // open structure file
    faidx_t *fhFas = nullptr;
    if (!fileStructure.empty()) {
        fhFas = fai_load(fileStructure.c_str());
        if (!fhFas) {
            std::cerr << "ribotools::cai::error, failed to laod FASTA structure file " << fileStructure << std::endl;
            return 1;
        }
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::cai::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // loop over BED record
    std::string line;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::stringstream iss(line);
        iss >> bed;

        // retrieve sequence
        double cai = 0.0;
        double cai_C = 0.0;
        double cai_E = 0.0;
        double cai_H = 0.0;

        int count = 0;
        int count_C = 0;
        int count_E = 0;
        int count_H = 0;

        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);

        char *structure = nullptr;
        int indexStructure = 0;
        if (fhFas) {
            int codonSpan = bed.cdsSpan / 3;
            structure = faidx_fetch_seq(fhFas, bed.name.c_str(), 0, codonSpan, &codonSpan);
        }

        for (int c = bed.cdsStart; c < bed.cdsEnd; c += 3) {
            char codon_seq[4];
            std::strncpy(codon_seq, &sequence[c], 3);
            codon_seq[3] = '\0';
            auto codon = std::string(codon_seq);
            double score = aatable.score(codon);
            if (score > 0.0) {
                cai += std::log(score);
                count++;

                // use structure
                if (structure) {
                    char AA = structure[indexStructure];
                    indexStructure++;

                    if (AA == 'C') {
                        cai_C += std::log(score);
                        count_C++;
                    }

                    if (AA == 'E') {
                        cai_E += std::log(score);
                        count_E++;
                    }

                    if (AA == 'H') {
                        cai_H += std::log(score);
                        count_H++;
                    }
                }


            }
        }

        // geometric mean
        cai = cai / count;
        cai = std::exp(cai);

        // result
        std::cout << bed.gene << "\t" << cai;
        if (structure) {

            if (count_C == 0) {
                cai_C = 0;
            }
            else {
                cai_C = cai_C / count_C;
                cai_C = std::exp(cai_C);
            }
            std::cout << "\t" << cai_C;

            if (count_E == 0) {
                cai_E = 0;
            }
            else {
                cai_E = cai_E / count_E;
                cai_E = std::exp(cai_E);
            }
            std::cout << "\t" << cai_E;

            if (count_H == 0) {
                cai_H = 0;
            }
            else {
                cai_H = cai_H / count_H;
                cai_H = std::exp(cai_H);
            }
            std::cout << "\t" << cai_H;
        }
        std::cout << std::endl;

        if (sequence)
            free(sequence);

        if (structure)
            free(structure);
    }

    // destructors
    fhBed.close();

    if (fhFai)
        fai_destroy(fhFai);

    if (fhFas)
        fai_destroy(fhFas);

    return 0;
}


/**
 * @brief calculateBackgroundFrequency
 * @param aatable amino acid table class
 * @param fileBed reference bed file for background frequency population
 * @param fileFasta DNA sequence
 * @param filesBam bam files to calculate coverage weight
 * @return 0 for success
 */

int calculateBackgroundFrequency(AminoAcidTable &aatable, const std::string &fileBed, const std::string &fileFasta, const std::vector<std::string> &filesBam)
{
    std::vector<std::shared_ptr<BamHandle>> handlesBam;

    // open BAM handles
    for (const auto &fileName : filesBam) {
        if (!fileName.empty()) {
            std::shared_ptr<BamHandle> handle = std::make_shared<BamHandle>(fileName, 255, 0);
            handlesBam.emplace_back(handle);
        }
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::cai::error, failed to load FASTA file " << fileFasta << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::cai::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // loop over BED record
    std::string line;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::stringstream iss(line);
        iss >> bed;

        // calculate record weight
        double weight = 1.0;
        if (!handlesBam.empty()) {
            weight = 0.0;
            for (auto handle : handlesBam)
                weight += static_cast<double>(handle->readsPerRegion(bed.transcript, bed.cdsStart, bed.cdsEnd));
            //weight = std::log(weight);
        }
        weight = weight / 1e6;
        //std::cout << bed.gene << "\t" << weight << std::endl;


        // retrieve sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        for (int c = bed.cdsStart; c < bed.cdsEnd; c += 3) {
            char codon_seq[4];
            std::strncpy(codon_seq, &sequence[c], 3);
            codon_seq[3] = '\0';
            auto codon = std::string(codon_seq);
            aatable.addParameters(codon, 1, weight);
        }

        if (sequence)
            free(sequence);
    }

    // destructors
    fhBed.close();
    fai_destroy(fhFai);

    return 0;
}
