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
int calculateCodonAdaptationIndex(AminoAcidTable &aatable, const std::string &fileBed, const std::string &fileFasta);

int main_cai(const int argc, const char *argv[])
{
    std::string fileBed_query;
    std::string fileBed_reference;
    std::string fileFasta;
    std::vector<std::string> filesBam;

    // arguments
    auto p = ArgumentParser("cai", std::string(VERSION), "calculates codon adaptation index (CAI)");
    p.addArgumentRequired("query").setKeyShort("-q").setKeyLong("--query").setHelp("BED file containing query transcript annotation");
    p.addArgumentRequired("reference").setKeyShort("-r").setKeyLong("--reference").setHelp("BED file containing reference transcript annotation");
    p.addArgumentRequired("sequence").setKeyShort("-f").setKeyLong("--fasta").setHelp("FASTA file containing transcript sequence");
    p.addArgumentOptional("weight").setCount(-1).setKeyShort("-b").setKeyLong("--bam").setHelp("list of BAM files to weight reference codons").setDefaultValue<std::string>("");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed_query = p.get<std::string>("query");
        fileBed_reference = p.get<std::string>("reference");
        fileFasta = p.get<std::string>("sequence");
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
    if (calculateCodonAdaptationIndex(aatable, fileBed_query, fileFasta) != 0) {
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
int calculateCodonAdaptationIndex(AminoAcidTable &aatable, const std::string &fileBed, const std::string &fileFasta)
{
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

        // retrieve sequence
        double cai = 0.0;
        int count = 0;
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        for (int c = bed.cdsStart; c < bed.cdsEnd; c += 3) {
            char codon_seq[4];
            std::strncpy(codon_seq, &sequence[c], 3);
            codon_seq[3] = '\0';
            auto codon = std::string(codon_seq);
            double score = aatable.score(codon);
            if (score > 0.0) {
                cai += log(score);
                count++;
            }
        }

        // geometric mean
        cai = cai / count;
        cai = std::exp(cai);

        // result
        std::cout << bed.gene << "\t" << cai << std::endl;

        if (sequence)
            free(sequence);
    }

    // destructors
    fhBed.close();
    fai_destroy(fhFai);

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
