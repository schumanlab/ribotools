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

int main_cai(const int argc, const char *argv[])
{
    std::string fileBed_query;
    std::string fileBed_reference;
    std::string fileFasta;
    std::vector<std::string> filesBam;
    std::vector<std::shared_ptr<BamHandle>> handlesBam;

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
    fhBed.open(fileBed_reference);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::cai::error, failed to open BED file " << fileBed_reference << std::endl;
        return 1;
    }

    // loop over BED record
    AminoAcidTable aa;
    //aa.test();
    std::string line;
    //int line_counter = 0;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::stringstream iss(line);
        iss >> bed;

        // calculate record weight
        int weight = 1;
        if (!handlesBam.empty()) {
            weight = 0;
            for (auto handle : handlesBam)
                weight += handle->readsPerRegion(bed.transcript, bed.cdsStart, bed.cdsEnd);
        }

        // retrieve sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        for (int c = bed.cdsStart; c < bed.cdsEnd; c += 3) {
            char codon_seq[4];
            std::strncpy(codon_seq, &sequence[c], 3);
            codon_seq[3] = '\0';
            auto codon = std::string(codon_seq);
            aa.addParameters(codon, weight);

        }

        if (sequence)
            free(sequence);

        //std::cout << bed.gene << "\t" << weight << std::endl;
    }

    aa.calculateRSCU();
    aa.write();


    // destructors
    fhBed.close();
    fai_destroy(fhFai);



    return 0;
}
