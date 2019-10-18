#include <iostream>
#include <memory>
#include <fstream>
#include <numeric>
#include <cmath>

#include <htslib/faidx.h>

#include "argumentparser.h"
#include "bamhandle.h"
#include "bedrecord.h"
#include "aminoacidtable.h"
#include "version.h"


int main_pausing(const int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::vector<std::string> filesBam;
    std::vector<std::shared_ptr<BamHandle>> handlesBam;
    int backgroundWindow_basic;
    int backgroundWindow_flank;
    int skipCodons;

    auto p = ArgumentParser("pausing", std::string(VERSION), "calculates z-score pausing score per codon");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentRequired("sequence").setKeyShort("-f").setKeyLong("--fasta").setHelp("FASTA file containing transcript sequence");
    p.addArgumentOptional("basic").setKeyShort("-b").setKeyLong("--basic").setDefaultValue<int>(150).setHelp("number of codons to calculate background score");
    p.addArgumentOptional("flank").setKeyShort("-f").setKeyLong("--flank").setDefaultValue<int>(5).setHelp("number of flanking codons to calculate background score");
    p.addArgumentOptional("skip").setKeyShort("-s").setKeyLong("--skip").setDefaultValue<int>(10).setHelp("number of codons to skip after/before start/stop codon");
    p.addArgumentPositional("alignment").setCount(-1).setHelp("list of RiboSeq BAM files");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        fileFasta = p.get<std::string>("sequence");
        backgroundWindow_basic = p.get<int>("basic");
        backgroundWindow_flank = p.get<int>("flank");
        skipCodons = p.get<int>("skip");
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

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::pausing::error, failed to load FASTA file " << fileFasta << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::pausing::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // print header
    std::cout << "#name\tcodon\tcount\tmin.pausing.score\tmax.pausing.score" << std::endl;

    // loop over BED record
    std::string line;
    //int line_counter = 0;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::stringstream iss(line);
        iss >> bed;

        int lengthCodons = bed.cdsSpan / 3;

        if (lengthCodons <= 2*skipCodons + 2*backgroundWindow_flank + 1) continue;

        // accumulate coverage from all files
        std::vector<int> codons(static_cast<std::size_t>(lengthCodons), 0);
        for (auto handle : handlesBam)
            handle->calculateSiteCoverage(codons, bed.transcript, 0, bed.span, bed.cdsStart, true);

        // fiter based on coverage and empty codons
        double background_average = static_cast<double>(std::accumulate(codons.begin(), codons.end(), 0)) / lengthCodons;
        if (background_average < 0.1) continue;

        // calculate basic background
        int backgroundWindow_offset = std::min(backgroundWindow_basic, lengthCodons - 2*skipCodons);
        double background_basic = static_cast<double>(std::accumulate(codons.begin() + skipCodons,
                                                                      codons.begin() + skipCodons + backgroundWindow_offset,
                                                                      0)) / backgroundWindow_offset;


        // retrieve sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        auto aa = AminoAcidTable();

        // count paused codons
        std::vector<int>::const_iterator it;
        std::vector<int>::const_iterator it_prev;
        std::vector<int>::const_iterator it_next;
        for (it = codons.begin() + skipCodons; it != (codons.end() - skipCodons); ++it) {

            // calculate previous background
            double background_prev = 0.0;
            it_prev = it - backgroundWindow_flank;
            if (codons.begin() <= it_prev)
                background_prev = static_cast<double>(std::accumulate(it_prev, it, background_prev)) / backgroundWindow_flank;

            // calculate next background
            double background_next = 0.0;
            it_next = it + backgroundWindow_flank + 1;
            if ((it + 1) < codons.end() && (it_next <= codons.end()))
                background_next = static_cast<double>(std::accumulate(it + 1, it_next, background_next)) / backgroundWindow_flank;

            // best background
            double background = std::max(background_prev, background_prev);
            background = std::max(background, background_basic);

            // calculate zscore
            double zscore = 0.0;
            if (background > 0.0)
                zscore = (*it - background) / std::sqrt(background);

            // current codon code
            int idx_codon = static_cast<int>(it - codons.begin());
            int idx_nucleotide = idx_codon * 3 + bed.cdsStart;
            char codon[4];
            std::strncpy(codon, &sequence[idx_nucleotide], 3);
            codon[3] = '\0';

            aa.addPausingScore(std::string(codon), zscore);


            //if (zscore != 0.0)
                //std::cout << bed.gene << "\t" << (it - codons.begin()) << "\t" << codon << "\t" << zscore << std::endl;
        }

        // output of AminoAcid map
        aa.log(bed.name);

        if (sequence)
            free(sequence);


    }


    // destructors
    fhBed.close();
    fai_destroy(fhFai);

    return 0;
    /*

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            BamHandle *handle = new BamHandle(fileNameNext, 255, 0);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::pausing::error, provide BAM file." << std::endl;
        return 1;
    }


    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::pausing::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::pausing::error, failed to load fasta reference " << fileFasta << std::endl;
        return 1;
    }

    // output header
    std::cout << "# name\tcodon\tn.codons\tn.asites\tsum.zscore" << std::endl;

    // loop over BED record
    std::string line;
    int line_counter = 0;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        int lengthCodons = bed.cdsSpan / 3;

        // check if codon length is too small
        if (lengthCodons <= 2*skipCodons + 2*backgroundWindow_flank + 1) continue;

        // calculate A-site codon coverage
        std::vector<int> codons(static_cast<size_t>(lengthCodons), 0);
        for (auto handle : handlesBam)
            handle->calculateSiteCoverage(codons, bed.transcript, 0, bed.span, bed.cdsStart, true);

        // count empty codons
        //int emptyCodons = static_cast<int>(std::count(codons.begin(), codons.end(), 0));
        //double emptyRatio = static_cast<double>(emptyCodons) / lengthCodons;

        // get average coverage
        double background_average = static_cast<double>(std::accumulate(codons.begin(), codons.end(), 0)) / lengthCodons;

        // fiter based on coverage and empty codons
        //if (background_average < 0.1 || emptyRatio > 0.75) continue;
        if (background_average < 0.1) continue;


        // calculate background
        int backgrounWindow_offset = std::min(backgroundWindow_basic, lengthCodons - 2*skipCodons);
        double background_basic = static_cast<double>(std::accumulate(codons.begin() + skipCodons, codons.begin() + backgrounWindow_offset, 0)) / backgrounWindow_offset;

        // retrieve sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        std::vector<int>::const_iterator it;
        std::vector<int>::const_iterator it_prev;
        std::vector<int>::const_iterator it_next;
        auto aa = AminoAcids();

        for (it = codons.begin() + skipCodons; it != (codons.end() - skipCodons); ++it)
        {
            // skip empty codons
            if (*it == 0) continue;

            // calculate previous background
            double background_prev = 0.0;
            it_prev = it - backgroundWindow_flank;
            if (codons.begin() <= it_prev)
                background_prev = static_cast<double>(std::accumulate(it_prev, it, background_prev)) / backgroundWindow_flank;

            // calculate next background
            double background_next = 0.0;
            it_next = it + backgroundWindow_flank + 1;
            if ((it + 1) < codons.end() && (it_next <= codons.end()))
                background_next = static_cast<double>(std::accumulate(it+1, it_next, background_next)) / backgroundWindow_flank;

            // background value
            double background = std::max(background_next, background_prev);
            background = std::max(background, background_basic);

            // skip if background is empty
            if (background == 0.0) continue;

            // zscore
            double zscore = (*it - background) / std::sqrt(background);

            // current codon code
            int idx_codon = static_cast<int>(it - codons.begin());
            int idx_nucleotide = idx_codon * 3 + bed.cdsStart;
            char codon[4];
            std::strncpy(codon, &sequence[idx_nucleotide], 3);
            codon[3] = '\0';

            aa.addTime(std::string(codon), zscore, static_cast<double>(*it));

            // debug output
            //std::cout << idx_codon << "\t" << codon << "\t" << *it << "\t" << background << "\t" << zscore << std::endl;
        }

        // output of AminoAcid map
        aa.log(bed.name);

        if (sequence)
            free(sequence);

        line_counter++;
    }

    std::cerr << "used genes: " << line_counter << std::endl;

    // destructors
    fhBed.close();
    fai_destroy(fhFai);
    for (auto handle : handlesBam)
        delete handle;


    return 0;
    */
}
