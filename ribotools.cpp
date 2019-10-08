#include <iostream>
#include "version.h"

int main_basefreq(const int argc, const char *argv[]);
int main_codonfreq(const int argc, const char *argv[]);
int main_codonrate(const int argc, const char *argv[]);
int main_count(const int argc, const char *argv[]);
int main_features(const int argc, const char *argv[]);
int main_gcratio(const int argc, const char *argv[]);
int main_gcref(const int argc, const char *argv[]);
int main_irate(const int argc, const char *argv[]);
int main_length(const int argc, const char *argv[]);
int main_metagene(const int argc, const char *argv[]);
int main_mtdr(const int argc, const char *argv[]);
int main_pausing(const int argc, const char *argv[]);
int main_poffset(const int argc, const char *argv[]);
int main_translate(const int argc, const char *argv[]);
int main_uorfs(const int argc, const char *argv[]);
int main_utrseq(const int argc, const char *argv[]);

int printUsage();
int printVersion();
int printContact();

int main(const int argc, const char *argv[]) {

        // make sure sub command is present
        if (argc < 2)
            return printUsage();

        // parse on subcommand
        const std::string subcommand = std::string(argv[1]);

        if (subcommand == "-h" || subcommand == "--help") return printUsage();

        else if (subcommand == "-v" || subcommand == "--version") return printVersion();

        else if (subcommand == "-c" || subcommand == "--contact") return printContact();

        else if (subcommand == "basefreq") return main_basefreq(argc - 1, argv + 1);

        else if (subcommand == "codonfreq") return main_codonfreq(argc - 1, argv + 1);

        else if (subcommand == "codonrate") return main_codonrate(argc - 1, argv + 1);

        else if (subcommand == "count") return main_count(argc - 1, argv + 1);

        else if (subcommand == "features") return main_features(argc - 1, argv + 1);

        else if (subcommand == "gcratio") return main_gcratio(argc - 1, argv + 1);

        else if (subcommand == "gcref") return main_gcref(argc - 1, argv + 1);

        else if (subcommand == "irate") return main_irate(argc - 1, argv + 1);

        else if (subcommand == "length") return main_length(argc - 1, argv + 1);

        else if (subcommand == "metagene") return main_metagene(argc - 1, argv + 1);

        else if (subcommand == "mtdr") return main_mtdr(argc - 1, argv + 1);

        else if (subcommand == "pausing") return main_pausing(argc - 1, argv + 1);

        else if (subcommand == "poffset") return main_poffset(argc - 1, argv + 1);

        else if (subcommand == "translate") return main_translate(argc - 1, argv + 1);

        else if (subcommand == "uorfs") return main_uorfs(argc - 1, argv + 1);

        else if (subcommand == "utrseq") return main_utrseq(argc - 1, argv + 1);

        else {
            std::cerr << "ribotools::error, unknown subcommand " << subcommand << std::endl;
            return printUsage();
        }
}


int printVersion()
{
    std::cout << "version " << VERSION << std::endl;
    return 0;
}


int printContact()
{
    std::cerr << "Scientific Computing Facility" << std::endl;
    std::cerr << "Max-Planck Institute For Brain Research" << std::endl;
    std::cerr << "Frankfurt am Main, Germany" << std::endl;
    std::cerr << "source code: https://gitlab.mpcdf.mpg.de/mpibr/schu/ribotools" << std::endl;
    std::cerr << "bug reports: sciclist@brain.mpg.de" << std::endl;
    std::cerr << "author: georgi a. tushev" << std::endl;
    return 0;
}


int printUsage()
{
    std::cerr << "ribotools" << std::endl;
    std::cerr << "a toolset to profile ana analyse ribosome footprint sequencing" << std::endl;
    printVersion();
    std::cerr << std::endl;
    std::cerr << "usage: ribotools <subcommand> [-h|--help -v|--version -c|--contact]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "[subcommands]" << std::endl;
    std::cerr << "    basefreq       calculates base frequency per ORF" << std::endl;
    std::cerr << "    codonfreq      calculates codon frequency per ORF" << std::endl;
    std::cerr << "    codonrate      calculates codon decoding rate" << std::endl;
    std::cerr << "    count          counts reads per ORF from BAM files" << std::endl;
    std::cerr << "    features       counts reads per gene features from BAM and BED files" << std::endl;
    std::cerr << "    gcratio        calculates GC-ratio in reads from BAM file" << std::endl;
    std::cerr << "    gcref          calculates GC-ratio in ORF from BED and FASTA files" << std::endl;
    std::cerr << "    irate          calculates initiation rate based on BAM files from RNASeq and RiboSeq" << std::endl;
    std::cerr << "    length         calculates reads length based on CIGAR string from BAM file" << std::endl;
    std::cerr << "    metagene       project footprints to a metagene histogram" << std::endl;
    std::cerr << "    mtdr           calculates mean transcript decoding rate (MTDR)" << std::endl;
    std::cerr << "    pausing        calculates z-score pausing score per codon" << std::endl;
    std::cerr << "    poffset        calculates P-site offset per read length" << std::endl;
    std::cerr << "    translate      translates transcripts based on BED and FASTA files" << std::endl;
    std::cerr << "    uorfs          screens for upstream open reading frame" << std::endl;
    std::cerr << "    utrseq         exports UTR sequence from BED and FASTA files" << std::endl;
    std::cerr << std::endl;
    std::cerr << "[options]" << std::endl;
    std::cerr << "    -h, --help     print this help message" << std::endl;
    std::cerr << "    -v, --version  what major.minor.build version of ribotools is used" << std::endl;
    std::cerr << "    -c, --contact  feature requests, bugs, mailing lists, etc." << std::endl;
    std::cerr << std::endl;

    return 0;
}



/*

std::string description();



int main(int argc, const char *argv[])
{
    auto p = ArgumentParser("ribotools", "0.0.1", description());
    p.addArgumentCommand("basefreq").setHelp("calculates basefrequency");
    p.addArgumentCommand("codonfreq").setHelp("calculates codon frequency");
    p.addArgumentCommand("codonrate").setHelp("calculates codon decoding rate");
    p.addArgumentCommand("count").setHelp("counts reads per CDS in BAM files");
    p.addArgumentCommand("features").setHelp("counts reads per gene features in BAM files");
    p.addArgumentCommand("gcratio").setHelp("calculate GC ratio in BAM reads");
    p.addArgumentCommand("gcref").setHelp("calculates GC ratio per gene features");
    p.addArgumentCommand("irate").setHelp("calculates initation rate based on RNASeq and RiboSeq");
    p.addArgumentCommand("length").setHelp("calculate read length based on cigar string");
    p.addArgumentCommand("metagene").setHelp("project footprints to a metagene histogram");
    p.addArgumentCommand("mtdr").setHelp("calculates mean transcript decoding rate (MTDR)");
    p.addArgumentCommand("pausing").setHelp("calculates z-score pausing score per codon");
    p.addArgumentCommand("poffset").setHelp("calculates P-site offset per read length");
    p.addArgumentCommand("translate").setHelp("translate transcripts based on BED and FASTA");
    p.addArgumentCommand("uorfs").setHelp("screens for upstream open reading frames");
    p.addArgumentCommand("utrseq").setHelp("exports UTR sequence from BED and FASTA");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version");

    try {
        p.parse(argc, argv);

        if (p.get<bool>("basefreq")) return main_basefreq(argc - 1, argv + 1, p.version());

        if (p.get<bool>("codonfreq")) return main_codonfreq(argc - 1, argv + 1, p.version());

        if (p.get<bool>("codonrate")) return main_codonrate(argc - 1, argv + 1, p.version());

        if (p.get<bool>("count")) return main_count(argc - 1, argv + 1, p.version());

        if (p.get<bool>("features")) return main_features(argc - 1, argv + 1, p.version());

        if (p.get<bool>("gcratio")) return main_gcratio(argc - 1, argv + 1, p.version());

        if (p.get<bool>("gcref")) return main_gcref(argc - 1, argv + 1, p.version());

        if (p.get<bool>("irate")) return main_irate(argc - 1, argv + 1, p.version());

        if (p.get<bool>("length")) return main_length(argc - 1, argv + 1, p.version());

        if (p.get<bool>("metagene")) return main_metagene(argc - 1, argv + 1, p.version());

        if (p.get<bool>("mtdr")) return main_mtdr(argc - 1, argv + 1, p.version());

        if (p.get<bool>("pausing")) return main_pausing(argc - 1, argv + 1, p.version());

        if (p.get<bool>("poffset")) return main_poffset(argc - 1, argv + 1, p.version());

        if (p.get<bool>("translate")) return main_translate(argc - 1, argv + 1, p.version());

        if (p.get<bool>("uorfs")) return main_uorfs(argc - 1, argv + 1, p.version());

        if (p.get<bool>("utrseq")) return main_utrseq(argc - 1, argv + 1, p.version());

        if (p.get<bool>("help")) {
            std::cerr << p.printHelp() << std::endl;
            return 1;
        }

        if (p.get<bool>("version")) {
            std::cerr << p.printVersion() << std::endl;
            return 1;
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}


std::string description() {
    std::stringstream oss;
    oss << "assembly of tools to profile and analyse" << std::endl;
    oss << "ribosome footprint sequencing" << std::endl;
    oss << "  Scientific Computing Facility" << std::endl;
    oss << "  Max-Planck Institute For Brain Research" << std::endl;
    oss << "  https://gitlab.mpcdf.mpg.de/mpibr/schu/ribotools" << std::endl;
    oss << "  bug reports: sciclist@brain.mpg.de" << std::endl;
    return oss.str();
}
*/
