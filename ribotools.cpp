#include <iostream>
#include "argumentparser.h"

std::string description();

int main_basefreq(int argc, const char *argv[]);
int main_codonfreq(int argc, const char *argv[]);
int main_codonrate(int argc, const char *argv[]);
int main_count(int argc, const char *argv[]);
int main_features(int argc, const char *argv[]);
int main_gcratio(int argc, const char *argv[]);
int main_gcref(int argc, const char *argv[]);
int main_irate(int argc, const char *argv[]);
int main_length(int argc, const char *argv[]);
int main_metagene(int argc, const char *argv[]);
int main_mtdr(int argc, const char *argv[]);
int main_pausing(int argc, const char *argv[]);
int main_poffset(int argc, const char *argv[]);
int main_translate(int argc, const char *argv[]);
int main_uorfs(int argc, const char *argv[]);
int main_utrseq(int argc, const char *argv[]);

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
    p.addArgumentCommand("translage").setHelp("translate transcripts based on BED and FASTA");
    p.addArgumentCommand("uorfs").setHelp("screens for upstream open reading frames");
    p.addArgumentCommand("utrseq").setHelp("exports UTR sequence from BED and FASTA");
    p.addArgumentFlag("help message").setKeyShort("-h").setKeyLong("--help");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version");

    try {
        p.parse(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }



    /*
    // parse on subcommand
    const std::string subcommand = std::string(argv[1]);

    if (subcommand == "-h" || subcommand == "--help") return usage();
    
    else if (subcommand == "-v" || subcommand == "--version") return version();

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
        std::cerr << "Error, unknown subcommand " << subcommand << std::endl;
        return usage();
    }
    */
    
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
