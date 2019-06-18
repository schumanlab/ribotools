#include <iostream>

#define VERSION 1.01

int version();
int usage();

int main_asite(int argc, const char *argv[]);
int main_codonfreq(int argc, const char *argv[]);
int main_codonrate(int argc, const char *argv[]);
int main_features(int argc, const char *argv[]);
int main_metagene(int argc, const char *argv[]);
int main_mtdr(int argc, const char *argv[]);
int main_uorfs(int argc, const char *argv[]);
int main_utrseq(int argc, const char *argv[]);

int main(int argc, const char *argv[])
{
    // make sure sub command is present
    if (argc < 2)
        return usage();

    // parse on subcommand
    const std::string subcommand = std::string(argv[1]);

    if (subcommand == "-h" || subcommand == "--help") return usage();
    
    else if (subcommand == "-v" || subcommand == "--version") return version();

    else if (subcommand == "asite") return main_asite(argc - 1, argv + 1);

    else if (subcommand == "codonfreq") return main_codonfreq(argc - 1, argv + 1);

    else if (subcommand == "codonrate") return main_codonrate(argc - 1, argv + 1);

    else if (subcommand == "features") return main_features(argc - 1, argv + 1);
    
    else if (subcommand == "metagene") return main_metagene(argc - 1, argv + 1);

    else if (subcommand == "mtdr") return main_mtdr(argc - 1, argv + 1);

    else if (subcommand == "uorfs") return main_uorfs(argc - 1, argv + 1);

    else if (subcommand == "utrseq") return main_utrseq(argc - 1, argv + 1);

    else {
        std::cerr << "Error, unknown subcommand " << subcommand << std::endl;
        return usage();
    }
    
    return 0;
}




int version()
{
    std::cout << "Ribotools " << VERSION << std::endl;
    std::cout << "  Scientific Computing Facility" << std::endl;
    std::cout << "  Max-Planck Institute For Brain Research" << std::endl;
    std::cout << "  https://gitlab.mpcdf.mpg.de/mpibr/schu/ribotools" << std::endl;
    std::cout << "  bug reports: sciclist@brain.mpg.de" << std::endl;
    return 0;
}


int usage()
{
    std::cout << "Ribotools help" << std::endl;
    return 0;
}