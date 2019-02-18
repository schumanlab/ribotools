#include <iostream>

#include "version.hpp"

int help();
int metageneMain(int argc, char const *argv[]);
int extractumi(int argc, char const *argv[]);
int countasite(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{

    // no sub-commands
    if (argc == 1)
        return help();

    // check first argument for subcommand
    std::string subCommand(argv[1]);

    if ((subCommand == "-h") || (subCommand == "--help")) {
        return help();
    }
    else if ((subCommand == "-v") || (subCommand == "--version")) {
        return version("ribotools");
    }
    else if ((subCommand == "-c") || (subCommand == "--contact")) {
        return contact();
    }
    else if (subCommand == "metagene") {
        return metageneMain(argc - 1, argv + 1);
    }
    else if (subCommand == "extractumi") {
        return extractumi(argc - 1, argv + 1);
    }
    else if (subCommand == "countasite") {
        return countasite(argc - 1, argv + 1);
    }
    else {
        std::cerr << "ribotools" << std::endl;
        std::cerr << '\t' << "Error:: unknown subcommand " << subCommand << std::endl;
        return help();
    }

    return 0;
}

int help()
{
    std::cerr << std::endl;
    std::cerr << "ribotools is a toolset for ribosome footprint analysis." << std::endl;
    version("Version: ");
    std::cerr << std::endl;
    std::cerr << "About: developed in the Scientific Computing Facility at Max-Planck Instittue For Brain Research." << std::endl;
    std::cerr << "Docs: https://ribotools.github.molgen.mpg.de" << std::endl;
    std::cerr << "Code: https://github.molgen.mpg.de/MPIBR-Bioinformatics/ribotools" << std::endl;
    std::cerr << "Mail: sciclist@brain.mpg.de" << std::endl << std::endl;
    std::cerr << "Usage: ribotools <subcommand> [options]" << std::endl << std::endl;
    std::cerr << "The ribotools sub-commands include:" << std::endl;
    std::cerr << '\t' << "metagene" << '\t' << "creates a normalized histogram of footprint coverage" << std::endl;
    std::cerr << std::endl;
    std::cerr << "General help:" << std::endl;
    std::cerr << '\t' << "-h | --help" << '\t' << "Print this help menu." << std::endl;
    std::cerr << '\t' << "-v | --version" << '\t' << "What version of ribotools are you using?" << std::endl;
    std::cerr << '\t' << "-c | --contact" << '\t' << "Feature requests, bugs, mailing lists, etc." << std::endl;
    std::cerr << std::endl;

    return 0;
}
