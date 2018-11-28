#include <iostream>

#include "version.hpp"
#include "metagene.hpp"

int help(void);

int main(int argc, char const *argv[])
{
    // check for sub-command
    if (argc < 2) 
        return help();
    
    // parse sub-command
    std::string subCommand(argv[1]);

    if ((subCommand == "-h") || (subCommand == "--help"))
        return help();

    else if ((subCommand == "-v") || (subCommand == "--version"))
        return version("ribotools");
        
    else if ((subCommand == "-c") || (subCommand == "--contact"))
        return contact();

    else if (subCommand == "metagene")
    {
        auto obj = MetaGene();
        if(!obj.parse(argc - 1, argv + 1))
            return 0;

        
    }
    else
    {
        std::cerr << "ribotools" << std::endl;
        std::cerr << '\t' << "Error:: unknown subcommand " << subCommand << std::endl;
        std::cerr << std::endl;
        return help();
    }

    return 0;
}

int help()
{
    std::cerr << "ribotools is a toolset for ribosome footprint analysis." << std::endl;
    version("Version: ");
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