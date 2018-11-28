#include <iostream>
#include "version.h"

int help_ribotools(void);
int error_ribotools(const std::string &message);
int contact_ribotools(void);

int MetaGeneMain(int argc, char const *argv[]);


int main(int argc, char const *argv[])
{
    // check for sub-command
    if (argc < 2) 
        return help_ribotools();
    
    // parse sub-command
    std::string subCommand(argv[1]);

    if ((subCommand == "-h") || (subCommand == "--help"))
        return help_ribotools();

    else if ((subCommand == "-v") || (subCommand == "--version"))
        return version("ribotools");
        
    else if ((subCommand == "-c") || (subCommand == "--contact"))
        return contact_ribotools();

    else if (subCommand == "metagene")
        return MetaGeneMain(argc - 1, argv + 1);

    else
        return error_ribotools("unknown command argument " + subCommand);

    return 0;
}

int help_ribotools()
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

int error_ribotools(const std::string &message)
{
    std::cerr << "ribotools" << std::endl;
    std::cerr << '\t' << "Error:: " << message << std::endl;
    std::cerr << std::endl;
    return help_ribotools();
}



int contact_ribotools()
{
    std::cerr << std::endl;
    std::cerr << "- for further help or to report a bug, please email:" << std::endl;
    std::cerr << '\t' << "sciclist@brain.mpg.de" << std::endl;

    std::cerr << std::endl;
    std::cerr << "- development repository can be found at:" << std::endl;
    std::cerr << '\t' << "https://github.molgen.mpg.de/MPIBR-Bioinformatics/ribotools" << std::endl;

    std::cerr << std::endl;
    std::cerr << "- stable releases can be found at:" << std::endl;
    std::cerr << '\t' << "https://github.molgen.mpg.de/MPIBR-Bioinformatics/ribotools/releases" << std::endl;

    std::cerr << std::endl;

    return 0;
}