#include <iostream>

static const std::string PROGRAM_NAME = "ribotools";
static const std::string PROGRAM_VERSION = "1.0.0";

int MetaGeneMain(int argc, char const *argv[]);
int help(void);
int error(const std::string &message);
int version(void);
int contact(void);

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
        return version();
        
    else if ((subCommand == "-c") || (subCommand == "--contact"))
        return contact();

    else if (subCommand == "metagene")
        return MetaGeneMain(argc - 1, argv + 1);

    else
        return error("unknown command argument " + subCommand);

    return 0;
}

int help()
{
    std::cerr << PROGRAM_NAME << " is a toolset for ribosome footprint analysis." << std::endl;
    std::cerr << "Version: " << PROGRAM_VERSION << std::endl;
    std::cerr << "About: developed in the Scientific Computing Facility at Max-Planck Instittue For Brain Research." << std::endl;
    std::cerr << "Docs: https://ribotools.github.molgen.mpg.de" << std::endl;
    std::cerr << "Code: https://github.molgen.mpg.de/MPIBR-Bioinformatics/ribotools" << std::endl;
    std::cerr << "Mail: sciclist@brain.mpg.de" << std::endl << std::endl;
    std::cerr << "Usage: " << PROGRAM_NAME << " <subcommand> [options]" << std::endl << std::endl;
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

int error(const std::string &message)
{
    std::cerr << PROGRAM_NAME << std::endl;
    std::cerr << '\t' << "Error:: " << message << std::endl;
    std::cerr << std::endl;
    return help();
}

int version()
{
    std::cerr << PROGRAM_NAME << " v" << PROGRAM_VERSION << std::endl;
    return 0;
}

int contact()
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