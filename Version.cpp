#include "version.hpp"

int version(const std::string &PROGRAM_NAME)
{
    std::cerr << PROGRAM_NAME << " v" << VERSION << std::endl;
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