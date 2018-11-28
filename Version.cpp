#include "version.h"

int version(const std::string &PROGRAM_NAME)
{
    std::cerr << PROGRAM_NAME << " v" << VERSION << std::endl;
    return 0;
}