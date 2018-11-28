#include <iostream>
#include "version.h"
#include "metagene.h"


int help_metagene(void);

int MetaGeneMain(int argc, char const *argv[])
{
    std::string fileBed;
    std::vector<std::string> filesBam;
    int nbins = 100;
    float depth = 20.0;
    float boundLower = -0.3;
    float boundUpper = 1.3;


    if (argc == 1)
        return help_metagene();

    std::cout << "METAGENE " << argc << std::endl;
    return 0;
}

int help_metagene(void)
{
    std::cerr << std::endl;
    std::cerr << "Tool: ribotools metagene" << std::endl;
    version("Version: ");
    std::cerr << "Summary: creates coverage histogram by binning transcripts length." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage: ribotools metagene [OPTIONS] -a <bed> -b <bams>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "[REQUIRED]" << std::endl;
    std::cerr << '\t' << "-a | --bed    annotation file in BED12 format." << std::endl;
    std::cerr << '\t' << "-b | --bams   alignment files in BAM format." << std::endl;
    std::cerr << std::endl;
    std::cerr << "[OPTIONS]" << std::endl;
    std::cerr << '\t' << "-n | --nbins  number of histogram bins [100]." << std::endl;
    std::cerr << '\t' << "-l | --lower  lower histogram bound [-0.3]." << std::endl;
    std::cerr << '\t' << "-u | --upper  upper histogram bound [1.3]." << std::endl;
    std::cerr << '\t' << "-d | --depth  minimum depth [20]." << std::endl;
    std::cerr << std::endl;
    return 0;
}
