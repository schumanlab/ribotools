#include <iostream>

#include <htslib/faidx.h>

#include "parserargv.h"

int parseRates(const std::string &fileName);

int main_codonrate(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::string fileRates;

    // parse command line arguments
    ParserArgv parser(argc, argv);
    if (!(parser.find("-bed") && parser.next(fileBed))) {
        std::cerr << "ribotools::codonrate::error, provide a BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("-fasta") && parser.next(fileFasta))) {
        std::cerr << "ribotools::codonrate::error, provide a BED file." << std::endl;
        return 1;
    }

    if (!(parser.find("-rates") && parser.next(fileRates))) {
        std::cerr << "ribotools::codonrate::error, provide a BED file." << std::endl;
        return 1;
    }


    return 0;
}


int parseRates(const std::string &fileName)
{
    return 0;
}