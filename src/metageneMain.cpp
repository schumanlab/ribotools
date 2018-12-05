#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>

#include "version.hpp"
#include "parserbed.hpp"
#include "bedrecord.hpp"

class MetaGeneConfig
{
public:
    MetaGeneConfig(int argc = 0, char const *argv[] = nullptr);

    bool errorFlag;
    int nbins;
    float boundLower;
    float boundUpper;
    std::string fileBed;
    std::vector<std::string> filesGbed;

    void help(void);
    void error(const std::string &message);
};


int metageneMain(int argc, char const *argv[])
{
    auto config = MetaGeneConfig(argc, argv);
    auto parser = ParserBed(config.fileBed);

    auto tic = std::chrono::high_resolution_clock::now();
    int lines = 0;
    while (parser.next() > 0)
    {
        auto ss = std::stringstream(parser.buffer.s);
        auto bed = BedRecord();
        ss >> bed;
        //std::cout << bed.name << std::endl;
        lines++;
    }

    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = toc - tic;
    std::cout << "ribotools metagene " << elapsed.count() << " s.\n";
    std::cout << "BEDLINES: " << lines << std::endl;

    return 0;
}


MetaGeneConfig::MetaGeneConfig(int argc, char const *argv[]) :
        errorFlag(false),
        nbins(100),
        boundLower(-0.3),
        boundUpper(1.3),
        fileBed("")
{
    // parse parameters
    bool requiredBed = false;
    bool requiredGbed = false;
    int p = 1;
    while (p < argc)
    {
        // next key
        std::string argumentKey(argv[p]);
        if (argumentKey.at(0) != '-')
        {
            error("provided argument '" + argumentKey + "' is not a valid key."); 
            return;   
        }

        // check for keys without a value
        if ((argumentKey == "-h") || (argumentKey == "--help"))
        {
            help();
            return;
        }
        else if ((argumentKey == "-v") || (argumentKey == "--version"))
        {
            version("ribotools metagene");
            return;
        }
        else if ((argumentKey == "-c") || (argumentKey == "--contact"))
        {
            contact();
            return;
        }

        // increment parameter index
        p++;
        if (p >= argc)
        {
            error("provided key '" + argumentKey + "' requires a value.");
            return;
        }

        // check for keys with value
        std::string argumentValue(argv[p]);
        if ((argumentKey == "-a") || (argumentKey == "--bed"))
        {
            fileBed = argumentValue;
            requiredBed = true;
        }
        else if ((argumentKey == "-g") || (argumentKey == "--gbed"))
        {
            filesGbed.push_back(argumentValue);
            int n = p + 1;
            while( n < argc)
            {
                std::string argumentValueNext(argv[n]);
                if (argumentValueNext.at(0) == '-')
                    break;
                filesGbed.push_back(argumentValueNext);
                n++;
            }
            p += filesGbed.size() - 1;
            requiredGbed = true;
        }
        else if ((argumentKey == "-n") || (argumentKey == "--nbins"))
        {
            nbins = std::stoi(argumentValue);
        }
        else if ((argumentKey == "-l") || (argumentKey == "--lower"))
        {
            boundLower = std::stof(argumentValue);
        }
        else if ((argumentKey == "-u") || (argumentKey == "--upper"))
        {
            boundUpper = std::stof(argumentValue);
        }
        else
        {
            error("unknown key '" + argumentKey + "'.");
            return;
        }
        p++;
    }

    // check required parameters
    if (requiredBed == false)
    {
        error("required annotation file is missing.");
        return;
    }

    if (requiredGbed == false)
    {
        error("required coverage file is missing.");
        return;
    }

};


void MetaGeneConfig::help()
{
    std::cerr << std::endl;
    std::cerr << "Tool: ribotools metagene" << std::endl;
    version("Version: ");
    std::cerr << std::endl;
    std::cerr << "Summary: creates coverage histogram by binning transcripts length." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage: ribotools metagene [OPTIONS] -a <bed> -g <gbeds>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required:" << std::endl;
    std::cerr << '\t' << "-a | --bed    annotation file in BED12 format." << std::endl;
    std::cerr << '\t' << "-g | --gbed   coverage files in GraphBed format." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << '\t' << "-n | --nbins  number of histogram bins [100]." << std::endl;
    std::cerr << '\t' << "-l | --lower  lower histogram bound [-0.3]." << std::endl;
    std::cerr << '\t' << "-u | --upper  upper histogram bound [1.3]." << std::endl;
    std::cerr << '\t' << "-d | --depth  minimum depth [20]." << std::endl;
    std::cerr << std::endl;
    std::cerr << "General help:" << std::endl;
    std::cerr << '\t' << "-h | --help" << '\t' << "Print this help menu." << std::endl;
    std::cerr << '\t' << "-v | --version" << '\t' << "What version of ribotools are you using?" << std::endl;
    std::cerr << '\t' << "-c | --contact" << '\t' << "Feature requests, bugs, mailing lists, etc." << std::endl;
    std::cerr << std::endl;
}


void MetaGeneConfig::error(const std::string &message)
{
    std::cerr << std::endl;
    std::cerr << "ribotools metagene" << std::endl;
    std::cerr << '\t' << "Error: " << message << std::endl;
    errorFlag = true;
    help();
}

