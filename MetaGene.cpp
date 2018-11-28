#include "metagene.hpp"

MetaGene::MetaGene() :
    m_nbins(100),
    m_depth(20.0),
    m_boundLower(-0.3),
    m_boundUpper(1.3)
{
    
}


MetaGene::~MetaGene()
{
    
}


bool MetaGene::parse(int argc, char const *argv[])
{
    auto parser = ArgumentParser();
    parser.addArgument("-h", "--help", 0);
    parser.addArgument("-v", "--version", 0);
    parser.addArgument("-c", "--contact", 0);
    parser.addArgument("-a", "--bed", 1);
    parser.addArgument("-g", "--gbed", '+');
    parser.addArgument("-n", "--nbins", 1);
    parser.addArgument("-l", "--lower", 1);
    parser.addArgument("-u", "--upper", 1);
    parser.addArgument("d", "--depth", 1);
    parser.appName("ribotools metagene");
    parser.ignoreFirstArgument(true);
    std::cout << parser.usage() << std::endl;
    parser.parse(argc, argv);

    m_fileBed = parser.retrieve<std::string>("bed");
    std::cout << m_fileBed << std::endl;

    /*
    if (argc < 2)
    {
        help();
        return false;
    }
    
    int a = 1;
    while (a < argc)
    {
        // get current argument key
        std::string argumentKey(argv[a]);
        a++;
        std::cout << argumentKey << std::endl;

        
        // check if is a help / version / contact
        if ((argumentKey == "-h") || (argumentKey == "--help"))
        {
            help();
            return false;
        }
        else if((argumentKey == "-v") || (argumentKey == "--version"))
        {
            version("ribotools metagene ");
            return false;
        }
        else if((argumentKey == "-c") || (argumentKey == "--contact"))
        {
            contact();
            return false;
        }

        // get value to current key
        a++;
        if (a >= argc)
        {
            error("missing value for argument " + argumentKey);
            return false;
        }
        std::string argumentValue(argv[a]);

        // add single key-value pairs
        if ((argumentKey == "-a") || (argumentKey == "--bed"))
        {
            m_fileBed = argumentValue;
        }
        else if ((argumentKey == "-n") || (argumentKey == "--nbins"))
        {
            m_nbins = std::stoi(argumentValue);
        }
        else if ((argumentKey == "-l") || (argumentKey == "--lower"))
        {
            m_boundLower = std::stof(argumentValue);
        }
        else if ((argumentKey == "-u") || (argumentKey == "--upper"))
        {
            m_boundUpper = std::stof(argumentValue);
        }
        else if ((argumentKey == "-d") || (argumentKey == "--depth"))
        {
            m_depth = std::stof(argumentValue);
        }
        else if ((argumentKey == "-g") || (argumentKey == "--gbeds"))
        {
            m_filesGbed.push_back(argumentValue);
            while (++a < argc)
            {
                std::string argumentNext(argv[a]);
                if (argumentNext.at(0) == '-')
                {
                    a--;
                    break;
                }
                m_filesGbed.push_back(argumentNext);
            }
        }
        
    }
    */

    return true;
}

void MetaGene::help()
{
    std::cerr << std::endl;
    std::cerr << "Tool: ribotools metagene" << std::endl;
    version("Version: ");
    std::cerr << "Summary: creates coverage histogram by binning transcripts length." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage: ribotools metagene [OPTIONS] -a <bed> -g <gbeds>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required:" << std::endl;
    std::cerr << '\t' << "-a | --bed    annotation file in BED12 format." << std::endl;
    std::cerr << '\t' << "-g | --gbeds   coverage files in GraphBed format." << std::endl;
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


void MetaGene::error(const std::string &errorMessage)
{
    std::cerr << "ribotools metagene" << std::endl;
    std::cerr << '\t' << "Error:: " << errorMessage << std::endl;
}