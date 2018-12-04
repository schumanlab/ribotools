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
    if (m_fhBed != nullptr)
    {
        if (bgzf_close(m_fhBed) != 0)
            error("failed to close BED file");
    }

    /*
    for (auto const &fhGbed: m_fhGbed)
    {
        if (fhGbed != nullptr)
        {
            if (bgzf_close(fhGbed) != 0)
                error("failed to close GBED file");
        }
    }

    for (auto const &fhTabix: m_fhTabix)
    {
        if (fhTabix != nullptr)
        {
            regidx_destroy(fhTabix);
        }
    }
    */

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
    parser.parse(argc, argv);

    m_fileBed = parser.retrieve<std::string>("bed");
    m_filesGbed = parser.retrieve<std::vector<std::string>>("gbed");


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


void MetaGene::open()
{
    // open annotation file
    m_fhBed = bgzf_open(m_fileBed.c_str(), "r");
    if (!m_fhBed)
    {
        error("failed to open BED file " + m_fileBed);
        return;
    }

    BGZF *fh = bgzf_open(m_filesGbed.at(0).c_str(), "r");
    if (fh == nullptr)
    {
        error("failed to open GBED file " + m_filesGbed.at(0));
        return;
    }
    bgzf_close(fh);

    regidx_t * fhx = regidx_init(m_filesGbed.at(0).c_str(), MetaGene::tabixParse, MetaGene::tabixFree, sizeof(char*), NULL);
    if (fhx == nullptr)
    {
        error("failed to open GBED index file ");
        return;
    }
    regidx_destroy(fhx);

/*
    // open coverage files
    for (auto const &fileGbed: m_filesGbed)
    {
        // GBED file
        BGZF *fhGbed = bgzf_open(fileGbed.c_str(), "r");
        if (fhGbed == nullptr)
        {
            error("failed to open GBED file " + fileGbed);
            return;
        }
        m_fhGbed.push_back(fhGbed);

        // TABIX file
        regidx_t *fhTabix = regidx_init(fileGbed.c_str(), MetaGene::tabixParse, MetaGene::tabixFree, sizeof(char*), NULL);
        if (fhTabix == nullptr)
        {
            error("failed to open GBED index file " + fileGbed);
            return;
        }
        m_fhTabix.push_back(fhTabix);
    }
*/

}


void MetaGene::pileup()
{
    open();
    
    //readBed();
}


int MetaGene::tabixParse(const char *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *usr)
{
    // Use the standard parser for CHROM,FROM,TO
    int i, ret = regidx_parse_tab(line,chr_beg,chr_end,reg,NULL,NULL);
    if ( ret!=0 ) return ret;

    // Skip the fields that were parsed above
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    for (i=0; i<3; i++)
    {
        while ( *ss && !isspace(*ss) ) ss++;
        if ( !*ss ) return -2;  // wrong number of fields
        while ( *ss && isspace(*ss) ) ss++;
    }
    if ( !*ss ) return -2;

    // Parse the payload
    char *se = ss;
    while ( *se && !isspace(*se) ) se++;
    char **dat = (char**) payload;
    *dat = (char*) malloc(se-ss+1);
    memcpy(*dat,ss,se-ss+1);
    (*dat)[se-ss] = 0;
    return 0;
}


void MetaGene::tabixFree(void *payload)
{
    char **dat = (char**)payload;
    free(*dat);
}

void MetaGene::readBed()
{
    kstring_t str = {0, 0, NULL};
    uint32_t counterLines = 0;

    while (bgzf_getline(m_fhBed, '\n', &str) > 0)
    {
        auto ss = std::stringstream(str.s);
        auto bed = BedRecord();
        ss >> bed;
        bed.parseExons();
        counterLines++;

        queryBed(bed);
        break;

    }
    free(ks_release(&str));
    std::cout << "Annotation records: " << counterLines << std::endl;
}


void MetaGene::queryBed(const BedRecord &bed)
{
    regitr_t tabixIterator;

    std::cout << bed.chrom << ":" << bed.chromStart << "-" << bed.chromEnd << std::endl;

    if (regidx_overlap(m_fhTabix.at(0), bed.chrom.c_str(), bed.chromStart, bed.chromEnd, &tabixIterator))
    {
        while (REGITR_OVERLAP(tabixIterator, bed.chromStart, bed.chromEnd))
        {
            std::cout << REGITR_START(tabixIterator) << " " << REGITR_END(tabixIterator) << " " << REGITR_PAYLOAD(tabixIterator, char*) << std::endl;
            tabixIterator.i++;
        }
        
    }
    else
    {
        std::cout << "nohit" << std::endl;
    }

    
}