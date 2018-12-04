#include <iostream>

#include <chrono>

#include "version.hpp"
#include "metagene.hpp"

#include <htslib/tbx.h>

int help(void);

int main(int argc, char const *argv[])
{
    std::string fileName = "/Users/tushevg/Desktop/RiboData/gbed/MonoVsPoly/Total_01.gbed.gz";
    std::string query = "chr18:56193977-56295869";



    auto tic = std::chrono::high_resolution_clock::now();

    htsFile *fp = hts_open(fileName.c_str(), "r");
    if (!fp)
    {
        std::cerr << "ERROR: could not read hts file." << std::endl;
        return 1;
    }

    tbx_t* tbx = tbx_index_load(fileName.c_str());
    if (!tbx)
    {
        std::cerr << "ERROR: could not read index." << std::endl;
        return 1;
    }

    // create iterator
    hts_itr_t *iter = tbx_itr_querys(tbx, query.c_str());
    if (!iter)
    {
        std::cerr << "ERROR: could not query the region." << std::endl;
        return 1;
    }

    kstring_t str;
    while (tbx_itr_next(fp, tbx, iter, &str) >= 0)
    {
        std::string line(str.s);
        //std::cout << "Line: " << line << std::endl;
        iter->i++;
    }
    

    free(ks_release(&str));
    hts_close(fp);
    tbx_itr_destroy(iter);
    tbx_destroy(tbx);
    //BGZF *fh = bgzf_open(fileName.c_str(), "r");
    //bgzf_close(fh); 

    

    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = toc - tic;
    std::cout << "Elapsed time: " << elapsed.count() << " s.\n";


    /*
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
        obj.pileup();

        
    }
    else
    {
        std::cerr << "ribotools" << std::endl;
        std::cerr << '\t' << "Error:: unknown subcommand " << subCommand << std::endl;
        std::cerr << std::endl;
        return help();
    }
    */

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