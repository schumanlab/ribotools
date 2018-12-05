#include <iostream>

#include <chrono>
#include <htslib/tbx.h>
#include <htslib/kstring.h>


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

    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = toc - tic;
    std::cout << "Elapsed time: " << elapsed.count() << " s.\n";

    return 0;
}
