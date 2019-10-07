#include <iostream>
#include <map>
#include <htslib/sam.h>

#include "parserargv.h"
#include "bamhandle.h"

int main_length(const int argc, const char *argv[])
{
    std::vector<BamHandle*> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);
    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            auto handle = new BamHandle(fileNameNext, 255, 0);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::length::error, provide BAM file." << std::endl;
        return 1;
    }

    // accumulate length
    bam1_t *b = bam_init1();
    for (auto handle : handlesBam) {
        int count = 0;
        std::map<int, int> qlen_map;
        while (handle->readBam(b) > 0) {
            int qlen = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
            auto node = qlen_map.find(qlen);
            if (node == qlen_map.end()) {
                qlen_map[qlen] = 1;
            }
            else {
                qlen_map[qlen]++;
            }
                
            count++;
        }
        std::cout << "# " << handle->name() << "\t" << count << std::endl;
        for (auto node : qlen_map)
            std::cout << node.first << "\t" << node.second << std::endl;
    }
    bam_destroy1(b);
    

    // destructors
    for (auto handle : handlesBam)
        delete handle;

    return 0;
}
