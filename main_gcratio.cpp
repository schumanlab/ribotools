#include <iostream>
#include <memory>

#include "parserargv.h"
#include "bamhandle.h"

int main_gcration(int argc, const char *argv[])
{
    std::vector<std::shared_ptr<BamHandle>> handlesBam;

    // parse command line parameters
    ParserArgv parser(argc, argv);

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            std::shared_ptr<BamHandle> handle = std::make_shared<BamHandle>(fileNameNext, 255, 0);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::gcratio::error, provide BAM file." << std::endl;
        return 1;
    }

    // calculate gc-content per BAM file
    std::cout << "# name\tread.count\tgc-ratio_mean\tgc-ratio_std" << std::endl;
    for (auto handle : handlesBam) {
        double gc_mean = 0.0;
        double gc_std = 0.0;
        int readCount = handle->calculateGCcontent(gc_mean, gc_std);
        std::cout << handle->name() << "\t" << readCount << "\t" << gc_mean << "\t" << gc_std << std::endl;
    }

    return 0;
}
