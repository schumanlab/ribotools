#include <iostream>
#include <memory>

#include "parserargv.h"
#include "bamhandle.h"

int main_basefreq(const int argc, const char *argv[])
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
    std::cout << "# name\tread.base\tA\tC\tG\tT\tN" << std::endl;
    for (auto handle : handlesBam) {

        const int readLength = 41;
        std::vector<int> base_A(readLength, 0);
        std::vector<int> base_C(readLength, 0);
        std::vector<int> base_G(readLength, 0);
        std::vector<int> base_T(readLength, 0);
        std::vector<int> base_N(readLength, 0);

        handle->calculateBaseContent(readLength, base_A, base_C, base_G, base_T, base_N);

        for (int i = 0; i < readLength; ++i) {
            std::cout << handle->name() << "\t" << i << "\t" <<
                         base_A.at(i) << "\t" <<
                         base_C.at(i) << "\t" <<
                         base_G.at(i) << "\t" <<
                         base_T.at(i) << "\t" <<
                         base_N.at(i) << "\t" << std::endl;
        }


    }

    return 0;
}
