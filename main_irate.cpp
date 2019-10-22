#include <iostream>
#include <memory>
#include <fstream>
#include <unordered_map>

#include "argumentparser.h"
#include "bedrecord.h"
#include "bamhandle.h"
#include "version.h"

void loadAndCountBamHandles(std::vector<std::shared_ptr<BamHandle>> &handlesBam, const std::vector<std::string> &filesBam)
{
    for (const auto &fileName : filesBam) {
        std::shared_ptr<BamHandle> handle = std::make_shared<BamHandle>(fileName, 255, 0);
        handle->countUniqueReads();
        std::cerr << handle->name() << "\t" << handle->reads() << std::endl;
        handlesBam.emplace_back(handle);
    }
}

int main_irate(const int argc, const char *argv[])
{
    std::string fileBed;
    std::vector<std::string> filesBam_rfp;
    std::vector<std::string> filesBam_rna;
    std::vector<std::shared_ptr<BamHandle>> handlesBam_rfp;
    std::vector<std::shared_ptr<BamHandle>> handlesBam_rna;

    auto p = ArgumentParser("irate", std::string(VERSION), "calculate initiation rate from Ribo Footprint and RNA coverage");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentRequired("RFP").setKeyShort("-f").setKeyLong("--rfp").setHelp("BAM files of ribosome footrpint coverage").setCount(-1);
    p.addArgumentRequired("RNA").setKeyShort("-n").setKeyLong("--rna").setHelp("BAM files of RNA coverage").setCount(-1);
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        filesBam_rfp = p.get<std::vector<std::string>>("RFP");
        filesBam_rna = p.get<std::vector<std::string>>("RNA");
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // open BAM handles
    loadAndCountBamHandles(handlesBam_rfp, filesBam_rfp);
    loadAndCountBamHandles(handlesBam_rna, filesBam_rna);

    // check if handles size is matching
    if (handlesBam_rfp.size() != handlesBam_rna.size()) {
        std::cerr << "ribotools::irate::error, provide equal number of BAM files for Ribosome Footprint(--rfp) and RNASeq(--rna)" << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::irate::rrror, failed to open BED reference " << fileBed << std::endl;
        return 1;
    }

    // loop over bed records
    std::string line;
    while (std::getline(fhBed, line)) {

        // parse bed line
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        // stats per file
        std::vector<std::shared_ptr<BamHandle>>::iterator iteratorRFP;
        std::vector<std::shared_ptr<BamHandle>>::iterator iteratorRNA;

        std::cout << bed.gene << "\t" << bed.span << "\t" << bed.cdsSpan;

        for (iteratorRFP = handlesBam_rfp.begin(), iteratorRNA = handlesBam_rna.begin();
             (iteratorRFP != handlesBam_rfp.end()) && (iteratorRNA != handlesBam_rna.end());
             ++iteratorRFP, ++iteratorRNA) {

            std::vector<int> rfpc(static_cast<size_t>(bed.cdsSpan), 0);
            int reads_rna = (*iteratorRNA)->readsPerRegion(bed.transcript, 0, bed.span);

            (*iteratorRFP)->calculateSiteCoverage(rfpc, bed.transcript, bed.cdsStart, bed.cdsEnd, bed.cdsStart, true);
            int reads_rfp = 0;
            int reads_ini = 0;
            for (size_t k = 3; k < rfpc.size(); ++k) {

                reads_rfp += rfpc.at(k);

                if (k < 33)
                    reads_ini += rfpc.at(k);

            }

            std::cout << "\t" << (*iteratorRFP)->reads() << "\t" << reads_rfp << "\t" << reads_ini << "\t" <<
                         (*iteratorRNA)->reads() << "\t" << reads_rna;
        }
        std::cout << std::endl;


    }


    // close bed file
    fhBed.close();

    return 0;
}
