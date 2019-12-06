#include <iostream>
#include <sstream>
#include <map>
#include <set>

#include "argumentparser.h"
#include "version.h"

#include <htslib/hts.h>
#include <htslib/kseq.h>

struct Exon {
    std::string chrom;
    int32_t chromStart;
    int32_t chromEnd;
    char strand;

    bool operator<(const Exon &rhs) const {
        if (strand != rhs.strand)
            return strand < rhs.strand;
        else if (chrom != rhs.chrom)
            return chrom < rhs.chrom;
        else if (chromStart != rhs.chromStart)
            return chromStart < rhs.chromStart;
        else
            return chromEnd < rhs.chromEnd;
    }
};



struct gtf_t {
    char score = 'x';
    char frame = 'x';
    char strand = 'x';
    int32_t start = 0;
    int32_t end = 0;
    std::string seqname = "";
    std::string source = "";
    std::string feature = "";
    std::string gene_id = "";
    std::string transcript_id = "";

    explicit gtf_t(const char *line = nullptr) {
        if (!line) return;

        std::istringstream iss(line);

        // parse stream
        iss >>
                seqname >>
                source >>
                feature >>
                start >>
                end >>
                score >>
                strand >>
                frame;


        gene_id = parseAttributes(line, "gene_id");
        transcript_id = parseAttributes(line, "transcript_id");
    }


    std::string parseAttributes(const std::string &buffer, const std::string &query) {
        std::size_t posStart = buffer.find(query) + query.size() + 2;
        std::size_t posEnd = buffer.find(";", posStart + 1) - 1;
        if (posStart == std::string::npos)
            return "posStart";

        if (posStart == std::string::npos)
            return "posEnd";



        return buffer.substr(posStart, posEnd - posStart);
    }

};




int main_gtftobed(const int argc, const char *argv[])
{
    std::string fileGtf;

    auto p = ArgumentParser("gtftobed", std::string(VERSION), "converts GTF file format to BED12 file format");
    p.addArgumentRequired("annotation").setKeyShort("-g").setKeyLong("--gtf").setHelp("annotation in GTF file format");

    try {
        p.parse(argc, argv);
        fileGtf = p.get<std::string>("annotation");
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }




    // open GTF file
    htsFile *fp = hts_open(fileGtf.c_str(), "r");
    if (!fp) {
        std::cerr << "ribotools::gtftobed::error, failed to open GTF file to read " << fileGtf << std::endl;
        return EXIT_FAILURE;
    }

    // read GTF file line by line
    std::map<std::string, std::set<Exon>> map_Exons;
    std::map<std::string, std::set<Exon>> map_CDS;

    int counter = 0;
    kstring_t line = {0, 0, nullptr};
    while (hts_getline(fp, KS_SEP_LINE, &line) >= 0) {

        auto record = gtf_t(line.s);
        std::string key = record.transcript_id + ";" + record.gene_id;
        Exon interval = {record.seqname, record.start - 1, record.end, record.strand};

        if (record.feature == "exon")
            map_Exons[key].insert(interval);


        if (record.feature == "CDS")
            map_CDS[key].insert(interval);

        //std::cout << record.feature << "\t" << record.gene_id << "\t" << record.transcript_id << std::endl;
        counter++;
    }

    if (line.s)
        free(line.s);

    std::cerr << "Lines: " << counter << std::endl;


    // close GTF file
    if (fp)
        hts_close(fp);


    // print map
    for (const auto &[name, value] : map_Exons) {

        std::string chrom = value.begin()->chrom;
        int32_t chromStart = value.begin()->chromStart;
        int32_t chromEnd = value.rbegin()->chromEnd;
        int32_t thickStart = chromStart;
        int32_t thickEnd = chromEnd;
        char strand = value.begin()->strand;
        std::string blocksStarts;
        std::string blocksSizes;

        auto itCDS = map_CDS.find(name);
        if (itCDS != map_CDS.end()) {
            thickStart = itCDS->second.begin()->chromStart;
            thickEnd = itCDS->second.rbegin()->chromEnd;
        }


        for (std::set<Exon>::const_iterator it = value.begin();
             it != value.end(); ++it) {

            blocksStarts += std::to_string(it->chromStart - chromStart) + ",";
            blocksSizes += std::to_string(it->chromEnd - it->chromStart) + ",";
            //std::cout << it->chrom << "\t" << it->chromStart << "\t" << it->chromEnd << "\t" << name << "\t" << it->strand << std::endl;
        }

        std::cout << chrom << "\t"
                  << chromStart << "\t"
                  <<  chromEnd << "\t"
                   << name << "\t"
                   << 0 << "\t"
                   << strand << "\t"
                   << thickStart << "\t"
                   << thickEnd << "\t"
                   << 0 << "\t"
                   << value.size() << "\t"
                   << blocksSizes << "\t"
                   << blocksStarts << std::endl;

    }


    return 0;
}

