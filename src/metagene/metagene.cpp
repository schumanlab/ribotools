#include "MetaGene.h"

MetaGene::MetaGene() :
    m_fhBed(nullptr),
    m_fhBam(nullptr),
    m_fhBai(nullptr),
    m_header(nullptr),
    m_bam(nullptr)
{
    
}


MetaGene::~MetaGene()
{
    if (m_fhBed != nullptr)
    {
        if (bgzf_close(m_fhBed) != 0)
            error("failed to close BED file");
    }

    if (m_fhBam != nullptr)
        sam_close(m_fhBam);

    if (m_fhBai != nullptr)
        hts_idx_destroy(m_fhBai);

    if (m_header != nullptr)
        bam_hdr_destroy(m_header);

    if (m_bam != nullptr)
        bam_destroy1(m_bam);
}

void MetaGene::open(const std::string &fileBed, const std::string &fileBam)
{
    // open bed file hande
    m_fhBed = bgzf_open(fileBed.c_str(), "r");
    if (!m_fhBed)
    {
        error("failed to open BED file" + fileBed);
        return;
    }

    // open bam file handle
    m_fhBam = hts_open(fileBam.c_str(), "r");
    if (m_fhBam == nullptr)
    {
        error("failed to open BAM file " + fileBam + " to read");
        return;
    }

    // read bam header
    m_header = sam_hdr_read(m_fhBam);
    if (m_header == nullptr)
    {
        error("failed to read BAM header");
        return;
    }

    // open bam index
    std::string fileBai = fileBam + ".bai";
    m_fhBai = sam_index_load(m_fhBam, fileBai.c_str());
    if (m_fhBai == nullptr)
    {
        error("failed to read BAM index" + fileBai);
        return;
    }

    // initialise bam record
    m_bam = bam_init1();
}


void MetaGene::pileup()
{
    kstring_t line = {0, 0, NULL};
    uint32_t bedLineCounter = 0;
    while (bgzf_getline(m_fhBed, '\n', &line) > 0)
    {
        // parse bed line
        auto ssline = std::stringstream(line.s);
        auto bed = BedLine();
        ssline >> bed;

        // parse exons
        bed.parseExons();
        if (bed.cdsSpan == 0)
            continue;

        // get reads per gene
        std::string region = bed.chrom + ':' + std::to_string(bed.thickStart) + '-' + std::to_string(bed.thickEnd);
        hts_itr_t *iterBam;
        if ((iterBam = bam_itr_querys(m_fhBai, m_header, region.c_str())) == 0)
            continue;
        
        // accumulate coverage
        khash_t(32) *hmap = kh_init(32);
        khiter_t hiter;
        int absent;
        int reads = 0;
        while (bam_itr_next(m_fhBam, iterBam, m_bam) >= 0)
        {
            int32_t readStart = m_bam->core.pos;
            //int32_t readLength = bam_calqlen(&m_bam->core, bam_get_cigar(m_bam));
            int32_t readOffset = 0;
            bed.toLinear(readOffset, readStart);
            int32_t readRelative = (bed.strand == '+') ? (readOffset - bed.cdsStart) : (bed.cdsEnd - readOffset);

            // accumulate coverage
            hiter = kh_put(32, hmap, readRelative, &absent);
            if (absent)
            {
                kh_value(hmap, hiter) = 1;
            }
            else
            {
                kh_value(hmap, hiter)++;
            }
            reads++;
        }
        hts_itr_destroy(iterBam);
        
        if (reads < (bed.cdsSpan / 30))
            continue;

        // accumulate result
        std::cout << bed.name << '\t' << bed.cdsSpan << '\t';
        for (hiter = kh_begin(hmap); hiter != kh_end(hmap); ++hiter)
        {
            if (kh_exist(hmap, hiter))
            {
                int32_t key = kh_key(hmap, hiter);
                uint32_t value = kh_value(hmap, hiter);
                std::cout << key << ',' << value << ',';
            }
        }
        std::cout << '\n';
        kh_destroy(32, hmap);

        
        bedLineCounter++;
    }
    free(ks_release(&line));
    std::cerr << "BedLines: " << bedLineCounter << std::endl;
}


void MetaGene::error(const std::string &errorMessage)
{
    std::cerr << "MetaGene::Error" << std::endl;
    std::cerr << '\t' << errorMessage << std::endl;
}

/* bam_calqlen
 * Calculate query alignment length
 */
inline uint32_t MetaGene::bam_calqlen(const bam1_core_t *core, const uint32_t *cigar)
{
    int i;
    uint32_t qlen = 0;
    
    for (i = 0; i < core->n_cigar; ++i)
    {
        int span = cigar[i] >> 4;
        int op = cigar[i] & 0xf;
        
        if ( (op == BAM_CMATCH) || (op == BAM_CDEL) || (op == BAM_CREF_SKIP))
        {
            qlen += span;
        }
    }
    
    return qlen;
}