#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <htslib/sam.h>

struct Region {
    int64_t start, end;
    std::string id;
};
struct Chromosome {
    std::vector<Region> regions;
    std::vector<int64_t> starts, ends;
    std::vector<int64_t> left_bins, right_bins;

    void prepare() {
        for (const auto &region : regions) {
            starts.push_back(region.start);
            ends.push_back(region.end);
        }
        std::sort(starts.begin(), starts.end());
        std::sort(ends.begin(), ends.end());
        starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
        ends.erase(std::unique(ends.begin(), ends.end()), ends.end());
        left_bins.resize(ends.size() + 1);
        right_bins.resize(starts.size() + 1);
    }

    void count(int64_t left, int64_t right) {
        if (left >= ends.back() || right <= starts.front()) return;
        ++left_bins[std::upper_bound(ends.begin(), ends.end(), left) - ends.begin()];
        ++right_bins[std::lower_bound(starts.begin(), starts.end(), right) - starts.begin()];
    }

    void report() {
        for (size_t i = 1; i < left_bins.size(); ++i) left_bins[i] += left_bins[i - 1];
        for (size_t i = 1; i < right_bins.size(); ++i) right_bins[i] += right_bins[i - 1];

        for (const auto &region : regions) {
            const size_t end = std::lower_bound(ends.begin(), ends.end(), region.end) - ends.begin();
            const size_t start = std::lower_bound(starts.begin(), starts.end(), region.start) - starts.begin();

            const auto count = left_bins[end] - right_bins[start];
            if (count < 0) throw std::runtime_error("Negative interval count");
            std::cout << region.id << '\t' << count << '\n';
        }
    }
};

int run(const char *bed_path) {
    std::ifstream bed(bed_path);
    if (!bed) throw std::runtime_error("Cannot read BED");
    std::map<std::string, Chromosome> regions;
    std::string chrom, id;
    int64_t start, end;
    while (bed >> chrom >> start >> end >> id) {
        if (start < 0 || end <= start) throw std::runtime_error("Invalid BED interval");
        regions[chrom].regions.push_back({start, end, id});
    }
    if (!bed.eof() || regions.empty()) throw std::runtime_error("Invalid or empty BED");

    auto bam = std::unique_ptr<samFile, decltype(&hts_close)>(sam_open("-", "rb"), hts_close);
    if (!bam) throw std::runtime_error("Cannot open BAM stream");
    auto header = std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)>(sam_hdr_read(bam.get()), sam_hdr_destroy);
    if (!header) throw std::runtime_error("Cannot read BAM header");
    if (sam_hdr_name2tid(header.get(), "chr1") < 0 ||
        sam_hdr_tid2len(header.get(), sam_hdr_name2tid(header.get(), "chr1")) != 248956422) {
        throw std::runtime_error("Unexpected BAM coordinate system; expected GRCh38 chr1");
    }
    std::vector<Chromosome *> indexed(sam_hdr_nref(header.get()), nullptr);
    for (auto &[name, group] : regions) {
        const int tid = sam_hdr_name2tid(header.get(), name.c_str());
        if (tid < 0) throw std::runtime_error("Region contig missing: " + name);
        group.prepare();
        if (group.ends.back() > sam_hdr_tid2len(header.get(), tid)) throw std::runtime_error("Region exceeds contig");
        indexed[tid] = &group;
    }

    auto record = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>(bam_init1(), bam_destroy1);
    if (!record) throw std::runtime_error("Cannot allocate BAM record");
    constexpr uint16_t rejected = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY;
    int status;
    uint64_t paired = 0, single = 0;
    while ((status = sam_read1(bam.get(), header.get(), record.get())) >= 0) {
        const auto &core = record->core;
        if ((core.flag & rejected) || core.tid < 0) continue;
        int64_t left = core.pos, right;
        if (core.flag & BAM_FPAIRED) {
            if (!(core.flag & BAM_FREAD1) || !(core.flag & BAM_FPROPER_PAIR) ||
                (core.flag & BAM_FMUNMAP) || core.mtid != core.tid || core.mpos < 0 || !core.isize) continue;
            left = std::min<int64_t>(core.pos, core.mpos);
            right = left + (core.isize > 0 ? core.isize : -core.isize);
            ++paired;
        } else {
            right = bam_endpos(record.get());
            ++single;
        }
        if (right <= left) throw std::runtime_error("Invalid mapped fragment");
        if (auto *group = indexed.at(core.tid)) group->count(left, right);
    }
    if (status != -1 || paired + single == 0 || (paired && single)) {
        throw std::runtime_error("Truncated BAM, zero usable data, or mixed read layouts");
    }
    std::cout << "#usable\t" << paired + single << '\t' << (paired ? "paired_fragment" : "single_read") << '\n';
    for (auto &[name, group] : regions) group.report();
    return 0;
}

int main(int argc, char **argv) {
    try {
        if (argc != 2) throw std::runtime_error("Usage: count_fragments regions.bed < input.bam");
        return run(argv[1]);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
