#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <htslib/sam.h>

struct Region {int64_t start,end; std::string id; char strand; uint64_t count=0;};
int main(int argc,char** argv) {
 try {
  if(argc!=2) throw std::runtime_error("Usage: count_eclip regions.bed < BAM");
  std::ifstream input(argv[1]);
  std::map<std::string,std::vector<Region>> regions;
  std::string chr,id; int64_t start,end; int score; char strand;
  while(input>>chr>>start>>end>>id>>score>>strand) {
   if(start<0 || end<=start || (strand!='+' && strand!='-')) throw std::runtime_error("Invalid BED6 region");
   regions[chr].push_back({start,end,id,strand,0});
  }
  if(!input.eof() || regions.empty()) throw std::runtime_error("Missing/invalid regions");
  auto bam=std::unique_ptr<samFile,decltype(&hts_close)>(sam_open("-","rb"),hts_close);
  if(!bam) throw std::runtime_error("Cannot open BAM");
  auto header=std::unique_ptr<sam_hdr_t,decltype(&sam_hdr_destroy)>(sam_hdr_read(bam.get()),sam_hdr_destroy);
  if(!header) throw std::runtime_error("Missing BAM header");
  const int chr1=sam_hdr_name2tid(header.get(),"chr1");
  if(chr1<0 || sam_hdr_tid2len(header.get(),chr1)!=248956422) throw std::runtime_error("Expected GRCh38 coordinate system");
  std::vector<std::vector<Region>*> indexed(sam_hdr_nref(header.get()),nullptr);
  for(auto& entry:regions) {
   auto& rs=entry.second;
   std::sort(rs.begin(),rs.end(),[](const Region&a,const Region&b){return a.start<b.start;});
   int tid=sam_hdr_name2tid(header.get(),entry.first.c_str());
   if(tid<0) throw std::runtime_error("Region contig absent");
   for(const auto&r:rs) if(r.end>sam_hdr_tid2len(header.get(),tid)) throw std::runtime_error("Region exceeds contig");
   indexed[tid]=&rs;
  }
  auto rec=std::unique_ptr<bam1_t,decltype(&bam_destroy1)>(bam_init1(),bam_destroy1);
  uint64_t records=0,usable=0,r1=0,r2=0,single=0,filtered=0,multi=0,nh_missing=0,spliced=0;
  const uint16_t rejected=BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP|BAM_FSUPPLEMENTARY;
  int status;
  while((status=sam_read1(bam.get(),header.get(),rec.get()))>=0) {
   ++records; const auto& c=rec->core;
   if((c.flag&rejected) || c.tid<0) {++filtered;continue;}
   if(c.flag&BAM_FPAIRED) {
    if(c.flag&BAM_FREAD1) {++r1;continue;}
    if(!(c.flag&BAM_FREAD2)) throw std::runtime_error("Paired alignment has no read-number flag");
   }
   auto* nh=bam_aux_get(rec.get(),"NH");
   if(nh && bam_aux2i(nh)>1) {++multi;continue;}
   if(!nh) ++nh_missing;
   ++usable; if(c.flag&BAM_FPAIRED) ++r2; else ++single;
   auto* rs=indexed.at(c.tid); if(!rs) continue;
   int64_t left=c.pos,right=bam_endpos(rec.get());
   auto stop=std::lower_bound(rs->begin(),rs->end(),right,[](const Region&r,int64_t v){return r.start<v;});
   if(stop==rs->begin()) continue;
   std::vector<std::pair<int64_t,int64_t>> blocks;
   int64_t position=left; auto* cigar=bam_get_cigar(rec.get()); bool splice=false;
   for(uint32_t j=0;j<c.n_cigar;++j) {
    uint32_t op=bam_cigar_op(cigar[j]),len=bam_cigar_oplen(cigar[j]);
    if(op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF) blocks.push_back({position,position+len});
    if(op==BAM_CREF_SKIP) splice=true;
    if(bam_cigar_type(op)&2) position+=len;
   }
   if(splice) ++spliced;
   char read_strand=(c.flag&BAM_FREVERSE)?'-':'+';
   for(auto it=rs->begin();it!=stop;++it) {
    if(it->end<=left || it->strand!=read_strand) continue;
    for(const auto& block:blocks) if(block.first<it->end && block.second>it->start) {++it->count;break;}
   }
  }
  if(status!=-1 || usable==0) throw std::runtime_error("Truncated BAM or no usable read-2/single alignments");
  std::cout<<"#meta\t{\"records\":"<<records<<",\"usable\":"<<usable<<",\"read1_skipped\":"<<r1
   <<",\"read2_used\":"<<r2<<",\"unpaired_used\":"<<single<<",\"filtered\":"<<filtered
   <<",\"multimappers_skipped\":"<<multi<<",\"NH_missing\":"<<nh_missing<<",\"spliced_region_contig_reads\":"<<spliced<<"}\n";
  for(const auto& entry:regions) for(const auto&r:entry.second) std::cout<<r.id<<'\t'<<r.count<<'\n';
  return 0;
 } catch(const std::exception&e) {std::cerr<<e.what()<<'\n';return 1;}
}
