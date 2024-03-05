//g++ -g /Users/mingzhiye/CppProjects/htslib/src/main.cpp -o /Users/mingzhiye/CppProjects/htslib/src/main -I/Users/mingzhiye/htslib/include -L/Users/mingzhiye/htslib/lib -lhts -lm
// cp /Users/mingzhiye/CppProjects/htslib/src/Makevars ~/.R/Makevars
#include<iostream>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <Rcpp.h>
#include <fstream>
#include "checkFile.h"

using namespace Rcpp;
using namespace std;

#define BCF_DT_ID       0


// [[Rcpp::export]]
Rcpp::DataFrame readVCF(CharacterVector fname0)
{
    if (fname0.size() == 0) {
        stop("CharacterVector is empty.");
    }
    std::string fname_str = Rcpp::as<std::string>(fname0[0]);
    const char* fname = fname_str.c_str();

    if(!fileExists(fname)){
        stop("Input VCF file not exist");
    }

    //const char* fname = Rcpp::as<const char*>(fname0);
    htsFile *fp    = hts_open(fname,"rb");

    bcf_hdr_t *hdr = bcf_hdr_read(fp);

    //int n=hdr->nhrec;
    //cout<<"asdasdasd"<<n<<"asdasdasd"<<endl;

    bcf1_t *rec    = bcf_init();
    vector<string> chrom;
    vector<int> pos;
    vector<string> id;
    vector<string> ref;
    vector<string> alt;
    vector<float> qual;
    vector<string> filter;
    while ( bcf_read(fp, hdr, rec)>=0 )
    {
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);


        //Process CHROM
        const char* c=bcf_hdr_id2name(hdr, rec->rid);
        string c_s=c;
        chrom.push_back(c_s);


        //Process POS
        pos.push_back(rec->pos+1);

        //Process QUAL
        qual.push_back(rec->qual);

        //Process ID
        string id_s=rec->d.id;
        id.push_back(id_s);

        //process REF
        string ref_s=rec->d.allele[0];
        ref.push_back(ref_s);

        //process ALT
        string alt_s;
        if(rec->n_allele > 1){alt_s=rec->d.allele[1];}
        else{alt_s="None";}
        string space=" ";

        //printf("%lu\n", (unsigned long)rec->n_allele);
        if(rec->n_allele > 2){
            for (int i=2; i<rec->n_allele; ++i){
                alt_s=alt_s+space;
                string alt_s2=rec->d.allele[i];
                alt_s=alt_s+alt_s2;
                }
        }
        alt.push_back(alt_s);
        //check if is snp
        //printf("ttt %d\n", bcf_is_snp(rec));

        //Process FILTER
        string filter_s;
        if (rec->d.n_flt<=0){filter_s=".";}
        else{filter_s="";}
        for(int i=0; i<rec->d.n_flt; ++i){
            if(i!=0){filter_s=filter_s+space;}
            int filter_i_id=rec->d.flt[i];
            string filter_i= bcf_hdr_int2id(hdr, BCF_DT_ID, filter_i_id);
            filter_s=filter_s+filter_i;

        }
        filter.push_back(filter_s);
    }



    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
    //Rcpp::DataFrame df = DataFrame::create(Named("CHROM") = chrom,Named("POS") = pos,Named("REF") = ref,Named("ALT") = alt,Named("QUAL") = qual,Named("FILTER") = filter);
    Rcpp::DataFrame df = DataFrame::create(Named("CHROM") = chrom, Named("POS") = pos,Named("ID") = id,Named("QUAL") = qual,Named("REF") = ref,Named("ALT") = alt,Named("FILTER") = filter);
    return df;
}
