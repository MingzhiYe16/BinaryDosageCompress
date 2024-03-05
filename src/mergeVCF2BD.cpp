#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <cmath> // for NAN
#include <Rcpp.h>
#include "readWriteVector.h"

using namespace Rcpp;
using namespace std;

typedef std::vector<float> GP_Vector; // Vector to hold GP values
typedef std::unordered_map<std::string, std::unordered_map<std::string, GP_Vector>> SNP_Matrix; // Matrix with SNP IDs as rows and sample names as columns
typedef std::vector<std::string> SNP_Info_Vector; // Vector to hold CHROM, POS, ID, REF, ALT
typedef unordered_set<std::string> StringSet; // Set to hold unique strings

// Helper function to parse genotype information for a SNP
void parseGenotypes(bcf_hdr_t *hdr, bcf1_t *rec, const std::string& snp_id, SNP_Matrix& snpMatrix) {
    int nsmpl = bcf_hdr_nsamples(hdr);

    float *ds_arr = nullptr, *gp_arr = nullptr;
    int nds = 0, ngp = 0;

    // Try to get DS values
    nds = bcf_get_format_float(hdr, rec, "DS", &ds_arr, &nds);
    // Try to get GP values
    ngp = bcf_get_format_float(hdr, rec, "GP", &gp_arr, &ngp);

    // Loop over all samples
    for (int i = 0; i < nsmpl; i++) {
        // Initialize with NAN
        GP_Vector genotypeData(4, NAN);

        // DS
        if (nds > 0 && i < nds) {
            genotypeData[0] = ds_arr[i];
        }

        // GP
        if (ngp > 0) {
            for (int j = 0; j < 3 && 3*i+j < ngp; j++) { // Assuming GP has three values per sample
                genotypeData[1 + j] = gp_arr[3*i+j];
            }
        }

        // Insert into matrix
        std::string sampleName = hdr->samples[i];
        snpMatrix[snp_id][sampleName] = genotypeData;
    }

    // Clean up
    if (gp_arr) free(gp_arr);
}


// Function to read all VCFs and transform into a matrix and a std::vector of vectors of strings
// [[Rcpp::export]]
void mergeVCF2BD(Rcpp::CharacterVector vcfFiles, Rcpp::CharacterVector bd_file_name0) {
    std::string bd_file_name_str = Rcpp::as<std::string>(bd_file_name0[0]);
    const char* bd_file_name = bd_file_name_str.c_str();
    SNP_Matrix snpMatrix;
    std::vector<SNP_Info_Vector> snpInfo;
    std::vector<std::string> snpVector;

    ofstream outFile(bd_file_name, ios::binary);
    if (!outFile)
    {
        cerr << "Error opening file for writing!" << endl;
        return;
    }
    StringSet snpSet; // To hold unique SNPs
    StringSet sampleSet; // To hold unique sample names

    // Parse each VCF file
    for (int file_index = 0; file_index < vcfFiles.size(); ++file_index) {
        std::string vcfFile = Rcpp::as<std::string>(vcfFiles[file_index]);
        htsFile *vcf_fp = bcf_open(vcfFile.c_str(), "r");
        bcf_hdr_t *bcf_hdr = bcf_hdr_read(vcf_fp);
        bcf1_t *bcf_rec = bcf_init();

        int numSamples = bcf_hdr_nsamples(bcf_hdr);
        for (int i = 0; i < numSamples; ++i) {
            std::string sample_name(bcf_hdr_int2id(bcf_hdr, BCF_DT_SAMPLE, i));
            sampleSet.insert(sample_name); // Insert sample name into sampleSet
        }

        while (bcf_read(vcf_fp, bcf_hdr, bcf_rec) == 0) {
            if (bcf_unpack(bcf_rec, BCF_UN_ALL) == 0) {
                // Construct the SNP identifier
                std::string chrom = bcf_seqname(bcf_hdr, bcf_rec);
                std::string pos = to_string(bcf_rec->pos+1);
                std::string ref = bcf_rec->d.allele[0];
                std::string alt;
                if (bcf_rec->n_allele > 1)
                {
                    alt = bcf_rec->d.allele[1];
                }
                else
                {
                    alt = "None";
                }
                std::string space = " ";

                // printf("%lu\n", (unsigned long)rec->n_allele);
                if (bcf_rec->n_allele > 2)
                {
                    for (int i = 2; i < bcf_rec->n_allele; ++i)
                    {
                        alt = alt + space;
                        std::string alt2 = bcf_rec->d.allele[i];
                        alt = alt + alt2;
                    }
                }
                std::string snp_id = chrom + ":" + pos + ":" + ref + ":" + alt;

                // Insert SNP info if it's new
                if (snpSet.insert(snp_id).second) {
                    snpInfo.push_back({chrom, pos, bcf_rec->d.id, ref, alt});
                    snpVector.push_back(snp_id);
                }

                // Parse genotypes for this SNP
                parseGenotypes(bcf_hdr, bcf_rec, snp_id, snpMatrix);

 
            }
        }

        // Clean up
        bcf_destroy(bcf_rec);
        bcf_hdr_destroy(bcf_hdr);
        bcf_close(vcf_fp);
    }

    std::vector<std::string> sampleVector;
    sampleVector.reserve(sampleSet.size());
    for (const auto& sampleName : sampleSet) {
        sampleVector.push_back(sampleName);
    }


    
    // Fill in missing data for each SNP and sample
    for (const auto &snp_id : snpVector) {
        for (const auto &sampleName : sampleVector) {
            // If the sample does not have this SNP, insert a std::vector with NANs
            if (snpMatrix[snp_id].find(sampleName) == snpMatrix[snp_id].end()) {
                snpMatrix[snp_id][sampleName] = GP_Vector(4, NAN);
            }
        }
    }
    cout<<"read completed"<<endl;

    std::vector<streampos> positions; // Store positions of each std::vector of strings

    positions.push_back(outFile.tellp());
    writeStringVector(outFile, sampleVector);


    int index = 0;
    for (const auto &snp_id : snpVector) 
    {
        positions.push_back(outFile.tellp());
        writeStringVector(outFile, snpInfo[index]);
        index++;

        // get data_in_samples 2D std::vector from snpMatrix[snp_id] for write into BD file
        std::vector<std::vector<float>> data_in_samples(4, std::vector<float>(snpMatrix[snp_id].size(), float(std::numeric_limits<float>::quiet_NaN())));
        int j=0;
        for (const auto &sample : sampleVector){
            for(int i= 0;i<4;i++){
                data_in_samples[i][j]=snpMatrix[snp_id][sample][i];
            }
            j++;
        }

        // write sample data into BD file
        positions.push_back(outFile.tellp());
        writeFloatVectors(outFile, data_in_samples);
    }
    // Now write positions at the end of the file
    size_t numPositions = positions.size();

    for (const auto &pos : positions)
    {
        outFile.write(reinterpret_cast<const char *>(&pos), sizeof(pos));
    }
    outFile.write(reinterpret_cast<const char *>(&numPositions), sizeof(numPositions));

    outFile.close();
}

