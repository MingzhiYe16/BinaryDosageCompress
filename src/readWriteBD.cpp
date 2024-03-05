#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <Rcpp.h>
#include "checkFile.h"
#include "readWriteVector.h"

using namespace Rcpp;
using namespace std;


std::vector<float> getFloatsFromTag(bcf_hdr_t *header, bcf1_t *line, const char *tag) {
  int nret = 0;
  float *values = nullptr;
  int n = bcf_get_format_float(header, line, tag, &values, &nret);

  std::vector<float> result;
  if (n > 0) {
    result.assign(values, values + nret);
  } else {
    result.assign(nret, std::numeric_limits<float>::quiet_NaN());
  }

  free(values);
  return result;
}

// [[Rcpp::export]]
void writeVectorsToFile(CharacterVector bd_file_name0, CharacterVector vcf_name0)
{

  if (bd_file_name0.size() == 0 || vcf_name0.size() == 0) {
    stop("CharacterVector is empty.");
  }
  std::string bd_file_name_str = Rcpp::as<std::string>(bd_file_name0[0]);
  const char* bd_file_name = bd_file_name_str.c_str();

  std::string vcf_name_str = Rcpp::as<std::string>(vcf_name0[0]);
  const char* vcf_name = vcf_name_str.c_str();
  if(!fileExists(vcf_name)){
    stop("Input VCF file not exist");
  }

  // const char* fname = "/Users/mingzhiye/CppProjects/htslib/example.VCF";
  htsFile *inf = bcf_open(vcf_name, "r");

  bcf_hdr_t *hdr = bcf_hdr_read(inf);

  // int n=hdr->nhrec;
  // cout<<"asdasdasd"<<n<<"asdasdasd"<<endl;

  bcf1_t *rec = bcf_init();
  int nsamples = bcf_hdr_nsamples(hdr);
  ofstream outFile(bd_file_name, ios::binary);
  if (!outFile)
  {
    cerr << "Error opening file for writing!" << endl;
    return;
  }
  vector<streampos> positions; // Store positions of each vector of strings

  // write sample names
  vector<string> sample_names;
  for (int i = 0; i < nsamples; ++i)
  {
    sample_names.push_back(hdr->samples[i]);
  }
  positions.push_back(outFile.tellp());
  writeStringVector(outFile, sample_names);


  while (bcf_read(inf, hdr, rec) == 0)
  {
    bcf_unpack(rec, BCF_UN_STR);
    bcf_unpack(rec, BCF_UN_INFO);
    vector<string> snp_profile;

    const char *c = bcf_hdr_id2name(hdr, rec->rid);
    string c_s = c;
    snp_profile.push_back(c_s);

    snp_profile.push_back(to_string(rec->pos+1));

    string id_s = rec->d.id;
    snp_profile.push_back(id_s);

    string ref_s = rec->d.allele[0];
    snp_profile.push_back(ref_s);

    // process ALT
    string alt_s;
    if (rec->n_allele > 1)
    {
      alt_s = rec->d.allele[1];
    }
    else
    {
      alt_s = "None";
    }
    string space = " ";

    // printf("%lu\n", (unsigned long)rec->n_allele);
    if (rec->n_allele > 2)
    {
      for (int i = 2; i < rec->n_allele; ++i)
      {
        alt_s = alt_s + space;
        string alt_s2 = rec->d.allele[i];
        alt_s = alt_s + alt_s2;
      }
    }
    snp_profile.push_back(alt_s);

    std::vector<float> gp_values = getFloatsFromTag(hdr, rec, "GP");
    std::vector<float> ds_values = getFloatsFromTag(hdr, rec, "DS");

    std::vector<std::vector<float>> snp_sample(4, std::vector<float>(nsamples, float(std::numeric_limits<float>::quiet_NaN())));

    for (size_t i = 0; i < nsamples; i++) {
      if ( i * 3 + 2 < gp_values.size()) {
        snp_sample[1][i] = gp_values[i * 3 + 0];
        snp_sample[2][i] = gp_values[i * 3 + 1];
        snp_sample[3][i] = gp_values[i * 3 + 2];
      }
      if (i < ds_values.size()) {
        snp_sample[0][i] = ds_values[i];
      }

    }
    // If the length of gp_for_record is less than 3, push_back NA values

    positions.push_back(outFile.tellp());
    writeStringVector(outFile, snp_profile);

    positions.push_back(outFile.tellp());
    writeFloatVectors(outFile, snp_sample);

    
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

// Function to read the positions of vectors of strings from binary file
vector<streampos> readVectorPositions(const string &bd_file_name)
{
  ifstream inFile(bd_file_name, ios::binary);
  if (!inFile)
  {
    Rcpp::stop("Unable to open file");
    return {};
  }

  // Go to the end to get the number of positions
  inFile.seekg(-static_cast<int>(sizeof(size_t)), ios::end);
  size_t numPositions;
  inFile.read(reinterpret_cast<char *>(&numPositions), sizeof(numPositions));

  // Move the file position pointer to the beginning of the positions list
  inFile.seekg(-static_cast<int>(sizeof(size_t) + numPositions * sizeof(streampos)), ios::end);

  // Read positions
  vector<streampos> positions(numPositions);
  for (size_t i = 0; i < numPositions; ++i)
  {
    inFile.read(reinterpret_cast<char *>(&positions[i]), sizeof(positions[i]));
  }

  inFile.close();
  return positions;
}


// [[Rcpp::export]]
Rcpp::NumericVector getStreamPositions(CharacterVector fname0) {
    std::string fname_str = Rcpp::as<std::string>(fname0[0]);
    const char* bd_file_name = fname_str.c_str();
    vector<streampos> positions = readVectorPositions(bd_file_name);
    ifstream inFile(bd_file_name, ios::binary);
    if (!inFile)
    {
      Rcpp::stop("Unable to open file");
      return {};
    }
    vector<streampos> positionsVector = readVectorPositions(bd_file_name);

    std::vector<double> doubleVector;

    for (size_t i = 2; i < positionsVector.size(); i+=2) {
        doubleVector.push_back(static_cast<double>(positionsVector[i]));
    }
    return Rcpp::wrap(doubleVector);
}

// [[Rcpp::export]]
Rcpp::DataFrame readSNP(CharacterVector fname0, int index0)
{
  int index = index0 * 2;
  std::string fname_str = Rcpp::as<std::string>(fname0[0]);
  const char* bd_file_name = fname_str.c_str();
  vector<streampos> positions = readVectorPositions(bd_file_name);
  if(index>=positions.size()){
    cerr << "Index out of range!" << endl;
    return {};
  }
  vector<vector<float>> ds = readFloatVectorsFromFileWithPosition(bd_file_name, positions[index]);
  vector<string> sample_name = readStringVectorFromFileWithPosition(bd_file_name, 0);
  Rcpp::DataFrame df = DataFrame::create(Named("Sample_Names") = sample_name, Named("DS") = ds[0], Named("GP1") = ds[1], Named("GP2") = ds[2], Named("GP3") = ds[3]);
  return df;
}

// [[Rcpp::export]]
Rcpp::DataFrame readSNPwithPosition(CharacterVector fname0, NumericVector pos0){
    std::string fname_str = Rcpp::as<std::string>(fname0[0]);
    const char* bd_file_name = fname_str.c_str();
    if (pos0.size() != 1) {
        Rcpp::stop("Input position invalid, input position length must be one.");
    }
    std::streampos pos = static_cast<std::streampos>(static_cast<std::int64_t>(pos0[0]));
    vector<vector<float>> ds = readFloatVectorsFromFileWithPosition(bd_file_name, pos);
    vector<string> sample_name = readStringVectorFromFileWithPosition(bd_file_name, 0);
    Rcpp::DataFrame df = DataFrame::create(Named("Sample_Names") = sample_name, Named("DS") = ds[0], Named("GP1") = ds[1], Named("GP2") = ds[2], Named("GP3") = ds[3]);
    return df;
}

// [[Rcpp::export]]
Rcpp::CharacterVector readSampleName(CharacterVector fname0){
    std::string fname_str = Rcpp::as<std::string>(fname0[0]);
    const char* bd_file_name = fname_str.c_str();
    return Rcpp::wrap(readStringVectorFromFileWithPosition(bd_file_name, 0));
}

// [[Rcpp::export]]
Rcpp::DataFrame getBDinfo(CharacterVector fname0)
{
  std::string fname_str = Rcpp::as<std::string>(fname0[0]);
  const char* bd_file_name = fname_str.c_str();
  vector<streampos> positions = readVectorPositions(bd_file_name);
  vector<string> chrom;
  vector<int> pos;
  vector<string> id;
  vector<string> ref;
  vector<string> alt;
  for (int i = 1; i < positions.size(); i += 2)
  {
    vector<string> profile = readStringVectorFromFileWithPosition(bd_file_name, positions[i]);
    chrom.push_back(profile[0]);
    pos.push_back(stoi(profile[1]));
    id.push_back(profile[2]);
    ref.push_back(profile[3]);
    alt.push_back(profile[4]);
  }
  Rcpp::DataFrame df = DataFrame::create(Named("CHROM") = chrom, Named("POS") = pos, Named("ID") = id, Named("REF") = ref, Named("ALT") = alt);
  return df;
}

