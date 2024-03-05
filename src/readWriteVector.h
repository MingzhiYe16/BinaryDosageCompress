#ifndef READ_WRITE_VECTOR_H
#define READ_WRITE_VECTOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void writeStringVector(std::ofstream& outFile, std::vector<std::string> string_vector);

void writeFloatVectors(std::ofstream& outFile, std::vector<std::vector<float>>  float_vectors);

std::vector<std::vector<float>> readFloatVectorsFromFileWithPosition(const std::string& filename, std::streampos position);

std::vector<std::string> readStringVectorFromFileWithPosition(const std::string &bd_file_name, std::streampos position);

#endif // READ_WRITE_VECTOR_H
