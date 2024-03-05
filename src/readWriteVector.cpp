#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "readWriteVector.h"
#include <zlib.h>
#include <cstdint>
using namespace std;

uint16_t BIT_FOR_VALUE = 10;
uint16_t VALUE_FOR_COMPLIMENT = 3;
uint16_t VALUE_FOR_REMOVE_BITS = 0x03FF; // 6 leftest bits be 0 and other 10(BIT_FOR_VALUE) bits be 1

std::vector<uint16_t> flatten(const std::vector<std::vector<float>>& vecOfVec) {
    size_t subVecSize = vecOfVec.empty() ? 0 : vecOfVec[0].size(); //vecOfVec is 4*subVecSize
    std::vector<uint16_t> flatVec;
    flatVec.reserve(subVecSize * 4); // Pre-allocate memory for efficiency

    size_t subvector_index=0;
    size_t entry_index=0;
    for (const auto& subVec : vecOfVec) {
        for (float value : subVec) {
            // Convert each float to a fixed-point representation and then to uint16_t
            int fixedPoint = static_cast<int>(value * 1000);
            uint16_t packedValue = static_cast<uint16_t>(fixedPoint);

            flatVec.push_back(packedValue);

            if(subvector_index == 3){
                uint16_t DS_compliment = ((flatVec[entry_index] - (flatVec[entry_index + 2 * subVecSize] + 2 * flatVec[entry_index + 3 * subVecSize])) + VALUE_FOR_COMPLIMENT) << BIT_FOR_VALUE; // for DS
                uint16_t GP0_compliment = (flatVec[entry_index + subVecSize] - (1000 - (flatVec[entry_index + 2 * subVecSize] + flatVec[entry_index + 3 * subVecSize])) + VALUE_FOR_COMPLIMENT) << BIT_FOR_VALUE; // for GP0
                flatVec[entry_index + 2 * subVecSize] |= DS_compliment;
                flatVec[entry_index + 3 * subVecSize] |= GP0_compliment;
            }

            entry_index++;
        }
        subvector_index++;
        entry_index=0;
    }
    if (subVecSize * 2 < flatVec.size()) {
        // Erase elements from the beginning up to the subVecSize-th element.
        // Note: begin() + subVecSize is the position of the subVecSize-th element (0-indexed).
        flatVec.erase(flatVec.begin(), flatVec.begin() + subVecSize * 2);
    }

    return flatVec;
}

std::vector<Bytef> compressData(const std::vector<uint16_t>& data) {
    uLongf compressedSize = compressBound(data.size() * sizeof(uint16_t));
    std::vector<Bytef> compressedData(compressedSize);

    compress(compressedData.data(), &compressedSize, 
             reinterpret_cast<const Bytef*>(data.data()), 
             data.size() * sizeof(uint16_t));

    compressedData.resize(compressedSize);
    return compressedData;
}

std::vector<std::vector<float>> reconstruct(const std::vector<uint16_t>& flatVec, size_t n) {
    if (flatVec.size() < 4 * n) {
        throw std::invalid_argument("flatVec does not contain enough elements.");
    }

    std::vector<std::vector<uint16_t>> vecOfVec(4, std::vector<uint16_t>(n, 0));

    for (size_t i = 0; i < n; ++i) {
        uint16_t DS_compliment = flatVec[i] >> BIT_FOR_VALUE; // for DS
        uint16_t GP0_compliment = flatVec[i + n] >> BIT_FOR_VALUE; // for GP0

        vecOfVec[2][i] = flatVec[i] & VALUE_FOR_REMOVE_BITS;
        vecOfVec[3][i] = flatVec[i + n]  & VALUE_FOR_REMOVE_BITS;
        vecOfVec[1][i] = std::max(vecOfVec[2][i] + vecOfVec[3][i] + VALUE_FOR_COMPLIMENT, 1000 + GP0_compliment) - vecOfVec[2][i] - vecOfVec[3][i] - VALUE_FOR_COMPLIMENT; //GP0
        vecOfVec[0][i] = std::max(int(VALUE_FOR_COMPLIMENT), vecOfVec[2][i] + 2 * vecOfVec[3][i] + DS_compliment) - VALUE_FOR_COMPLIMENT; //DS
    }

    std::vector<std::vector<float>> result(4, std::vector<float>(n));
    for (size_t i = 0; i < 4; ++i) {
        std::transform(vecOfVec[i].begin(), vecOfVec[i].end(), result[i].begin(),
                       [](uint16_t val) { return static_cast<float>(val) / 1000.0f; });
    }

    return result;
}


std::vector<uint16_t> decompressData(const std::vector<Bytef>& compressedData, size_t originalSize) {
    std::vector<uint16_t> decompressedData(originalSize / sizeof(uint16_t));
    uLongf decompressedSize = decompressedData.size() * sizeof(uint16_t);

    uncompress(reinterpret_cast<Bytef*>(decompressedData.data()), &decompressedSize,
               compressedData.data(), compressedData.size());

    return decompressedData;
}

void writeFloatVectors(std::ofstream& outFile, std::vector<std::vector<float>>  float_vectors) {
    std::vector<Bytef> compressedData=compressData(flatten(float_vectors));

    size_t dataSize = compressedData.size();
    outFile.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));

    size_t n = float_vectors[0].size();
    outFile.write(reinterpret_cast<const char*>(&n), sizeof(n));

    // Write the actual compressed data
    outFile.write(reinterpret_cast<const char*>(compressedData.data()), dataSize);
}

void writeStringVector(std::ofstream& outFile, std::vector<std::string> string_vector){
    size_t numStrings0 = string_vector.size();
    outFile.write(reinterpret_cast<const char *>(&numStrings0), sizeof(numStrings0)); // Write number of strings

    for (const auto &str : string_vector)
    {
        size_t length0 = str.length();
        outFile.write(reinterpret_cast<const char *>(&length0), sizeof(length0)); // Write string length
        outFile.write(str.c_str(), length0);                                      // Write string
  }
}

std::vector<std::vector<float>> readFloatVectorsFromFileWithPosition(const std::string& filename, std::streampos position) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
      cerr << "Error opening file for reading!" << endl;
      return {}; 
    }

    file.seekg(position);

    // Read the size of the compressed data
    size_t dataSize;
    file.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));

    // Read the length of the sub-vectors
    size_t n;
    file.read(reinterpret_cast<char*>(&n), sizeof(n));

    // Read the compressed data
    std::vector<Bytef> compressedData(dataSize);
    file.read(reinterpret_cast<char*>(compressedData.data()), dataSize);

    // Decompress the data
    size_t originalSize = n * 2 * sizeof(float); // Calculate the original size
    std::vector<uint16_t> decompressedData = decompressData(compressedData, originalSize);

    // Reconstruct the original vector of vectors
    return reconstruct(decompressedData, n);
}


// Function to read a vector of strings from the binary file using its position
std::vector<std::string> readStringVectorFromFileWithPosition(const std::string &bd_file_name, std::streampos position)
{
  ifstream inFile(bd_file_name, ios::binary);
  if (!inFile)
  {
    cerr << "Error opening file for reading!" << endl;
    return {};
  }

  inFile.seekg(position);

  // Read number of strings in the vector
  size_t numStrings;
  inFile.read(reinterpret_cast<char *>(&numStrings), sizeof(numStrings));

  // Read strings
  std::vector<std::string> result;
  for (size_t i = 0; i < numStrings; ++i)
  {
    size_t length;
    inFile.read(reinterpret_cast<char *>(&length), sizeof(length));

    char *buffer = new char[length];
    inFile.read(buffer, length);
    result.emplace_back(buffer, length);
    delete[] buffer;
  }

  inFile.close();
  return result;
}