#include <fstream>
#include "checkFile.h"

bool fileExists(const char* fname) {
    std::ifstream file(fname);
    return file.good();
}
