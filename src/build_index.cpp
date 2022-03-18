#include "btllib/util.hpp"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

int main (int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Wrong args.\n";
    std::exit(EXIT_FAILURE);
  }
  unsigned arg = 1;
  std::ifstream readsfile(argv[arg++]);
  std::ofstream indexfile(argv[arg++]);

  std::string line;
  std::string id;
  long i = 0, byte = 0, id_startbyte = 0, id_endbyte = 0;
  while (std::getline(readsfile, line)) {
    const auto endbyte = byte + line.size();
    if (i % 2 == 0) {
      id_startbyte = byte + 1;
      id_endbyte = endbyte;
      id = btllib::split(line, " ")[0].substr(1);
    } else {
      indexfile << id << '\t' << id_startbyte << '\t' << id_endbyte << '\t' << endbyte << '\n';
    }
    byte = endbyte + 1;
    i++;
  }
  return 0;
}