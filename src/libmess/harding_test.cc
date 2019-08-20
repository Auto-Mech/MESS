#include "harding.hh"
#include <iostream>
#include <fstream>

int main (int argc, char* argv[])
{
  if(argc != 2) {
    std::cerr << "usage: test input_file\n";
    return 1;
  }

  std::ifstream from(argv[1]);
  //
  if(!from) {
    std::cerr << "harding_test: cannot open " << argv[1] << " file\n";
    return 1;
  }

  harding_init_(from);
  return 0;
}
