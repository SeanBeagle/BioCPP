#include<string>
#include<iostream>
#include<fstream>
#include<algorithm>


char complement(char nucleotide);
std::string canonical(std::string fwd);

int main(int argc, char* argv[]) {
  std::string fwd = argv[1];
  std::cout << canonical(fwd) << std::endl;
  return EXIT_SUCCESS;
}

char complement(char nucleotide) {
  switch(nucleotide) {
    case 'a':
    case 'A': nucleotide = 'T';
              break;
    case 't':
    case 'T': nucleotide = 'A';
              break;
    case 'c':
    case 'C': nucleotide = 'G';
              break;
    case 'g':
    case 'G': nucleotide = 'C';
              break;     
  }
  return nucleotide;
}

std::string canonical(std::string fwd) {
    std::string rev = fwd;
    std::for_each(rev.begin(), rev.end(), complement);
    std::reverse(rev.begin(), rev.end());
    std::cout << fwd << " vs " << rev << std::endl;
    return std::min(fwd, rev); 
}
