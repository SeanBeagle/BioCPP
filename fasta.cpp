#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// FORWARD DECLARATION
struct Fasta;
struct Record;


// HEADERS
struct Fasta {
  Fasta(std::string input);
  std::string input;
  std::vector<unsigned long> residues = std::vector<unsigned long>(256, 0);
  unsigned long size = 0;
  unsigned contigs = 0;
  void print();
  static bool isHeader(std::string &line);
};

struct Record {
  Record();
  Record(std::string header);
  std::string header = "No Header";
  std::vector<unsigned long> residues = std::vector<unsigned long>(256, 0);
  unsigned long size = 0; 
  std::string getID();
  void print();
};


// MAIN
int main(int argc, char* argv[]) {
  if (argc <= 1) {
    // READ FROM: stdin
    Fasta fasta = Fasta("stdin");
    std::cerr << "Cannot Read From stdin... Yet!" << std::endl;
    exit(EXIT_FAILURE);
  } else {
      // READ FROM: file
      std::cout << "Reading File: " << argv[1] << std::endl;
      std::ifstream file_in (argv[1]);
      if (!file_in.is_open()) {
        std::cerr <<  "ERROR: Can't Open File: " << argv[1] << std::endl;
        exit(EXIT_FAILURE);
      }
      
      // Initialize objects
      Fasta fasta = Fasta(std::string(argv[1]));
      Record record;
      std::string line;

      // Validate that first line is a header
      std::getline(file_in, line);
      if (Fasta::isHeader(line)) {
        record = Record(line);
      } else {
          std::cerr << "ERROR: FastA does not start with >header" << std::endl;
          exit(EXIT_FAILURE);
      }
     
      std::cout << "{ \"records\" : [";
      // Continue reading lines
      while (std::getline(file_in, line)) {
        if (Fasta::isHeader(line)) {
          record.print();
          std::cout << ", ";
          record = Record(line);
        } else {
          fasta.size += line.length();
          record.size += line.length();
          for (auto residue: line) {
            // std::cout << c << std::endl;
            ++fasta.residues[residue];
            ++record.residues[residue];
          }
        }
      }
      record.print();
      std::cout << "], \"fasta\": ";
      fasta.print(); 
      std::cout << "}\n";
   }
}

/* Fasta Constructor */
Fasta::Fasta(std::string input): input(input) {}

/* Record Constructor */
Record::Record(std::string header): header(header) {}
Record::Record():header("No Header") {}

/* Print FastA JSON */
void Fasta::print() {
  std::cout << "{\"input\": \"" << input << "\"" 
            << ", \"bp\": " << size; 
  for (int i = 0; i < residues.size(); ++i) {
    if (residues[i] > 0) {
        std::cout << ", \""<< (char) i << "\": " << residues[i];  
    }
  }
  std::cout << "}";
}

/* Print Record JSON */
void Record::print() {
  std::cout << "{\"id\": \"" << getID() 
            << ", \"bp\": " << size; 
  for (int i = 0; i < residues.size(); ++i) {
    if (residues[i] > 0) {
        std::cout << ", \""<< (char) i << "\": " << residues[i];  
    }
  }
  std::cout << "}";
}

/* Helper Function: Headers start with the '>' character */
bool Fasta::isHeader(std::string &line) {
  return line[0] ==  '>';
}

std::string Record::getID() {
  return header.substr(1, header.find(' ') - 1);
}

