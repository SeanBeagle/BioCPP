/*
  FILE: fasta.cpp
AUTHOR: Sean Beagle
   URL: https://www.seanbeagle.com
        http://stronglab.org
DESC:   This program reads a FastA file and returns a JSON of file stats as seen
        below. [NOTE] Actual output is single line JSON.

        {
          "fasta": {"input": "file/path.fasta", "contigs": 2, "bp": 40,
                    "A": 10, "C": 10, "G": 10, "T": 10}, 
          "contigs": [
            {"id": "sequence1", "bp" 20, "A": 5, "C": 5, "G": 5, "T": 5},
            {"id": "sequence2", "bp" 20, "A": 5, "C": 5, "G": 5, "T": 5}
          ]
        }
        
        header: Is a line starting with the '>' character and denotes a new 
                sequence is on the following lines. The header format is assumed
                to be: >id description.

        residues: Are the characters that make up a sequence. 
                  A vector<unsigned long> of length 256 is used to keep tally of
                  all printable characters by ASCII value that may be included
*/

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
  static void readFile(std::string &file);
  static bool isHeader(std::string &line);
};

struct Record {
  Record();
  Record(std::string header);
  std::string header = "No Header";
  std::vector<unsigned long> residues = std::vector<unsigned long>(256, 0);
  unsigned long size = 0; 
  std::string id();
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
    Fasta::readFile(std::string(argv[1]);
  }
}

void Fasta::readFile(std::string &file) {
  // READ FROM: file
  std::ifstream file_in (file);
  if (!file_in.is_open()) {
    std::cerr <<  "ERROR: Can't Open File: " << argv[1] << std::endl;
    exit(EXIT_FAILURE);
  }
      
  // Initialize objects
  Fasta fasta = Fasta(file);
  Record record;
  std::string line;

  // Validate: First line is >header
  // TODO(seanbeagle): Consider single sequence FastA with no header?
  std::getline(file_in, line);
  if (Fasta::isHeader(line)) {
    ++fasta.contigs;
    record = Record(line);
  } else {
    std::cerr << "ERROR: FastA does not start with a >header" << std::endl;
    exit(EXIT_FAILURE);
  }
     
  // Continue reading lines and streaming JSON output
  std::cout << "{\"contigs\": [";
  while (std::getline(file_in, line)) {
    if (Fasta::isHeader(line)) {
      record.print();
      std::cout << ", ";
      ++fasta.contigs;
      record = Record(line);
    } else {
      fasta.size += line.length();
      record.size += line.length();
      for (auto residue: line) {
        ++fasta.residues[residue];
        ++record.residues[residue];
      }
    }
  }
  // Print final record + fasta information, and close JSON
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
            << ", \"bp\": " << size
            << ", \"contigs\": " << contigs; 
  for (int i = 0; i < residues.size(); ++i) {
    if (residues[i] > 0) {
        std::cout << ", \""<< (char) i << "\": " << residues[i];  
    }
  }
  std::cout << "}";
}

/* Print Record JSON */
void Record::print() {
  std::cout << "{\"id\": \"" << id() 
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

std::string Record::id() {
  return header.substr(1, header.find(' ') - 1);
}

