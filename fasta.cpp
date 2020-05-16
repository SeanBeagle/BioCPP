#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstddef>
#include <algorithm>

// FORWARD DECLARATION
struct Fasta;
struct Record;


// HEADERS
struct Fasta {
  Fasta(std::string input);
  std::string input;
  std::vector<Record> records;
  std::map<char, unsigned> residues;
  std::size_t size = 0;
  std::string toJSON();
  void addResidue(char a);
  void newRecord(std::string &header);
  void addSequence(std::string &sequence);
  static bool isHeader(std::string &line);
};

struct Record {
  Record(std::string header);
  std::string header;
  std::map<char, unsigned> residues;
  std::size_t size = 0; 
  std::string toJSON();
};


// MAIN
int main(int argc, char* argv[]) {
  if (argc < 1) {
    // READ FROM STDIN
    Fasta fa = Fasta("stdin");
    std::cerr << "Cannot Read From stdin Yet!\n";
    exit(EXIT_FAILURE);
  } else {
    // READ FROM FILE
      std::cout << "Reading File: " << argv[1] << std::endl;
      Fasta fa = Fasta(argv[1]);
      std::ifstream fasta (argv[1]);
      if (!fasta.is_open()) {
        std::cerr <<  "ERROR: Can't Open File: " << argv[1] << std::endl;
      }

      std::string line;
      while (std::getline(fasta, line)) {
        if (Fasta::isHeader(line)) {
          fa.newRecord(line);
        } else {
            fa.addSequence(line); 
        }
      }
      std::cout << fa.toJSON() << std::endl;
  }
}

/* Fasta Constructor */
Fasta::Fasta(std::string input) : input(input) {}

/* Record Constructor */
Record::Record(std::string header): header(header) {}

/* Add Record to Fasta */
void Fasta::newRecord(std::string &header) {
  records.emplace_back(header.substr(1));
}

/* Call addResidue for each char in sequence */
void Fasta::addSequence(std::string &sequence) {
  for (char c: sequence) {
    addResidue(c);  
  }
}

/* Increment residue count and size for Fasta and last Record */
void Fasta::addResidue(char a) {
  if (records.empty()) {
    records.emplace_back(Record("No Header"));
  }
  ++records.back().residues[a];
  ++records.back().size;
  ++residues[a];
  ++size;
}

std::string Fasta::toJSON() {
  std::string json;
  json += std::to_string(records.size()) + " records\n";  
  json += std::to_string(size) + " basepairs\n";
  // json += "fasta = {\"input\"=\"" + input + "\", \"basepairs\"=\"" + size + "\"}";
  return json;
}

/* Helper Function: FastA headers start with the '>' character */
bool Fasta::isHeader(std::string &line) {
  return line[0] ==  '>';
}

