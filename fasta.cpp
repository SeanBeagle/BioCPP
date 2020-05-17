#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
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
  unsigned long size = 0;
  std::string toJSON();
  void print();
  void addResidue(char c);
  void newRecord(std::string &header);
  void addSequence(std::string &sequence);
  static bool isHeader(std::string &line);
};

struct Record {
  Record(std::string header);
  std::string header;
  std::map<char, unsigned> residues;
  unsigned long size = 0; 
  std::string toJSON();
};


// MAIN
int main(int argc, char* argv[]) {
  if (argc < 1) {
    // READ FROM: stdin
    Fasta fa = Fasta("stdin");
    std::cerr << "Cannot Read From stdin Yet!\n";
    exit(EXIT_FAILURE);
  } else {
      // READ FROM: file
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
      fa.print();
  }
}

/* Fasta Constructor */
Fasta::Fasta(std::string input): input(input) {}

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
void Fasta::addResidue(char c) {
  if (records.empty()) {
    records.emplace_back(Record("No Header"));
  }
  ++records.back().residues[c];
  ++residues[c];
  ++records.back().size;
  ++size;
}

std::string Fasta::toJSON() {
  std::string json;
  json += "{";
  json += "\"input\": \"" + input + "\"";
  json += ", \"contigs\": " + std::to_string(records.size());  
  json += ", \"bp\": " + std::to_string(size);
  json += "}";
  // json += "fasta = {\"input\"=\"" + input + "\", \"basepairs\"=\"" + size + "\"}";
  return json;
}
void Fasta::print() {
  std::cout << "{" << std::endl;
  std::cout << "  \"fasta\": " << toJSON() << "," << std::endl;
  std::cout << "  \"records\": [" << std::endl; 
  for (int i = 0; i < records.size(); ++i) {
    std::cout << "    "<< records[i].toJSON();
    if (i < records.size()-1)
       std::cout << ','; 
    std::cout << std::endl;

  }
  std::cout << "  ]\n}";
}

std::string Record::toJSON() {
  std::string json;
  json += "{";
  json += "\"id\": \"" + header + "\"";
  json += ", \"bp\": " + std::to_string(size);
  json += "}";
  return json;
}

/* Helper Function: FastA headers start with the '>' character */
bool Fasta::isHeader(std::string &line) {
  return line[0] ==  '>';
}

