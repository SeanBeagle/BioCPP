/* matrix-metrics.cpp
AUTHOR: Sean Beagle
URL:    www.seanbeagle.com
        www.stronglab.org
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <string>
#include <vector>
#include <array>
#include <cctype>

// FORWARD DECLARATIONS
class RootMatrix;
class SeqRecord;
class Matrix;
class DataFrame;
class Row;


// HEADERS
class RootMatrix {
  public:
    size_t numRecords();
    size_t numPositions();
    std::string fasta();
    SeqRecord operator[](size_t record);
    char operator()(size_t record, size_t position);
    char at(size_t record, size_t position); // Replaced by operator[]
    static RootMatrix fromFasta(std::string fasta);
  private:
    size_t num_records_;
    size_t num_positions_;
    std::string fasta_;
    std::string matrix_;
    std::vector<SeqRecord> records_;
    size_t n = -1;
    RootMatrix(std::string fasta, size_t num_records, size_t num_positions);
    void newRecord(std::string &header);
    void addSequence(std::string &seq);
};

class SeqRecord {
  public:
    SeqRecord(std::string header, size_t index, RootMatrix* matrix);
    char operator[](size_t position);
    size_t operator[](char residue);
    std::string id();
    size_t index();
    std::string description();
    void countResidues(std::string &seq);
    std::string const json();
  private:
    std::array<size_t, 127> residues_ = {};
    std::string header_;
    size_t index_;
    RootMatrix* matrix_;
};

class Matrix {
  public:
    Matrix(size_t num_positions, std::vector<SeqRecord> records);
    Matrix(RootMatrix &root);
    std::string fasta;
    size_t const numRecords();
    size_t const numPositions();
    size_t const coreSize();
    size_t const numCoreSNP();
    size_t const numSNP();
    void operator+=(SeqRecord);
    SeqRecord operator[](size_t record);
    void toJSON();
    std::vector<SeqRecord>::const_iterator const begin();
    std::vector<SeqRecord>::const_iterator const end();
  private:
    size_t num_records_ = records_.size();
    size_t num_positions_;
    std::vector<SeqRecord> records_;
    DataFrame dataframe_;
};

class DataFrame {
  public:
    DataFrame(size_t size);
    bool isSNP(size_t position);
    bool isCoreSNP(size_t position);
  private:
    Matrix* matrix_;
    std::vector<Row> rows_;
    std::vector<bool> is_snp_;
    std::vector<bool> is_core_snp_;
    void calcPosStats();
    void calcOutlierStats();
};

class Row {
  public:
    Row(SeqRecord *record);
  private:
    SeqRecord* record_;
    size_t unique;
};


// MAIN ============================================
int main(int argc, char* argv[]) {
  if (argc > 1) {
    std::string fasta = std::string(argv[1]);
    RootMatrix root = RootMatrix::fromFasta(fasta);
    Matrix matrix = Matrix(root);
    
    for (auto record: matrix) {
      std::cout << record.json();
    }

    return EXIT_SUCCESS;
  }
}


/*******************************************************************************
 class RootMatrix
*******************************************************************************/

/* Private RootMatrix Constructor */
RootMatrix::RootMatrix(
  std::string fasta, size_t num_records, size_t num_positions): 
  num_positions_(num_positions), fasta_(fasta), num_records_(num_records) {
  // READ FASTA
  std::ifstream file_in (fasta);
  std::string line;
  while (std::getline(file_in, line)) {
    if (line[0] == '>') {
      newRecord(line);
    } else {
      addSequence(line);
    }
  }
}
                         
/* Static class method that Returns RootMatrix object from FastA file */                         
RootMatrix RootMatrix::fromFasta(std::string fasta) {
  std::ifstream file_in (fasta);
  // VALIDATE FILE
  if (!file_in.is_open()) {
    std::cerr << "ERROR: Can't open " << fasta << std::endl;
    exit(EXIT_FAILURE);
  }
  // GET MATRIX PARAMETERS
  size_t num_records=0, prev_size=0, this_size=0; 
  std::string line;
  while (std::getline(file_in, line)) {
    if (line[0] == '>') {
      if (prev_size > 0 && prev_size != this_size) {
        std::cerr << "ERROR: Not all sequences are the same length\n";
        exit(EXIT_FAILURE);
      }
      ++num_records;
      prev_size = this_size;
      this_size = 0;
    } else {
      this_size += line.length();
    }
  }
  return RootMatrix(fasta, num_records, this_size);
}

void RootMatrix::newRecord(std::string &header) {
  records_.emplace_back(SeqRecord(header, ++n, this));
}

void RootMatrix::addSequence(std::string &seq) {
  matrix_ += seq;
  records_.back().countResidues(seq);
}

SeqRecord RootMatrix::operator[](size_t record) {
  return records_[record];
}

char RootMatrix::operator()(size_t record, size_t position) {
  return matrix_[record*numPositions() + position];
}

size_t RootMatrix::numRecords() {
  return records_.size();
}

size_t RootMatrix::numPositions() {
  return num_positions_;
}

std::string RootMatrix::fasta() {
  return fasta_;
}


/*******************************************************************************
 class SeqRecord
 ******************************************************************************/

SeqRecord::SeqRecord(
  std::string header, size_t index, RootMatrix* matrix) 
  : header_(header), index_(index), matrix_(matrix) {}

void SeqRecord::countResidues(std::string &seq) {
  for (auto residue: seq)
    ++residues_[std::toupper(residue)];
}

std::string const SeqRecord::json() {
  std::stringstream ss;
  ss << "{\"id\": " << "\"" << id() << "\"";
  for (int i = 0; i < residues_.size(); ++i) {
    if (residues_[i] > 0) {
      ss << ", \"" << (char)i << "\": " << residues_[i];
    }
  }
  ss << "}" << std::endl;
  return ss.str();
}

char SeqRecord::operator[](size_t position) {
  return (*matrix_)(index_, position);
}

size_t SeqRecord::operator[](char residue) {
  return residues_[std::toupper(residue)];
}

size_t SeqRecord::index() {
  return index_;
}

std::string SeqRecord::id() {
  return header_.substr(1, header_.find(" ")-1);
}

std::string SeqRecord::description() {
  return header_.substr(header_.find(" ")+1);
}
  

/*******************************************************************************
 class Matrix
*******************************************************************************/

/* Construct Matrix from SeqRecord vector */
Matrix::Matrix(size_t num_positions, std::vector<SeqRecord> records) {}

/* Construct Matrix from all RootMatrix records */
Matrix::Matrix(RootMatrix &root): 
  num_positions_(root.numPositions()), num_records_(root.numRecords()) {
  for (int i = 0; i  < root.numRecords(); ++i) {
    records_.push_back(root[i]);
  }
}

void Matrix::toJSON() {
  std::cout << "{" << std::endl;
  for (auto &record: records_) {
    std::cout << "  " << record.json();
  }
  std::cout << "}";
}

std::vector<SeqRecord>::const_iterator const Matrix::begin() {
  return records_.begin();
}
std::vector<SeqRecord>::const_iterator const Matrix::end() {
  return records_.end();
}

SeqRecord Matrix::operator[](size_t record) {
  return records_[record];
}

size_t const Matrix::numRecords() {
  return records_.size();
}

size_t const Matrix::numPositions() {
  return num_positions_;
}