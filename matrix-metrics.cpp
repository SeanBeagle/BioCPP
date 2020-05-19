#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <vector>

// FORWARD DECLARATIONS
class RootMatrix;
class Matrix;
class SeqRecord;
// class DataFrame;  // TODO(seanbeagle): Use this for matrix stats...
// class Seq; // TODO(seanbeagle): Use this instead of string for sequences...


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
    void addHeader(std::string &header);
    void addSequence(std::string &seq);
};

class SeqRecord {
  public:
    SeqRecord(std::string header, size_t index, RootMatrix* matrix);
    char operator[](size_t position);
    std::string id();
    size_t index();
    std::string description();
  private:
    std::string header_;
    size_t index_;
    RootMatrix* matrix_;
};

class Matrix {
  public:
    Matrix(size_t num_positions, std::vector<SeqRecord*> records);
    Matrix(RootMatrix &root);
    std::string const fasta;
    size_t numRecords();
    size_t numPositions();
    void operator+=(SeqRecord);
    SeqRecord operator[](size_t record);
  private:
    size_t num_records_ = records_.size();
    size_t num_positions_;
    std::vector<SeqRecord*> records_;
};


// MAIN
int main(int argc, char* argv[]) {
  if (argc > 1) {
    std::string fasta = std::string(argv[1]);
    RootMatrix root = RootMatrix::fromFasta(fasta);
    std::cout << root.fasta() << " = " <<  root.numRecords() << " x " 
              << root.numPositions() << std::endl;

    Matrix matrix = Matrix(root);
    std::cout << matrix.numRecords() << "x" << matrix.numPositions() << std::endl;

    // print all SeqRecord position 0's
    for (int i = 0; i < root.numRecords(); ++i) {
      std::cout << root[i][0] << ",";
    }
    std::cout << std::endl;
    
    // print SeqRecord ID's
    for (int i = 0; i < root.numRecords(); ++i) {
      std::cout << root[i].id() << ", ";
    }
    std::cout << std::endl;

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
      addHeader(line);
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

void RootMatrix::addHeader(std::string &header) {
  records_.emplace_back(SeqRecord(header, ++n, this));
}

void RootMatrix::addSequence(std::string &seq) {
  matrix_ += seq;
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

/* SeqRecord Constructor */
SeqRecord::SeqRecord(
  std::string header, size_t index, RootMatrix* matrix) 
  : header_(header), index_(index), matrix_(matrix) {}

char SeqRecord::operator[](size_t position) {
  return (*matrix_)(index_, position);
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

Matrix::Matrix(size_t num_positions, std::vector<SeqRecord*> records) {}

Matrix::Matrix(RootMatrix &root): 
  num_positions_(root.numPositions()), num_records_(root.numRecords()) {
  std::cout << "new Matrix()" << std::endl;
  for (int i = 0; i  < root.numRecords(); ++i) {
    records_.push_back(&root[i]);
  }
}

size_t Matrix::numRecords() {
  return records_.size();
}

size_t Matrix::numPositions() {
  return num_positions_;
}


/* Matrix Constructor from FastA file format */  
// Matrix::Matrix(std::string fasta): fasta(fasta) {
//   this->matrix = (new std::string());
//   std::ifstream file(fasta);
//   std::string line;
//   while (std::getline(file, line)) {
//     if (line[0] == '>') {
//       records.emplace_back(SeqRecord(line, num_records++, this));
//     } else {
//       num_positions += line.length();
//       *matrix += line;
//     }
//   }
//   std::cout << "Finished loading " << matrix->size() << " nucleotides\n";
// } 






