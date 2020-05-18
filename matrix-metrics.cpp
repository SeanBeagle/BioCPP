#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <vector>
#include <memory>

// FORWARD DECLARATIONS
class RootMatrix;
class Matrix;
class SeqRecord;
class DataFrame;

// HEADERS
class RootMatrix {
  public:
    const size_t num_records;
    const size_t num_positions;
    const std::string fasta;
    SeqRecord operator[](unsigned index);
    static RootMatrix fromFasta(std::string fasta);
  private:
    std::vector<std::shared_ptr<SeqRecord>> records;
    RootMatrix(std::string fasta, unsigned num_records, unsigned num_positions);
}

class SeqRecord {
  public:
    SeqRecord(std::string header, unsigned index, 
              std::shared_pointer<RootMatrix> matrix);
    const std::string header;
    const unsigned index;
    char operator[](unsigned index);
  private:
    std::shared_pointer<RootMatrix> matrix;
};

class Matrix {
  public:
    Matrix(size_t num_positions, 
           std::vector<std::shared_pointer<SeqRecord>> records);
    std::string const fasta;
    unsigned numRecords();
    unsigned numPositions();
    void operator+=(SeqRecord);
    void operator+=(std::string sequence);
    SeqRecord operator[](unsigned index);
    char at(unsigned index);
  private:
    unsigned num_records = 0;
    unsigned num_positions = 0;
    std::shared_ptr<std::string> matrix;              //  LEARN TO USE THESE
    std::vector<std::shared_ptr<SeqRecord>> records;  //    SMART POINTERS
};


// MAIN
int main(int argc, char* argv[]) {
  if (argc > 1) {
    std::string fasta = argv[1];
    Matrix m = Matrix(fasta);
    std::cout << m.fasta << " = " <<  m.numRecords() << " x " 
              << m.numPositions() << std::endl;

    for (int i = 0; i < m.numRecords(); ++i) {
      std::cout << m[i].header << std::endl;
    }
    
    std::cout << "matrix[0][10] = " << m[0][10] << std::endl;
  }
}

/* RootMatrix Constructor */
RootMatrix::RootMatrix(std::string fasta, unsigned num_records, 
                       unsigned num_positions) 
                       : num_positions(num_positions), num_records(num_records), 
                         fasta(fasta) {}
                         
                         
static RootMatrix RootMatrix::fromFasta(std::string fasta) {
  std::ifstream file_in (fasta);
  if (!file_in.is_open()) {
    std::cerr << "ERROR: Can't open " << fasta << std::endl;
    exit(EXIT_FAILURE);
  }
 
  size_t records=0, prev_size=0, this_size=0; 
  std::string line;
  while (std::getline(file_in, line)) {
    if (line[0] == '>') {
      if (prev_size > 0 && prev_size != this_size) {
        std::cerr << "ERROR: Not all sequences are the same length\n";
        exit(EXIT_FAILURE);
      }
      ++records;
      prev_size = this_size;
      this_size = 0;
    } else {
      this_size += line.length();
    }
  }
}

/* SeqRecord Constructor */
SeqRecord::SeqRecord(std::string header, unsigned index, Matrix* matrix) 
  : header(header), index(index) {}

char SeqRecord::operator[](unsigned index) {
  return matrix->at(matrix->numPositions() * this->index + index);
}

/* Matrix Constructor from FastA file format */  
Matrix::Matrix(std::string fasta): fasta(fasta) {
  this->matrix = (new std::string());
  std::ifstream file(fasta);
  std::string line;
  while (std::getline(file, line)) {
    if (line[0] == '>') {
      records.emplace_back(SeqRecord(line, num_records++, this));
    } else {
      num_positions += line.length();
      *matrix += line;
    }
  }
  std::cout << "Finished loading " << matrix->size() << " nucleotides\n";
} 

/* Matrix getters */
char Matrix::at(unsigned index) { return matrix->at(index); }
SeqRecord Matrix::operator[](unsigned index) { return records[index]; }
unsigned Matrix::numRecords() { return num_records; }
unsigned Matrix::numPositions() { return num_positions; }

