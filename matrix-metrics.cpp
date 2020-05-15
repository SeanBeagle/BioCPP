#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

// FORWARD DECLARATIONS
class RootMatrix;
class Matrix;
class SeqRecord;
class DataFrame;

// HEADERS
class SeqRecord {
  public:
    SeqRecord(std::string header, unsigned index, Matrix* matrix);
    const std::string header;
    const unsigned index;
    char operator[](unsigned index);
  private:
    Matrix* matrix;
};

class Matrix {
  public:
    Matrix(std::string fasta);
    Matrix(std::shared_ptr<std::string> matrix, std::vector<SeqRecord> records);
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
    std::shared_ptr<std::string> matrix_;              //  LEARN TO USE THESE
    std::vector<std::shared_ptr<SeqRecord>> records_;  //    SMART POINTERS
    std::string* matrix;
    std::vector<SeqRecord> records;
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

