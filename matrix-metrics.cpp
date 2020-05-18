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
// class DataFrame;  // TODO(seanbeagle): Use this for matrix stats...
// class Seq; // TODO(seanbeagle): Use this instead of string for sequences...


// HEADERS
class RootMatrix {
  public:
    const size_t num_records;
    const size_t num_positions;
    const std::string fasta;
    SeqRecord operator[](size_t index);
    char at(size_t index);
    static RootMatrix fromFasta(std::string fasta);
  private:
    std::shared_ptr<RootMatrix> ptr;
    size_t n = -1;
    RootMatrix(std::string fasta, size_t num_records, size_t num_positions);
    std::vector<std::shared_ptr<SeqRecord>> records;
    void addHeader(std::string &header);
    void addSequence(std::string &seq);
    std::string matrix;
};

class SeqRecord {
  public:
    SeqRecord(std::string header, unsigned index, 
              std::shared_ptr<RootMatrix> matrix);
    const std::string header;
    const unsigned index;
    char operator[](unsigned index);
  private:
    std::shared_ptr<RootMatrix> matrix;
};

class Matrix {
  public:
    Matrix(size_t num_positions, 
           std::vector<std::shared_ptr<SeqRecord>> records);
    std::string const fasta;
    unsigned numRecords();
    unsigned numPositions();
    void operator+=(SeqRecord);
    void operator+=(std::string sequence);
    std::shared_ptr<SeqRecord> operator[](size_t index);
    char at(unsigned index);
  private:
    size_t num_records = 0;
    size_t num_positions = 0;
    std::shared_ptr<std::string> matrix;
    std::vector<std::shared_ptr<SeqRecord>> records;
};


// MAIN
int main(int argc, char* argv[]) {
  if (argc > 1) {
    std::string fasta = std::string(argv[1]);
    RootMatrix root = RootMatrix::fromFasta(fasta);
    std::cout << root.fasta << " = " <<  root.num_records << " x " 
              << root.num_positions << std::endl;

    // for (int i = 0; i < root.num_records; ++i) {
    //   std::cout << root[i].header << std::endl;
    // }

    std::cout << "root[10].header = " << root[10].header << std::endl;
    std::cout << "root.at(10) = " << root.at(10);

    return EXIT_SUCCESS;
    //std::cout << "matrix[0][10] = " << root[0][10] << std::endl;
  }
}


/*******************************************************************************
 class RootMatrix
 ******************************************************************************/

/* Private RootMatrix Constructor */
RootMatrix::RootMatrix(
  std::string fasta, size_t num_records, size_t num_positions): 
  num_positions(num_positions), fasta(fasta), num_records(num_records) {
  std::cout << "new RootMatrix(\"" << fasta << "\")\n";  // TODO: REMOVE LINE
  ptr = std::shared_ptr<RootMatrix> (this);
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
  // BUILD ROOT MATRIX
  return RootMatrix(fasta, num_records, this_size);
}

void RootMatrix::addHeader(std::string &header) {
  std::shared_ptr<SeqRecord> sp (new SeqRecord(header, ++n, ptr));
  records.push_back(sp);
}

void RootMatrix::addSequence(std::string &seq) {
  matrix += seq;
}

SeqRecord RootMatrix::operator[](size_t index) {
  return *records[index];
}

char RootMatrix::at(size_t index) {
  return 'x';
  // return matrix[index];
}


/*******************************************************************************
 class SeqRecord
 ******************************************************************************/

/* SeqRecord Constructor */
SeqRecord::SeqRecord(
  std::string header, unsigned index, std::shared_ptr<RootMatrix> matrix) 
  : header(header), index(index) {
    std::cout << "new SeqRecord(\"" << header << "\")\n"; // TODO: REMOVE LINE
  }


char SeqRecord::operator[](unsigned index) {
  return 'x';
  // RootMatrix root = *matrix;
  // std::cout << root.fasta;
  // std::cout << matrix->fasta;
  // return 'x';
  // return root.at(root.num_positions * this->index + index);
}

/*******************************************************************************
 class Matrix
 ******************************************************************************/

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




/* Matrix getters */
// char Matrix::at(unsigned index) { return matrix->at(index); }
// SeqRecord Matrix::operator[](unsigned index) { return records[index]; }
// unsigned Matrix::numRecords() { return num_records; }
// unsigned Matrix::numPositions() { return num_positions; }

