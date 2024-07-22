#pragma once
//#include <cstdlib>
//#include <iostream> // C++ I/O library header
#include "Random_BG.h"

class NodesInGraph{
 public:
  int numOfConnectionParityBit;
  int numOfConnectionSymbolBit;
  int maxDegParity;
  int *connectionParityBit;
  int *connectionSymbolBit;

  NodesInGraph(void);
  ~NodesInGraph(void);
  void setNumOfConnectionSymbolBit(int deg);
  void initConnectionParityBit(void);
  void initConnectionParityBit(int deg);
};

class BigGirth {
 public:
  int M, N,Z;
  int K;
  int EXPAND_DEPTH;
  int *(*H);

  int *localGirth;
  
  NodesInGraph *nodesInGraph;
  Random *myrandom;

  BigGirth(const int m, const int n,
	  const int *symbolDegSequence,  
	  const int sglConcent,
	  const int tgtGirth=100,
	  const int Z=1, const bool verbose_ = true, const int seed=-1, const int*checkDegSequence=0);
  BigGirth(void);



  void writeToFile_Hcompressed(const char *filename, const int Z=1);
  void writeToFile_Hmatrix(const char *filename );
  void writeToFile(const char *filename );
  void BuildMatrices(
	  int &maxrow, //! Maximum number of rows
	  int* &generator_compressed,   //!< pointer to compressed generator matrix
	  int &maxcol,					//
	  int* &parityCheck_compressed  //!< pointer to compressed parity check matrix matrix
	  );
  void loadH(void);
  int* Build_Hcompressed(int& ncheck);

  ~BigGirth(void);

 private:
  int selectParityConnect(int kthSymbol, int mthConnection, int Z, int & cycle);
  void updateConnection(int kthSymbol);
  bool verbose;

};
