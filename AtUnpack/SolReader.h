#ifndef SOLREADER_H
#define SOLREADER_H

#include <time.h> // time in nano-sec

#include "Hit.h"
#include <stdio.h> /// for FILE
#include <unistd.h>

#include <cstdlib>
#include <string>
#include <vector>

class SolReader {
private:
   FILE *inFile;
   unsigned int inFileSize;
   unsigned int filePos;
   unsigned int totNumBlock;

   unsigned short blockStartIdentifier;
   unsigned int numBlock;
   bool isScanned;

   void init();

   std::vector<unsigned int> blockPos;

public:
   SolReader();
   SolReader(std::string fileName, unsigned short dataType = 0);
   ~SolReader();

   void OpenFile(std::string fileName);
   int ReadNextBlock(int isSkip = 0); // opt = 0, noraml, 1, fast
   int ReadBlock(unsigned int index, bool verbose = false);

   void ScanNumBlock();

   bool IsEndOfFile() const { return (filePos >= inFileSize ? true : false); }
   unsigned int GetBlockID() const { return numBlock - 1; }
   unsigned int GetNumBlock() const { return numBlock; }
   unsigned int GetTotalNumBlock() const { return totNumBlock; }
   unsigned int GetFilePos() const { return filePos; }
   unsigned int GetFileSize() const { return inFileSize; }

   void RewindFile();

   Hit *hit;
};

#endif
