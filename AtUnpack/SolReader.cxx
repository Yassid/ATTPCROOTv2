#include "SolReader.h"

void SolReader::init()
{
   inFileSize = 0;
   numBlock = 0;
   filePos = 0;
   totNumBlock = 0;
   hit = new Hit();

   isScanned = false;

   blockPos.clear();
}

SolReader::SolReader()
{
   init();
}

SolReader::SolReader(std::string fileName, unsigned short dataType)
{
   init();
   OpenFile(fileName);
   hit->SetDataType(dataType, DPPType::PHA);
}

SolReader::~SolReader()
{
   if (!inFile)
      fclose(inFile);
   delete hit;
}

void SolReader::OpenFile(std::string fileName)
{
   inFile = fopen(fileName.c_str(), "r");
   if (inFile == NULL) {
      printf("Cannot open file : %s \n", fileName.c_str());
   } else {
      fseek(inFile, 0L, SEEK_END);
      inFileSize = ftell(inFile);
      rewind(inFile);
   }
}

int SolReader::ReadBlock(unsigned int index, bool verbose)
{
   if (isScanned == false)
      return -1;
   if (index >= totNumBlock)
      return -1;
   fseek(inFile, 0L, SEEK_SET);

   if (verbose)
      printf("Block index: %u, File Pos: %u byte\n", index, blockPos[index]);

   fseek(inFile, blockPos[index], SEEK_CUR);

   filePos = blockPos[index];

   numBlock = index;

   return ReadNextBlock();
}

int SolReader::ReadNextBlock(int isSkip)
{
   if (inFile == NULL)
      return -1;
   if (feof(inFile))
      return -1;
   if (filePos >= inFileSize)
      return -1;

   fread(&blockStartIdentifier, 2, 1, inFile);

   if ((blockStartIdentifier & 0xAA00) != 0xAA00) {
      printf("header fail.\n");
      return -2;
   }

   if ((blockStartIdentifier & 0xF) == DataFormat::Raw) {
      hit->SetDataType(DataFormat::Raw, ((blockStartIdentifier >> 1) & 0xF) == 0 ? DPPType::PHA : DPPType::PSD);
   }
   hit->dataType = blockStartIdentifier & 0xF;
   hit->DPPType = ((blockStartIdentifier >> 4) & 0xF) == 0 ? DPPType::PHA : DPPType::PSD;

   if (hit->dataType == DataFormat::ALL) {
      if (isSkip == 0) {
         fread(&hit->channel, 1, 1, inFile);
         fread(&hit->charge, 2, 1, inFile);
         if (hit->DPPType == DPPType::PSD)
            fread(&hit->energy_short, 2, 1, inFile);
         fread(&hit->timestamp, 6, 1, inFile);
         fread(&hit->fine_timestamp, 2, 1, inFile);
         fread(&hit->flags_high_priority, 1, 1, inFile);
         fread(&hit->flags_low_priority, 2, 1, inFile);
         fread(&hit->downSampling, 1, 1, inFile);
         fread(&hit->board_fail, 1, 1, inFile);
         fread(&hit->flush, 1, 1, inFile);
         fread(&hit->trigger_threashold, 2, 1, inFile);
         fread(&hit->event_size, 8, 1, inFile);
         fread(&hit->aggCounter, 4, 1, inFile);
      } else {
         fseek(inFile, hit->DPPType == DPPType::PHA ? 31 : 33, SEEK_CUR);
      }
      fread(&hit->traceLenght, 8, 1, inFile);
      if (isSkip == 0) {
         fread(hit->analog_probes_type, 2, 1, inFile);
         fread(hit->digital_probes_type, 4, 1, inFile);
         fread(hit->analog_probes[0], hit->traceLenght * 4, 1, inFile);
         fread(hit->analog_probes[1], hit->traceLenght * 4, 1, inFile);
         fread(hit->digital_probes[0], hit->traceLenght, 1, inFile);
         fread(hit->digital_probes[1], hit->traceLenght, 1, inFile);
         fread(hit->digital_probes[2], hit->traceLenght, 1, inFile);
         fread(hit->digital_probes[3], hit->traceLenght, 1, inFile);
      } else {
         fseek(inFile, 6 + hit->traceLenght * (12), SEEK_CUR);
      }

   } else if (hit->dataType == DataFormat::OneTrace) {
      if (isSkip == 0) {
         fread(&hit->channel, 1, 1, inFile);
         fread(&hit->charge, 2, 1, inFile);
         if (hit->DPPType == DPPType::PSD)
            fread(&hit->energy_short, 2, 1, inFile);
         fread(&hit->timestamp, 6, 1, inFile);
         fread(&hit->fine_timestamp, 2, 1, inFile);
         fread(&hit->flags_high_priority, 1, 1, inFile);
         fread(&hit->flags_low_priority, 2, 1, inFile);
      } else {
         fseek(inFile, hit->DPPType == DPPType::PHA ? 14 : 16, SEEK_CUR);
      }
      fread(&hit->traceLenght, 8, 1, inFile);
      if (isSkip == 0) {
         fread(&hit->analog_probes_type[0], 1, 1, inFile);
         fread(hit->analog_probes[0], hit->traceLenght * 4, 1, inFile);
      } else {
         fseek(inFile, 1 + hit->traceLenght * 4, SEEK_CUR);
      }

   } else if (hit->dataType == DataFormat::NoTrace) {
      if (isSkip == 0) {
         fread(&hit->channel, 1, 1, inFile);
         fread(&hit->charge, 2, 1, inFile);
         if (hit->DPPType == DPPType::PSD)
            fread(&hit->energy_short, 2, 1, inFile);
         fread(&hit->timestamp, 6, 1, inFile);
         fread(&hit->fine_timestamp, 2, 1, inFile);
         fread(&hit->flags_high_priority, 1, 1, inFile);
         fread(&hit->flags_low_priority, 2, 1, inFile);
      } else {
         fseek(inFile, hit->DPPType == DPPType::PHA ? 14 : 16, SEEK_CUR);
      }

   } else if (hit->dataType == DataFormat::MiniWithFineTime) {
      if (isSkip == 0) {
         fread(&hit->channel, 1, 1, inFile);
         fread(&hit->charge, 2, 1, inFile);
         if (hit->DPPType == DPPType::PSD)
            fread(&hit->energy_short, 2, 1, inFile);
         fread(&hit->timestamp, 6, 1, inFile);
         fread(&hit->fine_timestamp, 2, 1, inFile);
      } else {
         fseek(inFile, hit->DPPType == DPPType::PHA ? 11 : 13, SEEK_CUR);
      }

   } else if (hit->dataType == DataFormat::Minimum) {
      if (isSkip == 0) {
         fread(&hit->channel, 1, 1, inFile);
         fread(&hit->charge, 2, 1, inFile);
         if (hit->DPPType == DPPType::PSD)
            fread(&hit->energy_short, 2, 1, inFile);
         fread(&hit->timestamp, 6, 1, inFile);
      } else {
         fseek(inFile, hit->DPPType == DPPType::PHA ? 9 : 11, SEEK_CUR);
      }

   } else if (hit->dataType == DataFormat::Raw) {
      fread(&hit->dataSize, 8, 1, inFile);
      if (isSkip == 0) {
         fread(hit->data, hit->dataSize, 1, inFile);
      } else {
         fseek(inFile, hit->dataSize, SEEK_CUR);
      }
   }

   numBlock++;
   filePos = ftell(inFile);
   return 0;
}

void SolReader::RewindFile()
{
   rewind(inFile);
   filePos = 0;
   numBlock = 0;
}

void SolReader::ScanNumBlock()
{
   if (inFile == NULL)
      return;
   if (feof(inFile))
      return;

   numBlock = 0;
   blockPos.clear();

   blockPos.push_back(0);

   while (ReadNextBlock(1) == 0) {
      blockPos.push_back(filePos);
      printf("%u, %.2f%% %u/%u\n\033[A\r", numBlock, filePos * 100. / inFileSize, filePos, inFileSize);
   }

   totNumBlock = numBlock;
   numBlock = 0;
   isScanned = true;
   printf("\nScan complete: number of data Block : %u\n", totNumBlock);
   rewind(inFile);
   filePos = 0;

   // for( int i = 0; i < totNumBlock; i++){
   //  printf("%7d | %u \n", i, blockPos[i]);
   //}
}
