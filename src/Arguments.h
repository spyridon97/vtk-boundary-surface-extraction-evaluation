//
// Created by spiros.tsalikis on 9/24/24.
//

#ifndef CLI_H
#define CLI_H

#include <string>

struct Arguments
{
  std::string InputFileName;
  unsigned int NumberOfThreads = 1;
  std::string DeviceName = "TBB";
  unsigned int NumberOfTrials = 1;
  bool Randomize = false;

  bool HashDistribution = false;
  bool SClassifier = false;
  bool SHash = false;
  bool PClassifier = false;
  bool PHash = false;
  bool PHashFight = false;
  bool PHashSort = false;
  bool PHashCount = false;

  int HashFunction = 0;

  /**
   * @brief Parse command line arguments.
   *
   * @return The command line arguments
   */
  void ParseArguments(int argc, char** argv);
};

#endif // CLI_H
