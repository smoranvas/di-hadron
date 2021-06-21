#ifndef HEADERS_HXX
#include "Headers.hxx"
#endif

#ifndef UI_HXX
#define UI_HXX

/*** Global Variables ***/

TString gDataKind = "";

TString gTargetOption;
TString gRunNumber;

// input
TString gInputFile;
TString gOutputFile;

/*** Input/output functions ***/

void printUsage() {
  std::cout << "GetSimpleTuple_" << gDataKind << " program. Usage is:" << std::endl;
  std::cout << std::endl;
  std::cout << "./GetSimpleTuple_" << gDataKind << " -h" << std::endl;
  std::cout << "    prints usage and exit program" << std::endl;
  std::cout << std::endl;
  std::cout << "./GetSimpleTuple_" << gDataKind << " -t[target]" << std::endl;
  if (gDataKind == "data")
    std::cout << "    selects target, which can be: C, Fe, Pb" << std::endl;
  else
    std::cout << "    selects target, which can be: D, C, Fe, Pb" << std::endl;
  std::cout << std::endl;
  std::cout << "./GetSimpleTuple_" << gDataKind << " -r[run number]" << std::endl;
  std::cout << "    selects run number" << std::endl;
  if (gDataKind == "data")
    std::cout << "    numbering scheme for data files = clas_<run number>_*.pass2.root" << std::endl;
  else
    std::cout << "    numbering scheme for sim files = recsis<target>_<run number>.root" << std::endl;
  std::cout << std::endl;
}

void parseCommandLine(int argc, char *argv[]) {
  int c;
  if (argc == 1) {
    std::cerr << "Empty command line. Execute ./bin/GetSimpleTuple_" << gDataKind << " -h to print usage." << std::endl;
    exit(1);
  }
  while ((c = getopt(argc, argv, "ht:r:")) != -1) {
    switch (c) {
      case 'h':
        printUsage();
        exit(0);
        break;
      case 't':
        gTargetOption = optarg;
        break;
      case 'r':
        gRunNumber = optarg;
        break;
      default:
        std::cerr << "Unrecognized argument. Execute ./bin/GetSimpleTuple_" << gDataKind << " -h to print usage." << std::endl;
        exit(0);
        break;
    }
  }
}

void printOptions() {
  std::cout << "Executing GetSimpleTuple_" << gDataKind << " program. The chosen parameters are: " << std::endl;
  std::cout << "  TargetOption   = " << gTargetOption << std::endl;
  std::cout << "  RunNumber      = " << gRunNumber << std::endl;
  std::cout << std::endl;
}

#endif
