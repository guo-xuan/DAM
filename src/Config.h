/*
 * Config.h
 *
 *  Created on: Apr 23, 2016
 *      Author: Naux
 */

#ifndef CONFIG_H_
#define CONFIG_H_

// C headers:
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

// C++ headers:
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <utility>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <streambuf>
#include <map>
#include <exception>

using namespace std;

typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned int UINT32;
typedef int INT32;
typedef unsigned long long UINT64;
typedef long long INT64;

class Config {
public:
	static vector<string> vsInputVariantFiles;
	static string sOutputFolder;
	static UINT32 iNumVariants;
	static UINT32 iNumGroups;
	static UINT32 iNumChains;
	static double iNumIndependentVariants;
	static double iNumDependentVariants;
	static UINT32 iMaxVariants;

	Config();
	virtual ~Config();

	static void printHelp();
	static bool setConfig(int _argc, char ** _argv);
};

#endif /* CONFIG_H_ */
