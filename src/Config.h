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

// C++ headers:
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#elif _LINUX
#include <stdarg.h>
#include <sys/stat.h>
#endif

#ifdef _WIN32
#define ACCESS _access
#define MKDIR(a) _mkdir((a))
#elif _LINUX
#define ACCESS access
#define MKDIR(a) mkdir((a),0755)
#endif

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

	Config();
	virtual ~Config();

	static void printHelp();
	static bool setConfig(int _argc, char ** _argv);
};

#endif /* CONFIG_H_ */