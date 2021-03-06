/*
 * Config.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: Naux
 */

#include "Config.h"

Config::Config() {

}

Config::~Config() {

}

vector<string> Config::vsInputVariantFiles;
string Config::sOutputFolder = "";
UINT32 Config::iNumVariants = 0;
UINT32 Config::iNumGroups = 0;
UINT32 Config::iNumChains = 1;
double Config::iNumIndependentVariants = 5;
double Config::iNumDependentVariants = 5;
UINT32 Config::iMaxVariants = 20;
UINT32 Config::iMaxVariantsInInteraction = 5;

void Config::printHelp() {

}

bool Config::setConfig(int _argc, char ** _argv) {

	vector<string> vsArgumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for (int i = 0; i < _argc; i++) {
		cout << _argv[i] << ' ';
	}
	cout << endl;

	while (_argc--) {
		vsArgumentsList.push_back(*_argv++);
	}

	if (vsArgumentsList.size() == 1) {
		Config::printHelp();
		return false;
	}

	for (size_t i = 1; i <= vsArgumentsList.size() - 1; i++) {
		if (vsArgumentsList[i] == "-h" || vsArgumentsList[i] == "--help") {
			Config::printHelp();
			return false;
		} else if (vsArgumentsList[i] == "-i"
				|| vsArgumentsList[i] == "--input") {
			string sInputFilenames = vsArgumentsList[++i];
			stringstream ss(sInputFilenames);
			string sItem;
			while (getline(ss, sItem, ',')) {
				Config::vsInputVariantFiles.push_back(sItem);
			}
		} else if (vsArgumentsList[i] == "-o"
				|| vsArgumentsList[i] == "--out") {
			string outputFilename = vsArgumentsList[++i];
			Config::sOutputFolder = outputFilename;
		} else {
			Config::printHelp();
			return false;
		}
	}

	struct stat info;
	if (stat(sOutputFolder.c_str(), &info) != 0) {
		printf("cannot access %s\n", sOutputFolder.c_str());
		return false;
	} else if ((info.st_mode & S_IFDIR) == 0) { // S_ISDIR() doesn't exist on my windows
		printf("%s is no directory\n", sOutputFolder.c_str());
		return false;
	}

	return true;
}
