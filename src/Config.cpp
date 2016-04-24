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
		} else if (vsArgumentsList[i] == "-i" || vsArgumentsList[i] == "--input") {
			string sInputFilenames = vsArgumentsList[++i];
			stringstream ss(sInputFilenames);
			string sItem;
			while (getline(ss, sItem, ',')) {
				Config::vsInputVariantFiles.push_back(sItem);
			}
		} else if (vsArgumentsList[i] == "-o" || vsArgumentsList[i] == "--out") {
			string outputFilename = vsArgumentsList[++i];
			Config::sOutputFolder = outputFilename;
		} else {
			Config::printHelp();
			return false;
		}
	}

	return true;
}
