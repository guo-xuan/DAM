//============================================================================
// Name        : DAM.cpp
// Author      : Xuan Guo
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <omp.h>

//DAM header
#include "Config.h"

using namespace std;


int main(int argc, char **argv) {

	// parse command line options:
	if (!Config::setConfig(argc, argv)) {
		cout << "Please follow the above help information." << endl;
		return false;
	}

	return 0;
}
