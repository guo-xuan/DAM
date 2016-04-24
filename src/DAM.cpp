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

	int nProcessors = omp_get_max_threads();

	std::cout << nProcessors << std::endl;

	omp_set_num_threads(nProcessors);
	std::cout << omp_get_num_threads() << std::endl;

#pragma omp parallel for
	for (int i = 0; i < 5; i++) {
		int tid = omp_get_thread_num();
		std::cout << tid << "\t tid" << std::endl;
		int nThreads = omp_get_num_threads();
		std::cout << nThreads << "\t nThreads" << std::endl;
	}

	return 0;
}
