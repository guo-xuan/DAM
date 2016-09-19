/*
 * IndependentModule.h
 *
 *  Created on: Aug 9, 2016
 *      Author: naux
 */

#ifndef INDEPENDENTMODULE_H_
#define INDEPENDENTMODULE_H_

#include "Config.h"
#include "GwasData.h"

/**
 * This module has the following information:
 * 1. frequency of each genotype of each variant in each group of samples
 * 2. posterior probability of each variant with different association types
 */
class IndependentModule {
public:

	// row -> variant, column -> association type, element -> posterior
	double * vdIndependentPosteriorPerVariant;

	UINT32 iNumAssociationTypes;
	UINT32 iNumVariants;

	IndependentModule();
	~IndependentModule();

	double getIndependentPosterior(UINT32 _iVariantIndex, UINT32 _iAssociationType);
	void initial(GwasData * _pGwasData, vector<vector<vector<int>>> & vvviHyperGroup);
};

#endif /* INDEPENDENTMODULE_H_ */
