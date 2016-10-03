/*
 * IndependentModule.cpp
 *
 *  Created on: Aug 9, 2016
 *      Author: naux
 */

#include "IndependentModule.h"

IndependentModule::IndependentModule() {
	vdIndependentPosteriorPerVariant = NULL;
	iNumAssociationTypes = 0;
	iNumVariants = 0;
}

IndependentModule::~IndependentModule() {
	if(vdIndependentPosteriorPerVariant != NULL) {
		delete[] vdIndependentPosteriorPerVariant;
	}
}

double IndependentModule::getIndependentPosterior(UINT32 _iVariantIndex, UINT32 _iAssociationType) {
	if(_iVariantIndex < iNumVariants && _iAssociationType < iNumAssociationTypes) {
		return vdIndependentPosteriorPerVariant[_iVariantIndex * iNumAssociationTypes + _iAssociationType];
	} else {
		cout << "Independent Posterior Out of Range." << endl;
		return 0.0;
	}
}

void IndependentModule::initial(GwasData * _pGwasData, vector<vector<vector<int>>>& vvviHyperGroup) {
	iNumAssociationTypes = vvviHyperGroup.size();
	iNumVariants = _pGwasData->iNumVariants;
	vdIndependentPosteriorPerVariant = new double[iNumVariants*iNumAssociationTypes];
#pragma omp parallel
		{
			vector<vector<vector<double>>> * vvviHyperGroupFrequency;
			vector<vector<vector<UINT32>>> * vvviHyperGroupOccurrence;
			double posterior = 0;
			double ni = 0, alphai = 0, N = 0, alpha = 0;
			UINT32 iVariantNumTypes = 0;
#pragma omp for schedule(guided)
			for(UINT32 i = 0;i<iNumVariants;i++) {
				iVariantNumTypes = _pGwasData->getVariantNumTypes(i);
				_pGwasData->setVariantHyperFrequency(i, vvviHyperGroup);
				vvviHyperGroupFrequency = _pGwasData->getVariantHyperGroupFrequency(i);
				vvviHyperGroupOccurrence = _pGwasData->getVariantHyperGroupOccurrence(i);
				for(UINT32 j=0;j<iNumAssociationTypes;j++) {
					posterior = 0;
					for(size_t k=0;k<vvviHyperGroup.at(j).size();k++) {
						alpha = 1;
						N = 0;
						for(UINT32 g=0;g<iVariantNumTypes;g++) {
							ni = vvviHyperGroupOccurrence->at(j).at(k).at(g);
							alphai = vvviHyperGroupFrequency->at(j).at(k).at(g);
							posterior += (lgamma(ni + alphai) - lgamma(alphai));
							N += ni;
						}
						posterior += (lgamma(alpha) - lgamma(N+alpha));
					}
					vdIndependentPosteriorPerVariant[i*iNumAssociationTypes + j] = posterior;
				}
			}
		}
	}
