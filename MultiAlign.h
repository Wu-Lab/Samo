
#ifndef __MULTIALIGN_H
#define __MULTIALIGN_H


#include "PairAlign.h"


class MultiAlign {
	ProteinChain **m_chain, m_consensus;
	PairAlign *m_pair_align;
	int m_chain_num;
	double m_rmsd;
	int m_align_num;
	AlignParams m_params;

public:
	MultiAlign(int chain_num = 0);
	~MultiAlign();

	void setChainNum(int i);
	void setChain(int i, ProteinChain *chain);
	void setParams(const AlignParams &params) { m_params = params; }

	void align();

	void updateConsensus();

	void writePDBFile(const string &filename) const;
	void writeSolutionFile(const string &filename) const;

protected:
	void _modifyConsensus();
};


#endif // __MULTIALIGN_H
