
#ifndef __MULTIALIGN_H
#define __MULTIALIGN_H


#include "PairAlign.h"


class MultiAlign {
	ProteinChain **m_chain, m_consensus;
	PairAlign *m_pair_align;
	int m_chain_num;
	double m_rmsd;
	int m_align_num;
	double m_lambda, m_annealing_initial, m_annealing_rate;
	bool m_sequential_order, m_annealing;
	int m_heuristic_start;

public:
	MultiAlign(int chain_num = 0);
	~MultiAlign();

	void setChainNum(int i);
	void setChain(int i, ProteinChain *chain);
	void setLambda(double lambda) { m_lambda = lambda; }
	void setSequentialOrder(bool enable) { m_sequential_order = enable; }
	void setHeuristicStart(int start) { m_heuristic_start = start; }
	void setAnnealing(bool enable) { m_annealing = enable; }
	void setAnnealingInitial(double annealing_initial) { m_annealing_initial = annealing_initial; }
	void setAnnealingRate(double annealing_rate) { m_annealing_rate = annealing_rate; }

	void align();

	void updateConsensus();

	void writePDBFile(const char *filename);

protected:
	void _modifyConsensus();
};


#endif // __MULTIALIGN_H
