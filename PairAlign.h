
#ifndef __PAIRALIGN_H
#define __PAIRALIGN_H


#include "AlignParams.h"
#include "ProteinChain.h"


class PairAlign {
	ProteinChain *m_chain_a, *m_chain_b;
	int m_length_a, m_length_b;

	vector<int> m_alignment;
	double m_translation[3], m_rotation[3][3];
	int m_align_num, m_break_num, m_permu_num;
	double m_rmsd, m_score, m_sequence_identity;

	AlignParams m_params;

	vector<vector<double> > m_weights;

public:
	PairAlign(ProteinChain *chain_a = NULL, ProteinChain *chain_b = NULL);

	double rmsd() const { return m_rmsd; }
	int align_num() const { return m_align_num; }
	int alignment(int i) const { return m_alignment[i]; }

	void setChain(int i, ProteinChain *chain);
	void setParams(const AlignParams &params) { m_params = params; }

	double align();
	double alignBNB();
	double alignITER();

	double continueAlign();
	void postProcess();
	double postAlign(bool seq_order = false);

	bool getStart(int index, vector<int> &alignment);
	void initWeights();

	double evaluate(const string &filename);
	double improve(const string &filename);

	void setSolution(const double translation[3], const double rotation[3][3], const vector<int> &alignment);

	bool solveLeastSquare(double translation[3], double rotation[3][3], vector<int> &alignment);
	double solveMaxMatch(double translation[3], double rotation[3][3], vector<int> &alignment, double lambda);
	double solveMaxAlign(double translation[3], double rotation[3][3], vector<int> &alignment, double lambda);

	void writePDBFile(const string &filename) const;
	void writeSolutionFile(const string &filename) const;

private:
	int _getAlignNum(const vector<int> &alignment);
	int _getBreakNum(const vector<int> &alignment);
	int _getPermuNum(const vector<int> &alignment);
	double _getSequenceIdentity(const vector<int> &alignment);

	friend class MultiAlign;
};


#endif // __PAIRALIGN_H
