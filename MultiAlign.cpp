
#include "MultiAlign.h"


MultiAlign::MultiAlign(int chain_num)
{
	m_chain_num = chain_num;
	if (m_chain_num > 0) {
		m_chain = new ProteinChain * [m_chain_num];
		m_pair_align = new PairAlign [m_chain_num];
	}
	else {
		m_chain = NULL;
		m_pair_align = NULL;
	}
}

MultiAlign::~MultiAlign()
{
	delete[] m_pair_align;
	delete[] m_chain;
}

void MultiAlign::setChainNum(int i)
{
	if (i != m_chain_num) {
		m_chain_num = i;
		delete[] m_pair_align;
		delete[] m_chain;
		m_chain = new ProteinChain * [m_chain_num];
		m_pair_align = new PairAlign [m_chain_num];
	}
}

void MultiAlign::setChain(int i, ProteinChain *chain)
{
	m_chain[i] = chain;
}

void MultiAlign::align()
{
	int i, n, max_len;
	max_len = 0;
	n = 0;
	for (i=0; i<m_chain_num; i++) {
		if (m_chain[i]->length() > max_len) {
			max_len = m_chain[i]->length();
			n = i;
		}
	}
	m_consensus = *m_chain[n];
	m_align_num = 0;
	m_rmsd = 0;
	for (i=0; i<m_chain_num; i++) {
		m_pair_align[i].setChain(0, m_chain[i]);
		m_pair_align[i].setChain(1, &m_consensus);
		m_pair_align[i].setLambda(m_lambda);
		m_pair_align[i].setSequentialOrder(m_sequential_order);
		m_pair_align[i].setHeuristicStart(m_heuristic_start);
		m_pair_align[i].setAnnealing(m_annealing);
		m_pair_align[i].setAnnealingInitial(m_annealing_initial);
		m_pair_align[i].setAnnealingRate(m_annealing_rate);
		m_pair_align[i].align();
		m_align_num += m_pair_align[i].align_num();
		m_rmsd += m_pair_align[i].rmsd();
	}
	m_align_num /= m_chain_num;
	m_rmsd /= m_chain_num;
	Logger::info("Multiple Aligned: %d, RMSD: %f\n", m_align_num, m_rmsd);
	for (n=0; n<5; n++) {
		updateConsensus();
		m_align_num = 0;
		m_rmsd = 0;
		for (i=0; i<m_chain_num; i++) {
			m_pair_align[i].continueAlign();
			m_align_num += m_pair_align[i].align_num();
			m_rmsd += m_pair_align[i].rmsd();
		}
		m_align_num /= m_chain_num;
		m_rmsd /= m_chain_num;
		Logger::info("Multiple Aligned: %d, RMSD: %f\n", m_align_num, m_rmsd);
	}
}

void MultiAlign::updateConsensus()
{
	int i, j, k;
	double **matrix;
	int *alignment_number;

	alignment_number = new int [m_consensus.length()];
	for (i=0; i<m_consensus.length(); i++) {
		m_consensus[i][0] = 0;
		m_consensus[i][1] = 0;
		m_consensus[i][2] = 0;
		alignment_number[i] = 0;
	}
	for (i=0; i<m_chain_num; i++) {
		matrix = m_chain[i]->getMatrix(m_pair_align[i].m_translation, m_pair_align[i].m_rotation);
		for (j=0; j<m_chain[i]->length(); j++) {
			k = m_pair_align[i].alignment(j);
			if (k >= 0 && k < m_consensus.length()) {
				m_consensus[k][0] += matrix[j][0];
				m_consensus[k][1] += matrix[j][1];
				m_consensus[k][2] += matrix[j][2];
				alignment_number[k] ++;
			}
		}
		Matrix<double>::free(matrix);
	}
	for (i=0; i<m_consensus.length(); i++) {
		if (alignment_number[i] > 1) {
			m_consensus[i][0] /= alignment_number[i];
			m_consensus[i][1] /= alignment_number[i];
			m_consensus[i][2] /= alignment_number[i];
		}
	}
	delete[] alignment_number;
}

void MultiAlign::writePDBFile(const char *filename)
{
	FILE *fp;
	char buffer[81];
	int i, n;

	if ((fp = fopen(filename, "w")) == NULL) {
		Logger::error("Can not open the file: %s\n", filename);
		exit(1);
	}

	fprintf(fp, "HEADER    %-40s%30c\n", "MULTIPLE PROTEIN STRUCTURE ALIGNMENT", ' ');

	strcpy(buffer, "ALIGNMENT RESULT OF ");
	i = 0;
	n = 1;
	do {
		strcat(buffer, m_chain[i]->name());
		if (i != m_chain_num-1) {
			strcat(buffer, ", ");
		}
		if (i%8 == 4 || i == m_chain_num-1) {
			if (n == 1) {
				fprintf(fp, "TITLE     %-60s%10c\n", buffer, ' ');
			}
			else {
				fprintf(fp, "TITLE   %2d%-60s%10c\n", n, buffer, ' ');
			}
			strcpy(buffer, " ");
			n++;
		}
	} while (i++ < m_chain_num-1);

	for (i=0; i<m_chain_num; i++) {
		m_chain[i]->writePDBModel(fp, i+1, true, m_pair_align[i].m_translation, m_pair_align[i].m_rotation);
	}

	fprintf(fp, "END   %74c\n", ' ');

	fclose(fp);
}
