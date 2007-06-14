
#include <map>
#include <numeric>

#include "Utils.h"
#include "Matrix.h"
#include "PairAlign.h"
#include "SVD.h"
#include "FibHeap.h"

#define square(d) ((d)*(d))
#define get_square_distance(x, y) (square((x)[0]-(y)[0])+square((x)[1]-(y)[1])+square((x)[2]-(y)[2]))
#define get_distance(x, y) sqrt(get_square_distance(x, y))
#define get_outer_distance(x, y, i, j) get_distance((x)[i], (y)[j])
#define get_inner_distance(x, i, j) get_distance((x)[i], (x)[j])

PairAlign::PairAlign(ProteinChain *chain_a, ProteinChain *chain_b)
{
	int i, j;
	m_chain_a = chain_a;
	m_chain_b = chain_b;
	for (i=0; i<3; i++) {
		m_translation[i] = 0;
		for (j=0; j<3; j++) {
			m_rotation[i][j] = 0;
		}
		m_rotation[i][i] = 1.0;
	}
	if (m_chain_a != NULL) {
		m_length_a = m_chain_a->length();
		m_alignment.resize(m_length_a, -1);
	}
	if (m_chain_b != NULL) {
		m_length_b = m_chain_b->length();
	}
	m_lambda = 6.0;
	m_branch_and_bound = false;
	m_sequential_order = false;
	m_annealing = false;
	m_annealing_initial = 60.0;
	m_annealing_rate = 0.4;
	m_heuristic_start = 0;
	m_min_fragment_length = 8;
	m_fragment_threshold0 = 3;
	m_fragment_threshold1 = 4;
}

void PairAlign::setChain(int i, ProteinChain *chain)
{
	if (i == 0) {
		m_chain_a = chain;
		if (m_chain_a != NULL && m_length_a != m_chain_a->length()) {
			m_length_a = m_chain_a->length();
			m_alignment.resize(m_length_a);
		}
	}
	else {
		m_chain_b = chain;
		if (m_chain_b != NULL && m_length_b != m_chain_b->length()) {
			m_length_b = m_chain_b->length();
		}
	}
}

double PairAlign::align()
{
	m_align_num = 0;
	m_break_num = 0;
	m_permu_num = 0;
	m_sequence_identity = 0;
	m_rmsd = 0;
	if (m_length_a == 0 || m_length_b == 0)
	{
		Logger::warning("Attmpt to align empty chain!");
		return m_rmsd;
	}
	if (m_branch_and_bound) {
		alignBNB();
	}
	else {
		alignITER();
	}

// 	m_chain_a->writeChainCode(stderr);
// 	Logger::info("");
// 
// 	for (int i=0; i<m_length_a; i++) {
// 		if (m_alignment[i] >= 0) fprintf(stderr, "%c", (*m_chain_b)[m_alignment[i]].resCode());
// 		else fprintf(stderr, "%c", '.');
// 	}
// 	Logger::info("");
// 
// 	m_chain_b->writeChainCode(stderr);
// 	Logger::info("");

	return m_rmsd;
}

double PairAlign::alignBNB()
{
	double translation[3], rotation[3][3], lambda2, rmsd, score;
	int align_num;
	vector<int> alignment(m_length_a);
	vector<int> occupied(m_length_b);

	int index, i;

	lambda2 = m_lambda * m_lambda;

	alignITER();
	
	for (i=0; i<m_length_a; i++) {
		alignment[i] = -1;
	}
	for (i=0; i<m_length_b; i++) {
		occupied[i] = false;
	}

	index = 0;
	while (index >= 0) {
		i = alignment[index];
		if (i < -1) {
			// backtracking
			alignment[index] = -1;
			index--;
			continue;
		}
		else if (i >=0 && i < m_length_b) {
			occupied[i] = false;
		}
		do {
			i++;
		} while(i >= 0 && i < m_length_b && occupied[i]);
		if (i >= m_length_b) {
			i = -2;
		}
		alignment[index] = i;
		if (i >= 0 && i < m_length_b) {
			occupied[i] = true;
			// test branch
			align_num = _getAlignNum(alignment);
			if (align_num >= 1 && solveLeastSquare(translation, rotation, alignment)) {
				rmsd = m_chain_a->getRMSD(*m_chain_b, translation, rotation, alignment);
				score = (rmsd * rmsd - lambda2) * align_num;
//				Logger::debug("\tAligned: %d, RMSD: %f, Score: %f, Best: %f", align_num, rmsd, score, best_score);
				if (score < m_score) {
					// incumbent, candidate of best solution
					m_score = score;
					m_align_num = align_num;
					m_rmsd = rmsd;
					setSolution(translation, rotation, alignment);
					Logger::debug("\tScore: %f, Aligned: %d, RMSD: %f", score, align_num, rmsd);
				}
				else if (score-lambda2*(m_length_a-index-1) >= m_score) {
					// bad branch, pruning
					continue;
				}
			}
		}
		if (index < m_length_a-1) {
			// go on, next position
			index++;
			continue;
		}
		else {
			// last position, next value
			continue;
		}
	}

	return m_rmsd;
}

double PairAlign::alignITER()
{
	double translation[3], rotation[3][3], lambda, rmsd;
	int align_num;
	double score_old, score_new;
	int start_index;
	vector<int> alignment(m_length_a);
	m_score = HUGE_VAL;
	start_index = 0;
	while (getStart(start_index++, alignment)) {
		solveLeastSquare(translation, rotation, alignment);
		Logger::info("\tInitial solution: %f", m_chain_a->getRMSD(*m_chain_b, translation, rotation, alignment));
		if (m_annealing) {
			lambda = m_annealing_initial;
			do {
				score_new = HUGE_VAL;
				do {
					score_old = score_new;
					solveLeastSquare(translation, rotation, alignment);
					score_new = solveMaxMatch(translation, rotation, alignment, m_lambda+lambda);
					Logger::debug("\t%f, %f", m_lambda+lambda, score_new);
					if (score_new > score_old) {
						Logger::warning("Not convergent!");
						break;
					}
				} while (fabs(score_new - score_old) > 0.01);
				lambda *= m_annealing_rate;
			} while(lambda > 0.01);
		}
		else {
			score_new = HUGE_VAL;
			do {
				score_old = score_new;
				solveLeastSquare(translation, rotation, alignment);
				score_new = solveMaxMatch(translation, rotation, alignment, m_lambda);
				Logger::debug("\t%f", score_new);
				if (score_new > score_old) {
					Logger::warning("Not convergent!");
					break;
				}
			} while (fabs(score_new - score_old) > 0.01);
		}
		score_new = solveMaxMatch(translation, rotation, alignment, m_lambda);
		align_num = _getAlignNum(alignment);
		rmsd = m_chain_a->getRMSD(*m_chain_b, translation, rotation, alignment);
		Logger::info("\tScore: %f, Aligned: %d, RMSD: %f", score_new, align_num, rmsd);
		Logger::info("===============================================================================");
		if (score_new < m_score) {
			m_score = score_new;
			m_align_num = align_num;
			m_rmsd = rmsd;
			setSolution(translation, rotation, alignment);
		}
	}

	return m_rmsd;
}

// using distance matrix to setup initial solutions

double PairAlign::alignITER_dmstart()
{
	int i, j, k, l, m, n;
	double s;
	int length_a, length_b;
	double **inner_distance_a, **inner_distance_b, **similarity;
	int **gapless_aligned_fragment;

	multimap<int, int, greater<int> > candidates;
	multimap<int, int, greater<int> >::iterator ci;
	double translation[3], rotation[3][3], lambda, rmsd;
	int align_num;
	double score_old, score_new;
	vector<int> alignment(m_length_a);

	length_a = m_length_a - m_min_fragment_length + 1;
	length_b = m_length_b - m_min_fragment_length + 1;
	inner_distance_a = Matrix<double>::alloc(m_length_a, m_length_a);
	inner_distance_b = Matrix<double>::alloc(m_length_b, m_length_b);
	similarity = Matrix<double>::alloc(length_a, length_b);
	gapless_aligned_fragment = Matrix<int>::alloc(length_a, length_b);

	for (i=0; i<m_length_a; i++) {
		for (j=0; j<=i; j++) {
			inner_distance_a[i][j] = get_inner_distance(*m_chain_a, i, j);
			inner_distance_a[j][i] = inner_distance_a[i][j];
		}
	}
	for (i=0; i<m_length_b; i++) {
		for (j=0; j<=i; j++) {
			inner_distance_b[i][j] = get_inner_distance(*m_chain_b, i, j);
			inner_distance_b[j][i] = inner_distance_b[i][j];
		}
	}
	for (i=0; i<length_a; i++) {
		for (j=0; j<length_b; j++) {
			s = 0;
			for (k=0; k<m_min_fragment_length; k++) {
				for (l=0; l<m_min_fragment_length; l++) {
					s += fabs(inner_distance_a[i+k][i+l] - inner_distance_b[j+k][j+l]);
				}
			}
			similarity[i][j] = s / (m_min_fragment_length * m_min_fragment_length);
			if (similarity[i][j] < m_fragment_threshold0) {
				gapless_aligned_fragment[i][j] = m_min_fragment_length;
			}
			else {
				gapless_aligned_fragment[i][j] = 0;
			}
		}
	}
	for (i=0; i<length_a; i++) {
		for (j=0; j<length_b; j++) {
			if (gapless_aligned_fragment[i][j] > 0) {
				m = min(length_a-1-i, length_b-1-j);
				for (k=1; k<m; k++) {
					if (similarity[i+k][j+k] < m_fragment_threshold0) {
						n = m_min_fragment_length - 1 + k;
						s = 0;
						for (l=0; l<n; l++) {
							s += fabs(inner_distance_a[i+l][i+n] - inner_distance_b[j+l][j+n]);
						}
						if (s/n < m_fragment_threshold1) {
							gapless_aligned_fragment[i+k][j+k] = 0;
// 							s = (similarity[i][j] * n * n + s * 2);
// 							n++;
// 							similarity[i][j] = s / (n * n);
							gapless_aligned_fragment[i][j]++;
							continue;
						}
					}
					break;
				}
			}
		}
	}

	m = length_a+length_b;
	for (i=0; i<length_a; i++) {
		for (j=0; j<length_b; j++) {
			if (gapless_aligned_fragment[i][j] > 0) {
				candidates.insert(make_pair(gapless_aligned_fragment[i][j], i*m+j));
			}
		}
	}

	Logger::info("candidates = %d", candidates.size());

	m_score = HUGE_VAL;

	ci = candidates.begin();
	for (l=0; l<10; l++) {
		if (ci == candidates.end()) break;

		k = ci->second;
		i = k / m;
		j = k % m;

		Logger::info("%d, (%d, %d) %d", k, i, j, gapless_aligned_fragment[i][j]);

		for (k=0; k<m_length_a; k++) {
			alignment[k] = -1;
		}
		for (k=0; k<gapless_aligned_fragment[i][j]; k++) {
			alignment[i+k] = j+k;
		}

		solveLeastSquare(translation, rotation, alignment);
		Logger::info("\tInitial solution: %f", m_chain_a->getRMSD(*m_chain_b, translation, rotation, alignment));
		if (m_annealing) {
			lambda = m_annealing_initial;
			do {
				score_new = HUGE_VAL;
				do {
					score_old = score_new;
					solveLeastSquare(translation, rotation, alignment);
					score_new = solveMaxMatch(translation, rotation, alignment, m_lambda+lambda);
					Logger::debug("\t%f, %f", m_lambda+lambda, score_new);
					if (score_new > score_old) {
						Logger::warning("Not convergent!");
 						break;
					}
				} while (fabs(score_new - score_old) > 0.01);
				lambda *= m_annealing_rate;
			} while(lambda > 0.01);
		}
		else {
			score_new = HUGE_VAL;
			do {
				score_old = score_new;
				solveLeastSquare(translation, rotation, alignment);
				score_new = solveMaxMatch(translation, rotation, alignment, m_lambda);
				Logger::debug("\t%f", score_new);
				if (score_new > score_old) {
					Logger::warning("Not convergent!");
 					break;
				}
			} while (fabs(score_new - score_old) > 0.01);
		}
		score_new = solveMaxMatch(translation, rotation, alignment, m_lambda);
		align_num = _getAlignNum(alignment);
		rmsd = m_chain_a->getRMSD(*m_chain_b, translation, rotation, alignment);
		Logger::info("\tScore: %f, Aligned: %d, RMSD: %f", score_new, align_num, rmsd);
		Logger::info("===============================================================================");
		if (score_new < m_score) {
			m_score = score_new;
			m_align_num = align_num;
			m_rmsd = rmsd;
			setSolution(translation, rotation, alignment);
		}

	}

	Matrix<int>::free(gapless_aligned_fragment);
	Matrix<double>::free(similarity);
	Matrix<double>::free(inner_distance_a);
	Matrix<double>::free(inner_distance_b);
	return m_rmsd;
}

double PairAlign::continueAlign()
{
	double score_old, score_new;
	score_new = HUGE_VAL;
	do {
		score_old = score_new;
		solveLeastSquare(m_translation, m_rotation, m_alignment);
		score_new = solveMaxMatch(m_translation, m_rotation, m_alignment, m_lambda);
	} while (fabs(score_new - score_old) > 0.01);
	m_align_num = _getAlignNum(m_alignment);
	m_rmsd = m_chain_a->getRMSD(*m_chain_b, m_translation, m_rotation, m_alignment);
	return m_rmsd;
}

void PairAlign::postProcess()
{
	m_break_num = _getBreakNum(m_alignment);
	m_permu_num = _getPermuNum(m_alignment);
	Logger::info("PairAlign: %s (size=%d) vs %s (size=%d)\n\tAligned = %d, RMSD = %f\n\tBreak/Permutation = %d/%d, SeqId = %5.3f",
		m_chain_a->raw_name(), m_length_a, m_chain_b->raw_name(), m_length_b,
		m_align_num, m_rmsd, m_break_num, m_permu_num, _getSequenceIdentity(m_alignment));

	if (m_sequential_order) {
		postAlignWithSequentialOrder();
		Logger::info("PostAlign: %s (size=%d) vs %s (size=%d)\n\tAligned = %d, RMSD = %f\n\tBreak/Permutation = %d/%d, SeqId = %5.3f",
			m_chain_a->raw_name(), m_length_a, m_chain_b->raw_name(), m_length_b,
			m_align_num, m_rmsd, m_break_num, m_permu_num, _getSequenceIdentity(m_alignment));
	}
}

double PairAlign::postAlignWithSequentialOrder()
{
	solveMaxAlign(m_translation, m_rotation, m_alignment, m_lambda);
	m_align_num = _getAlignNum(m_alignment);
	m_rmsd = m_chain_a->getRMSD(*m_chain_b, m_translation, m_rotation, m_alignment);
	m_break_num = _getBreakNum(m_alignment);
	m_permu_num = _getPermuNum(m_alignment);
	return m_rmsd;
}

bool PairAlign::getStart(int index, vector<int> &alignment)
{
	int i, j, k;
	int block_length, block_number_a, block_number_b;
	int block_a, block_b;

	for (k=0; k<m_length_a; k++) {
		alignment[k] = -1;
	}
	if (m_heuristic_start > 0) {
		block_length = (max(m_length_a, m_length_b) + m_heuristic_start - 1) / m_heuristic_start;
		block_number_a = (m_length_a + block_length - 1) / block_length;
		block_number_b = (m_length_b + block_length - 1) / block_length;

		block_a = index / block_number_b;
		block_b = index % block_number_b;

		if (block_a < block_number_a) {
			for (k=0; k<block_length; k++) {
				i = block_a * block_length + k;
				j = block_b * block_length + k;
				if (i < m_length_a && j < m_length_b) {
					alignment[i] = j;
				}
			}
			return true;
		}
	}
	else {
		if (index == 0) {
			for (k=0; k<m_length_a; k++) {
				alignment[k] = k;
			}
			return true;
		}
	}
	return false;
}

double PairAlign::evaluate(const string &filename)
{
	FILE *fp;
	char buffer[100];
	int i, j;
	bool has_alignment, has_tranlation, has_rotation;

	has_alignment = false;
	has_tranlation = false;
	has_rotation = false;

	for (i=0; i<m_length_a; i++) {
		m_alignment[i] = -1;
	}
	for (i=0; i<3; i++) {
		m_translation[i] = 0;
		for (j=0; j<3; j++) {
			m_rotation[i][j] = 0;
		}
		m_rotation[i][i] = 1.0;
	}

	if ((fp = fopen(filename.c_str(), "r")) == NULL) {
		Logger::error("Can not open the file: %s\n", filename.c_str());
		exit(1);
	}

	while (!feof(fp)) {
		fgets(buffer, 100, fp);

		if (strncmp(buffer, "[Alignment]", 11) == 0) {
			for (i=0; i<m_length_a; i++) {
				fscanf(fp, "%d", &m_alignment[i]);
			}
			has_alignment = true;
		}
		else if (strncmp(buffer, "[Translation]", 13) == 0) {
			for (i=0; i<3; i++) {
				fscanf(fp, "%lf", &m_translation[i]);
			}
			has_tranlation = true;
		}
		else if (strncmp(buffer, "[Rotation]", 10) == 0) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					fscanf(fp, "%lf", &m_rotation[i][j]);
				}
			}
			has_rotation = true;
		}
	}
	
	fclose(fp);

	if (has_alignment) {
		solveLeastSquare(m_translation, m_rotation, m_alignment);
	}
	else if (has_tranlation && has_rotation) {
		solveMaxMatch(m_translation, m_rotation, m_alignment, m_lambda);
	}
	else if (has_tranlation) {
		Logger::error("Rotation matrix is missing!");
		exit(1);
	}
	else {
		Logger::error("Translation matrix is missing!");
		exit(1);
	}

	m_align_num = _getAlignNum(m_alignment);
	m_rmsd = m_chain_a->getRMSD(*m_chain_b, m_translation, m_rotation, m_alignment);
	Logger::debug("\tEvaluated: %d, RMSD: %f", m_align_num, m_rmsd);

	return m_rmsd;
}

double PairAlign::improve(const string &filename)
{
	double score_old, score_new;

	evaluate(filename);
	score_new = HUGE_VAL;
	do {
		score_old = score_new;
		solveLeastSquare(m_translation, m_rotation, m_alignment);
		score_new = solveMaxMatch(m_translation, m_rotation, m_alignment, m_lambda);
		Logger::debug("\t%f", score_new);
		if (score_new > score_old) {
			Logger::warning("Not convergent!");
			break;
		}
	} while (fabs(score_new - score_old) > 0.01);
	solveMaxMatch(m_translation, m_rotation, m_alignment, m_lambda);

	m_align_num = _getAlignNum(m_alignment);
	m_rmsd = m_chain_a->getRMSD(*m_chain_b, m_translation, m_rotation, m_alignment);
	Logger::debug("\tImproved: %d, RMSD: %f", m_align_num, m_rmsd);

	return m_rmsd;
}

void PairAlign::setSolution(const double translation[3], const double rotation[3][3], const vector<int> &alignment)
{
	int i, j;
	for (i=0; i<3; i++) {
		m_translation[i] = translation[i];
		for (j=0; j<3; j++) {
			m_rotation[i][j] = rotation[i][j];
		}
	}
	for (i=0; i<m_length_a; i++) {
		m_alignment[i] = alignment[i];
	}
}

bool PairAlign::solveLeastSquare(double translation[3], double rotation[3][3], vector<int> &alignment)
{
	double weight, center_a[3], center_b[3], min_sv_val;
	double **matrix_u, **matrix_v, vector_d[4], det;
	bool success;
	int min_sv_ind, i, j, k, l;
// 	double w;

	matrix_u = Matrix<double>::alloc(1, 3, 1, 3);
	matrix_v = Matrix<double>::alloc(1, 3, 1, 3);

	for (i=0; i<3; i++) {
		center_a[i] = 0;
		center_b[i] = 0;
	}
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			matrix_u[i+1][j+1] = 0;
		}
	}

// 	if (m_localization) {
// 		weight = 0;
// 		for (k=0; k<m_length_a; k++) {
// 			l = alignment[k];
// 			if (l >= 0) {
// 				w = m_local_weight_a[k];
// 				for (i=0; i<3; i++) {
// 					center_a[i] += (*m_chain_a)[k][i] * w;
// 					center_b[i] += (*m_chain_b)[l][i] * w;
// 				}
// 				for (i=0; i<3; i++) {
// 					for (j=0; j<3; j++) {
// 						matrix_u[i+1][j+1] += (*m_chain_a)[k][i] * (*m_chain_b)[l][j] * w;
// 					}
// 				}
// 				weight += w;
// 			}
// 		}
// 	}
// 	else {
		weight = _getAlignNum(alignment);
		for (k=0; k<m_length_a; k++) {
			l = alignment[k];
			if (l >= 0) {
				for (i=0; i<3; i++) {
					center_a[i] += (*m_chain_a)[k][i];
					center_b[i] += (*m_chain_b)[l][i];
				}
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						matrix_u[i+1][j+1] += (*m_chain_a)[k][i] * (*m_chain_b)[l][j];
					}
				}
			}
		}
// 	}

	for (i=0; i<3; i++) {
		center_a[i] /= weight;
		center_b[i] /= weight;
	}
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			matrix_u[i+1][j+1] = matrix_u[i+1][j+1] - center_a[i] * center_b[j] * weight;
		}
	}

	svdcmp(matrix_u, 3, 3, vector_d, matrix_v);

	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			rotation[i][j] = 0;
			for (k=0; k<3; k++) {
				rotation[i][j] += matrix_v[i+1][k+1] * matrix_u[j+1][k+1];
			}
		}
	}
	det = rotation[0][0] * rotation[1][1] * rotation[2][2];
	det += rotation[0][1] * rotation[1][2] * rotation[2][0];
	det += rotation[0][2] * rotation[1][0] * rotation[2][1];
	det -= rotation[0][2] * rotation[1][1] * rotation[2][0];
	det -= rotation[0][1] * rotation[1][0] * rotation[2][2];
	det -= rotation[0][0] * rotation[1][2] * rotation[2][1];

	success = true;
	if (det < 0) {
		min_sv_val = HUGE_VAL;
		for (i=1; i<=3; i++) {
			if (vector_d[i] < min_sv_val) {
				min_sv_val = vector_d[i];
				min_sv_ind = i;
			}
		}
		if (min_sv_val > 1e-6) {
			Logger::debug("Degenerate case, SVD method fails!");
			Logger::debug("\t\tD = (%f, %f, %f)", vector_d[1], vector_d[2], vector_d[3]);
			success = false;
		}
		for (i=1; i<=3; i++) {
			matrix_v[i][min_sv_ind] = - matrix_v[i][min_sv_ind];
		}
		for (i=0; i<3; i++) {
			for (j=0; j<3; j++) {
				rotation[i][j] = 0;
				for (k=0; k<3; k++) {
					rotation[i][j] += matrix_v[i+1][k+1] * matrix_u[j+1][k+1];
				}
			}
		}
	}
	
	for (i=0; i<3; i++) {
		translation[i] = center_b[i];
		for (j=0; j<3; j++) {
			translation[i] -= rotation[i][j] * center_a[j];
		}
	}

	Matrix<double>::free(matrix_u, 1, 1);
	Matrix<double>::free(matrix_v, 1, 1);
	return success;
}


// for using FibHeap

class HeapNode: public FibHeapNode
{
	double m_label;
	int m_index;

public:
	HeapNode() : FibHeapNode() { m_label = 0; };

	virtual void operator =(FibHeapNode& RHS);
	virtual int  operator ==(FibHeapNode& RHS);
	virtual int  operator <(FibHeapNode& RHS);

	virtual void operator =(double key);

	double getKeyValue() { return m_label; }
	void setKeyValue(double key) { m_label = key; }

	int getIndexValue() { return m_index; }
	void setIndexValue(int i) { m_index = i; }
};

void HeapNode::operator =(double label)
{
	HeapNode temp;
	temp.m_label = m_label = label;
	FHN_assign(temp);
}

void HeapNode::operator =(FibHeapNode& RHS)
{
	FHN_assign(RHS);
	m_label = ((HeapNode&) RHS).m_label;
}

int  HeapNode::operator ==(FibHeapNode& RHS)
{
	if (FHN_compare(RHS)) return 0;
	return m_label == ((HeapNode&) RHS).m_label ? 1 : 0;
}

int  HeapNode::operator <(FibHeapNode& RHS)
{
	int X;

	if ((X=FHN_compare(RHS)) != 0) return X < 0 ? 1 : 0;
	return m_label < ((HeapNode&) RHS).m_label ? 1 : 0;
}


double PairAlign::solveMaxMatch(double translation[3], double rotation[3][3], vector<int> &alignment, double lambda)
{
	double score, reduced_cost, **weight, label_t, current_value;
	bool *match_free, *label_free;
	int *active_num, **active_index, *backtrack, backtrack_t;
	HeapNode *label, *min_label, temp_label;
	FibHeap heap;
	double d;
	int i, j, k;

	weight = Matrix<double>::alloc(m_length_b+1, m_length_a);
	match_free = new bool [m_length_b+m_length_a];
	label_free = new bool [m_length_b+m_length_a];
	active_num = new int [m_length_b];
	active_index = Matrix<int>::alloc(m_length_b, m_length_a);
	backtrack = new int [m_length_b+m_length_a];
	label = new HeapNode [m_length_b+m_length_a];

	double **matrix = m_chain_a->getMatrix(translation, rotation);
	for (i=0; i<m_length_b; i++) {
		for (j=0; j<m_length_a; j++) {
			weight[i][j] = 0;
			for (k=0; k<3; k++) {
				d = matrix[j][k] - (*m_chain_b)[i][k];
				weight[i][j] += d * d;
			}
		}
	}
	Matrix<double>::free(matrix);

	lambda = lambda * lambda;
	for (i=0; i<m_length_b; i++) {
		active_num[i] = 0;
		for (j=0; j<m_length_a; j++) {
			weight[i][j] -= lambda;
			if (weight[i][j] <= 0) {
				active_index[i][active_num[i]] = j;
				active_num[i]++;
			}
// 			if (m_localization) {
// 				weight[i][j] *= m_local_weight_a[j];
// 			}
		}
	}
	
	for (i=0; i<m_length_b+m_length_a; i++) {
		label[i].setIndexValue(i);
		match_free[i] = true;
	}
	for (j=0; j<m_length_a; j++) {
		alignment[j] = -1;
		weight[m_length_b][j] = 0;
	}

	reduced_cost = 0;
	for (i=0; i<m_length_b; i++) {
		for (j=0; j<active_num[i]; j++) {
			k = active_index[i][j];
			if (weight[i][k] < reduced_cost) {
				reduced_cost = weight[i][k];
			}
		}
	}
	for (i=0; i<m_length_b; i++) {
		for (j=0; j<active_num[i]; j++) {
			k = active_index[i][j];
			weight[i][k] -= reduced_cost;
		}
	}
	score = 0;

	while (true) {
		// Dijkstra's algorithm
		for (i=0; i<m_length_b+m_length_a; i++) {
			if (i < m_length_b && match_free[i]) {
				label[i].setKeyValue(0);
			}
			else {
				label[i].setKeyValue(HUGE_VAL);
			}
			heap.insert(&label[i]);
			label_free[i] = true;
			backtrack[i] = -1;
		}
		label_t = HUGE_VAL;
		backtrack_t = -1;
		while (true) {
			min_label = (HeapNode *)heap.extractMin();
			if (min_label == NULL || min_label->getKeyValue() >= HUGE_VAL) break;
			i = min_label->getIndexValue();
			current_value = min_label->getKeyValue();
			label_free[i] = false;
			if (i < m_length_b) {
				for (j=0; j<active_num[i]; j++) {
					k = active_index[i][j];
					if (label_free[m_length_b+k] && alignment[k] != i) {
						d = current_value + weight[i][k];
						if (d < label[m_length_b+k].getKeyValue()) {
							temp_label = label[m_length_b+k];
							temp_label.setKeyValue(d);
							heap.decreaseKey(&label[m_length_b+k], temp_label);
							backtrack[m_length_b+k] = i;
						}
					}
				}
			}
			else {
				k = alignment[i-m_length_b];
				if (k >= 0 && label_free[k]) {
					d = current_value - weight[k][i-m_length_b];
					if (d < label[k].getKeyValue()) {
						temp_label = label[k];
						temp_label.setKeyValue(d);
						heap.decreaseKey(&label[k], temp_label);
						backtrack[k] = i;
					}
				}
				if (match_free[i]) {
					d = current_value + weight[m_length_b][i-m_length_b];
					if (d < label_t) {
						label_t = d;
						backtrack_t = i;
					}
				}
			}
		}
		while (min_label != NULL) {
			min_label = (HeapNode *)heap.extractMin();
		}
		reduced_cost += label_t;
		if (reduced_cost >= 0) break;
		score += reduced_cost;

		// Modify the weights
		for (i=0; i<m_length_b; i++) {
			for (j=0; j<active_num[i]; j++) {
				k = active_index[i][j];
				weight[i][k] += label[i].getKeyValue() - label[m_length_b+k].getKeyValue();
			}
		}
		for (j=0; j<m_length_a; j++) {
			weight[m_length_b][j] += label[m_length_b+j].getKeyValue() - label_t;
		}

		// Reverse path
		k = backtrack_t;
		while (k >= 0 && backtrack[k] >= 0) {
			alignment[k-m_length_b] = backtrack[k];
			match_free[k] = false;
			k = backtrack[k];
			match_free[k] = false;
			k = backtrack[k];
		}
	}

	delete[] label;
	delete[] backtrack;
	Matrix<int>::free(active_index);
	delete[] active_num;
	delete[] label_free;
	delete[] match_free;
	Matrix<double>::free(weight);

	return score;
}

double PairAlign::solveMaxAlign(double translation[3], double rotation[3][3], vector<int> &alignment, double lambda)
{
	double score, **weight, **score_matrix;
	int **backtrack_b, **backtrack_a;
	bool **align_matrix, match;
	double d, dd, da, db;
	int i, j, k;

	weight = Matrix<double>::alloc(m_length_b, m_length_a);
	score_matrix = Matrix<double>::alloc(m_length_b, m_length_a);
	backtrack_b = Matrix<int>::alloc(m_length_b, m_length_a);
	backtrack_a = Matrix<int>::alloc(m_length_b, m_length_a);
	align_matrix = Matrix<bool>::alloc(m_length_b, m_length_a);

	lambda = lambda * lambda;
	double **matrix = m_chain_a->getMatrix(translation, rotation);
	for (i=0; i<m_length_b; i++) {
		for (j=0; j<m_length_a; j++) {
			weight[i][j] = 0;
			for (k=0; k<3; k++) {
				d = matrix[j][k] - (*m_chain_b)[i][k];
				weight[i][j] += d * d;
			}
			weight[i][j] -= lambda;
		}
	}
	Matrix<double>::free(matrix);

	if (weight[0][0] <= 0) {
		score_matrix[0][0] = weight[0][0];
		align_matrix[0][0] = true;
	}
	else {
		score_matrix[0][0] = 0;
		align_matrix[0][0] = false;
	}
	backtrack_b[0][0] = -1;
	backtrack_a[0][0] = -1;
	for (i=1; i<m_length_b; i++) {
		if (weight[i][0] <= score_matrix[i-1][0]) {
			score_matrix[i][0] = weight[i][0];
			align_matrix[i][0] = true;
			backtrack_b[i][0] = -1;
			backtrack_a[i][0] = -1;
		}
		else {
			score_matrix[i][0] = score_matrix[i-1][0];
			align_matrix[i][0] = false;
			backtrack_b[i][0] = i-1;
			backtrack_a[i][0] = 0;
		}
	}
	for (j=1; j<m_length_a; j++) {
		if (weight[0][j] <= score_matrix[0][j]) {
			score_matrix[0][j] = weight[0][j];
			align_matrix[0][j] = true;
			backtrack_b[0][j] = -1;
			backtrack_a[0][j] = -1;
		}
		else {
			score_matrix[0][j] = score_matrix[0][j-1];
			align_matrix[0][j] = false;
			backtrack_b[0][j] = 0;
			backtrack_a[0][j] = j-1;
		}
	}

	for (i=1; i<m_length_b; i++) {
		for (j=1; j<m_length_a; j++) {
			if (weight[i][j] <= 0) {
				dd = score_matrix[i-1][j-1] + weight[i][j];
				match = true;
			}
			else {
				dd = score_matrix[i-1][j-1];
				match = false;
			}
			da = score_matrix[i-1][j];
			db = score_matrix[i][j-1];
			if (dd <= da && dd <= db) {
				score_matrix[i][j] = dd;
				align_matrix[i][j] = match;
				backtrack_b[i][j] = i-1;
				backtrack_a[i][j] = j-1;
			}
			else if (da <= db) {
				score_matrix[i][j] = da;
				align_matrix[i][j] = false;
				backtrack_b[i][j] = i-1;
				backtrack_a[i][j] = j;
			}
			else {
				score_matrix[i][j] = db;
				align_matrix[i][j] = false;
				backtrack_b[i][j] = i;
				backtrack_a[i][j] = j-1;
			}
		}
	}

	for (i=0; i<m_length_a; i++) {
		alignment[i] = -1;
	}
	i = m_length_b - 1;
	j = m_length_a - 1;
	score = score_matrix[i][j];
	while (i >= 0 && j >= 0) {
		if (align_matrix[i][j]) {
			alignment[j] = i;
		}
		k = backtrack_b[i][j];
		j = backtrack_a[i][j];
		i = k;
	}

	Matrix<bool>::free(align_matrix);
	Matrix<int>::free(backtrack_a);
	Matrix<int>::free(backtrack_b);
	Matrix<double>::free(score_matrix);
	Matrix<double>::free(weight);

	return score;
}

void PairAlign::writePDBFile(const string &filename) const
{
	FILE *fp;
	char buffer[80];

	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
		Logger::error("Can not open the file: %s\n", filename.c_str());
		exit(1);
	}

	fprintf(fp, "HEADER    %-40s%30c\n", "PAIRWISE PROTEIN STRUCTURE ALIGNMENT", ' ');

	strcpy(buffer, "ALIGNMENT RESULT OF ");
	strcat(buffer, m_chain_a->raw_name());
	strcat(buffer, " AND ");
	strcat(buffer, m_chain_b->raw_name());
	fprintf(fp, "TITLE     %-60s%10c\n", buffer, ' ');

	m_chain_a->writePDBModel(fp, 1, true, m_translation, m_rotation);
	m_chain_b->writePDBModel(fp, 2, true);

	fprintf(fp, "END   %74c\n", ' ');

	fclose(fp);
}

void PairAlign::writeSolutionFile(const string &filename) const
{
	FILE *fp;
	int i, j;

	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
		Logger::error("Can not open the file: %s\n", filename.c_str());
		exit(1);
	}

	fprintf(fp, "[Abstact]\n");
	fprintf(fp, "Alignment result of %s AND %s\n\n", m_chain_a->raw_name(), m_chain_b->raw_name());

	fprintf(fp, "[Result]\n");
	fprintf(fp, " %f %d\n\n", m_rmsd, m_align_num);

	fprintf(fp, "[Alignment]\n");
	for (i=0; i<m_length_a; i++) {
		fprintf(fp, " %3d", m_alignment[i]);
		if (i % 20 == 19) {
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n\n");
	
	fprintf(fp, "[Translation]\n");
	for (i=0; i<3; i++) {
		fprintf(fp, " %lf", m_translation[i]);
	}
	fprintf(fp, "\n\n");

	fprintf(fp, "[Rotation]\n");
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			fprintf(fp, " %lf", m_rotation[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fclose(fp);
}

int PairAlign::_getAlignNum(const vector<int> &alignment)
{
	int align_num, i;
	align_num = 0;
	for (i=0; i<m_length_a; i++) {
		if (alignment[i] >= 0) {
			align_num++;
		}
	}
	return align_num;
}

int PairAlign::_getBreakNum(const vector<int> &alignment)
{
	int break_num, i;
	bool in_break;
	break_num = 0;
	if (alignment[0] < 0) {
		in_break = true;
	}
	else {
		in_break = false;
	}
	for (i=1; i<m_length_a; i++) {
		if (in_break) {
			if (alignment[i] >= 0) in_break = false;
		}
		else if (alignment[i] < 0) {
			in_break = true;
			break_num++;
		}
		else if (alignment[i] != alignment[i-1]+1 && alignment[i] != alignment[i-1]-1) {
			break_num++;
		}
	}
	if (in_break && break_num > 0) break_num--;
	return break_num;
}

int PairAlign::_getPermuNum(const vector<int> &alignment)
{
	int permu_num, i, j;
	permu_num = 0;
	j = -1;
	for (i=0; i<m_length_a; i++) {
		if (alignment[i] >= 0) {
			if (alignment[i] < j) {
				permu_num++;
			}
			j = alignment[i];
		}
	}
	return permu_num;
}

double PairAlign::_getSequenceIdentity(const vector<int> &alignment)
{
	int i, s, n;
	double sequence_identity;
	s = 0;
	n = 0;
	for (i=0; i<m_length_a; i++) {
		if (alignment[i] >= 0) {
			if ((*m_chain_a)[i].isIdenticalResidue((*m_chain_b)[alignment[i]])) s++;
			n++;
		}
	}
	sequence_identity = (double)s / n;
	return sequence_identity;
}
