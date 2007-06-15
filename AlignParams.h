
#ifndef __ALIGNPARAMS_H
#define __ALIGNPARAMS_H


struct AlignParams {
	AlignParams();

	double lambda;
	int heuristic_start;
	bool branch_and_bound;
	bool sequential_order;
	bool annealing;
	double annealing_initial;
	double annealing_rate;
	int flagment_length;
};

inline AlignParams::AlignParams()
{
	lambda = 6.0;
	heuristic_start = 0;
	branch_and_bound = false;
	sequential_order = false;
	annealing = false;
	annealing_initial = 60.0;
	annealing_rate = 0.4;
	flagment_length = 5;
}


#endif // __ALIGNPARAMS_H