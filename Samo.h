
#ifndef __SAMO_H
#define __SAMO_H


#include <boost/program_options.hpp>

#include "Utils.h"
#include "PDB.h"
#include "ProteinChain.h"
#include "PairAlign.h"
#include "MultiAlign.h"


namespace po = ::boost::program_options;


class Samo {
	po::options_description m_options;
	po::options_description m_visible_options;
	po::variables_map m_args;

	int m_chain_num;
	vector<string> m_filenames;
	vector<PDB> m_pdbs;
	vector<ProteinChain> m_chains;

	static const char *m_version;
	static const char *m_year;

public:
	Samo(int argc = 0, char *argv[] = NULL);
	~Samo();

	void copyright();
	void usage();
	void printOptions();
	void parseOptions();

	void run();

	template<class A>
	void output(const A &align);

	void parseFileNames();
	void parseChainID(int i, char *token);
	void parsePocketID(int i, char *token);
};


template<class A>
void Samo::output(const A &align)
{
	if (m_args.count("output-solution")) {
		align.writeSolutionFile(m_args["output-solution"].as<string>());
	}
	if (m_args.count("output-pdb")) {
		align.writePDBFile(m_args["output-pdb"].as<string>());
	}
}


#endif // __SAMO_H
