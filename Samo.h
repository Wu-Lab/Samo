
#ifndef __SAMO_H
#define __SAMO_H


#include "Utils.h"
#include "PDB.h"
#include "ProteinChain.h"
#include "PairAlign.h"
#include "MultiAlign.h"


class Samo {
	ArgParser m_args;
	int m_chain_num;
	PDB *m_pdbs;
	ProteinChain *m_chains;

	static const char *m_version;
	static const char *m_year;

public:
	Samo(int argc = 0, char *argv[] = NULL);
	~Samo();

	void copyright();
	void usage();
	void defineOptions();
	void parseOptions(int argc = 0, char *argv[] = NULL);

	void run();

	void clear();

	void parseChainInfo();
	void parseChainID(int i, char *token);
	void parsePocketID(int i, char *token);
};


#endif // __SAMO_H
