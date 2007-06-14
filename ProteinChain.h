
#ifndef __PROTEINCHAIN_H
#define __PROTEINCHAIN_H


#include <vector>
#include <cstdio>

#include "PDB.h"


class ProteinChain {
	PDB *m_pdb;
	char m_raw_name[256];
	char m_chain_id;								// Chain identifier
	int m_pocket_id;
	int m_range[2];
	bool m_is_backbone;
	char m_name[6];									// Name of chain
	char m_id_code[5];								// This identifier is unique within PDB
	char m_dep_date[10];							// Deposition date
	char m_classification[41];						// Classifies the molecule(s)
	vector<PDBAtom> m_atoms;						// Atoms

public:
	ProteinChain();
	~ProteinChain();

	PDBAtom &operator [](int i) { return m_atoms[i]; }

	const char *raw_name() const { return m_raw_name; }
	char chain_id() const { if (m_chain_id == ' ') return '_'; else return m_chain_id; }
	int pocket_id() const { return m_pocket_id; }
	int range(int i) const { return m_range[i]; }
	int length() const { return m_atoms.size(); }
	const char *name() const { return m_name; }
	const char *id_code() const { return m_id_code; }
	const char *dep_date() const { return m_dep_date; }
	const char *classification() const { return m_classification; }

	void setPDB(PDB *pdb);
	void setRawName(const char *rname);
	void setChainID(char cid);
	void setPocketID(int pid);
	void setRange(int start, int end);
	void setBackbone(bool enable);

	void getChain(PDB *pdb = NULL, char cid = 0);
	void getAllChains(PDB *pdb = NULL);

	void getPocketChain(PDB *pdb = NULL, int pid = 0);

	void clearData();

	void writeChainFile(const char *filename);
	void writeChainCode(FILE *fp);
	void writePDBModel(FILE *fp, int model = 0, bool fullchain = false, const double translation[3] = NULL, const double rotation[3][3] = NULL);

	double **getMatrix();
	double **getMatrix(const double translation[3], const double rotation[3][3]);
	double getRMSD(ProteinChain *chain, double translation[3], double rotation[3][3], int *alignment);

protected:
	void _writePDBModel(FILE *fp, int model = 0, const double translation[3] = NULL, const double rotation[3][3] = NULL);
	bool _filterChain(const PDBAtom &atom) const;
	bool _filterPocket(const PDBAtom &atom) const;
};


#endif // __PROTEINCHAIN_H
