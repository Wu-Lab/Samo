
#ifndef __PDB_H
#define __PDB_H


/////////////////////////////////////////////////////////////////////////////////////
// Protein Data Bank
// Atomic Coordinate Entry Format Description
// Version 2.2, December 20, 1996
// http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
/////////////////////////////////////////////////////////////////////////////////////


#include "Utils.h"


class PDBAtom {
	int m_serial;			// Atom serial number
	char m_name[5];			// Atom name
	char m_res_name[4];		// Residue name
	char m_chain_id;		// Chain identifier
	int m_pocket_id;		// Pocket identifier
	int m_res_seq;			// Residue sequence number
	double m_coord[3];		// Coordinators of atom

public:
	PDBAtom();

	int serial() { return m_serial; }
	const char *name() { return m_name; }
	const char *res_name() { return m_res_name; }
	const char chain_id() { return m_chain_id; }
	const int pocket_id() { return m_pocket_id; }
	int res_seq() { return m_res_seq; }

	bool isCa();
	bool isIdenticalAtom(PDBAtom &atom);
	bool isIdenticalResidue(PDBAtom &atom);

	const char resCode() { return _getResidueCode(m_res_name); }

	double &operator [](int i) { return m_coord[i]; }
	void operator =(PDBAtom &atom) { _assign(atom); }

protected:
	void _assign(PDBAtom &atom);
	const char _getResidueCode(const char *name);

	friend class PDB;
};


class PDB {
	char m_filename[256];
	char m_id_code[5];								// This identifier is unique within PDB
	char m_dep_date[10];							// Deposition date
	char m_classification[41];						// Classifies the molecule(s)
	StringList m_compound_list;						// Description of the molecular components
	int m_num_coord;								// Number of atomic coordinate records (ATOM+HETATM)
	List<PDBAtom> m_atom_list;						// List of atoms

public:
	PDB(const char *filename = NULL);

	const char *filename() { return m_filename; }
	const char *id_code() { return m_id_code; }
	const char *dep_date() { return m_dep_date; }
	const char *classification() { return m_classification; }
	List<PDBAtom> &atom_list() { return m_atom_list; }

	void setFilename(const char *filename);

	int getPocketID(int index);
	char getChainID(int index);

	void readFile(const char *filename = NULL);				// read data from PDB file
	void writeFile(const char *filename);					// write data to PDB file

	void readPocket(const char *filename = NULL);

	void clearData();
};


#endif // __PDB_H
