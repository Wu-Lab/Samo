
#ifndef __PDB_H
#define __PDB_H


/////////////////////////////////////////////////////////////////////////////////////
// Protein Data Bank
// Atomic Coordinate Entry Format Description
// Version 2.2, December 20, 1996
// http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
/////////////////////////////////////////////////////////////////////////////////////


#include <vector>
#include <string>


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

	double &operator [](int i) { return m_coord[i]; }
	const double &operator [](int i) const { return m_coord[i]; }

	int serial() const { return m_serial; }
	const char *name() const { return m_name; }
	const char *res_name() const { return m_res_name; }
	const char chain_id() const { return m_chain_id; }
	const int pocket_id() const { return m_pocket_id; }
	int res_seq() const { return m_res_seq; }

	bool isCa() const;
	bool isIdenticalAtom(const PDBAtom &atom) const;
	bool isIdenticalResidue(const PDBAtom &atom) const;

	const char resCode() const { return _getResidueCode(m_res_name); }

protected:
	static const char _getResidueCode(const char *name);

	friend class PDB;
};


class PDB {
	string m_filename;
	char m_id_code[5];								// This identifier is unique within PDB
	char m_dep_date[10];							// Deposition date
	char m_classification[41];						// Classifies the molecule(s)
	vector<string> m_components;					// Description of the molecular components
	int m_num_coord;								// Number of atomic coordinate records (ATOM+HETATM)
	vector<PDBAtom> m_atoms;						// List of atoms

public:
	PDB(const char *filename = NULL);

	const char *filename() const { return m_filename.c_str(); }
	const char *id_code() const { return m_id_code; }
	const char *dep_date() const { return m_dep_date; }
	const char *classification() const { return m_classification; }
	vector<PDBAtom> &atoms() { return m_atoms; }

	void setFilename(const string &filename) { m_filename = filename; }

	int getPocketID(int index);
	char getChainID(int index);

	void readFile(const string &filename = string());		// read data from PDB file
	void writeFile(const string &filename);					// write data to PDB file

	void readPocket(const string &filename = string());

	void clearData();
};


#endif // __PDB_H
