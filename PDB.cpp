
#include <set>
#include <string>
#include <cstdio>

#include "Utils.h"
#include "PDB.h"


////////////////////////////////
//
// class PDBAtom

PDBAtom::PDBAtom()
{
	int i;

	m_serial = 0;
	for (i=0; i<5; i++) m_name[i] = 0;
	for (i=0; i<4; i++) m_res_name[i] = 0;
	m_chain_id = ' ';
	m_pocket_id = 0;
	m_res_seq = 0;
	m_coord[0] = 0;
	m_coord[1] = 0;
	m_coord[2] = 0;
}

bool PDBAtom::isCa() const
{
	return (strcmp(m_name, " CA ") == 0);
}

bool PDBAtom::isIdenticalAtom(const PDBAtom &atom) const
{
	return (strcmp(m_name, atom.m_name) == 0);
}

bool PDBAtom::isIdenticalResidue(const PDBAtom &atom) const
{
	return (strcmp(m_res_name, atom.m_res_name) == 0);
}

const char PDBAtom::_getResidueCode(const char *name)
{
	char code;
	if (strncmp(name, "ALA", 3) == 0) code = 'A';
	else if (strncmp(name, "ASX", 3) == 0) code = 'B';
	else if (strncmp(name, "CYS", 3) == 0) code = 'C';
	else if (strncmp(name, "ASP", 3) == 0) code = 'D';
	else if (strncmp(name, "GLU", 3) == 0) code = 'E';
	else if (strncmp(name, "PHE", 3) == 0) code = 'F';
	else if (strncmp(name, "GLY", 3) == 0) code = 'G';
	else if (strncmp(name, "HIS", 3) == 0) code = 'H';
	else if (strncmp(name, "ILE", 3) == 0) code = 'I';
	else if (strncmp(name, "LYS", 3) == 0) code = 'K';
	else if (strncmp(name, "LEU", 3) == 0) code = 'L';
	else if (strncmp(name, "MET", 3) == 0) code = 'M';
	else if (strncmp(name, "ASN", 3) == 0) code = 'N';
	else if (strncmp(name, "PRO", 3) == 0) code = 'P';
	else if (strncmp(name, "GLN", 3) == 0) code = 'Q';
	else if (strncmp(name, "ARG", 3) == 0) code = 'R';
	else if (strncmp(name, "SER", 3) == 0) code = 'S';
	else if (strncmp(name, "THR", 3) == 0) code = 'T';
	else if (strncmp(name, "VAL", 3) == 0) code = 'V';
	else if (strncmp(name, "TRP", 3) == 0) code = 'W';
	else if (strncmp(name, "TYR", 3) == 0) code = 'Y';
	else if (strncmp(name, "GLZ", 3) == 0) code = 'Z';
	else code = '-';
	return code;
}


////////////////////////////////
//
// class PDB

PDB::PDB(const char *filename)
{
	if (filename != NULL) setFilename(filename);
	else setFilename("unkown");
	clearData();
}

void PDB::setFilename(const char *filename)
{
	strcpy(m_filename, filename);
}

int PDB::getPocketID(int index)
{
	set<int> pockets;
	int i, j, pid;

	pid = -2;
	for (i=0; i<m_atoms.size(); ++i) {
		j = m_atoms[i].pocket_id();
		if (!pockets.count(j)) {
			pockets.insert(j);
			if (pockets.size() >= index) {
				pid = j;
				break;
			}
		}
	}
	return pid;
}

char PDB::getChainID(int index)
{
	set<char> chains;
	int i, j;
	char cid;

	cid = -2;
	for (i=0; i<m_atoms.size(); ++i) {
		j = m_atoms[i].chain_id();
		if (!chains.count(j)) {
			chains.insert(j);
			if (chains.size() >= index) {
				cid = j;
				break;
			}
		}
	}
	return cid;
}

void PDB::readFile(const char *filename)
{	
	FILE *fp;
	char buffer[100];
	int i;

	if (filename != NULL) setFilename(filename);

	if ((fp = fopen(m_filename, "r")) == NULL) {
		Logger::error("Can not open the file: %s\n", m_filename);
		exit(1);
	}

	clearData();
	Logger::debug("Read PDB file: %s", m_filename);

	i = 0;
	while (!feof(fp)) {
		fgets(buffer, 100, fp);

/////////////////////////////////////////////////////////////////////////////////////
// The HEADER record uniquely identifies a PDB entry through the idCode field. This
// record also provides a classification for the entry. Finally, it contains the
// date the coordinates were deposited at the PDB. 
//
// COLUMNS        DATA TYPE       FIELD           DEFINITION
// ----------------------------------------------------------------------------------
//  1 -  6        Record name     "HEADER"
// 11 - 50        String(40)      classification  Classifies the molecule(s)
// 51 - 59        Date            depDate         Deposition date.  This is the date
//                                                the coordinates were received by
//                                                the PDB
// 63 - 66        IDcode          idCode          This identifier is unique within PDB
/////////////////////////////////////////////////////////////////////////////////////

		if (strncmp(buffer, "HEADER", 6) == 0) {
			sscanf(buffer, "%*10c%40c%9c%*3c%4c\n", m_classification, m_dep_date, m_id_code);
		}

/////////////////////////////////////////////////////////////////////////////////////
// The COMPND record describes the macromolecular contents of an entry. Each
// macromolecule found in the entry is described by a set of token: value pairs,
// and is referred to as a COMPND record component. Since the concept of a molecule
// is difficult to specify exactly, PDB staff may exercise editorial judgment in
// consultation with depositors in assigning these names.
// For each macromolecular component, the molecule name, synonyms, number assigned
// by the Enzyme Commission (EC), and other relevant details are specified. 
//
// COLUMNS        DATA TYPE         FIELD          DEFINITION
// ----------------------------------------------------------------------------------
//  1 -  6        Record name       "COMPND"
//  9 - 10        Continuation      continuation   Allows concatenation of multiple records.
// 11 - 70        Specification     compound       Description of the molecular
//                list                             components.
/////////////////////////////////////////////////////////////////////////////////////

		else if (strncmp(buffer, "COMPND", 6) == 0) {
			char *buf = new char [61];
			buf[60] = 0;
			sscanf(buffer, "%*10c%60c\n", buf);
			m_components.push_back(buf);
		}

/////////////////////////////////////////////////////////////////////////////////////
// The MASTER record is a control record for bookkeeping. It lists the number of
// lines in the coordinate entry or file for selected record types. 
//
// COLUMNS       DATA TYPE      FIELD         DEFINITION
// ----------------------------------------------------------------------------------
//  1 -  6       Record name    "MASTER"
// 11 - 15       Integer        numRemark     Number of REMARK records
// 16 - 20       Integer        "0"
// 21 - 25       Integer        numHet        Number of HET records
// 26 - 30       Integer        numHelix      Number of HELIX records
// 31 - 35       Integer        numSheet      Number of SHEET records
// 36 - 40       Integer        numTurn       Number of TURN records
// 41 - 45       Integer        numSite       Number of SITE records
// 46 - 50       Integer        numXform      Number of coordinate transformation
//                                            records (ORIGX+SCALE+MTRIX)
// 51 - 55       Integer        numCoord      Number of atomic coordinate records
//                                            (ATOM+HETATM)
// 56 - 60       Integer        numTer        Number of TER records
// 61 - 65       Integer        numConect     Number of CONECT records
// 66 - 70       Integer        numSeq        Number of SEQRES records
/////////////////////////////////////////////////////////////////////////////////////

		else if (strncmp(buffer, "MASTER", 6) == 0)	{
			sscanf(buffer, "%*50c%5d\n", &m_num_coord);
		}

/////////////////////////////////////////////////////////////////////////////////////
// The ATOM records present the atomic coordinates for standard residues. They also
// present the occupancy and temperature factor for each atom. Heterogen coordinates
// use the HETATM record type. The element symbol is always present on each ATOM
// record; segment identifier and charge are optional. 
//
// COLUMNS        DATA TYPE       FIELD         DEFINITION
// ---------------------------------------------------------------------------------
//  1 -  6        Record name     "ATOM  "
//  7 - 11        Integer         serial        Atom serial number.
// 13 - 16        Atom            name          Atom name.
// 17             Character       altLoc        Alternate location indicator.
// 18 - 20        Residue name    resName       Residue name.
// 22             Character       chainID       Chain identifier.
// 23 - 26        Integer         resSeq        Residue sequence number.
// 27             AChar           iCode         Code for insertion of residues.
// 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
//                                              Angstroms.
// 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
//                                              Angstroms.
// 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
//                                              Angstroms.
// 55 - 60        Real(6.2)       occupancy     Occupancy.
// 61 - 66        Real(6.2)       tempFactor    Temperature factor.
// 73 - 76        LString(4)      segID         Segment identifier, left-justified.
// 77 - 78        LString(2)      element       Element symbol, right-justified.
// 79 - 80        LString(2)      charge        Charge on the atom.
/////////////////////////////////////////////////////////////////////////////////////

		else if (strncmp(buffer, "ATOM  ", 6) == 0) {
			PDBAtom atom;
			sscanf(buffer, "%*6c%5d%*c%4c%*c%3c%*c%c%4d%*c%*3c%8lf%8lf%8lf\n",
				&atom.m_serial,
				atom.m_name,
				atom.m_res_name,
				&atom.m_chain_id,
				&atom.m_res_seq,
				&atom.m_coord[0],
				&atom.m_coord[1],
				&atom.m_coord[2]);
//			Logger::debug("%5d, %4s, %3s, %c, %4d, %8.3f, %8.3f, %8.3f\n",
//				atom->m_serial,
//				atom->m_name,
//				atom->m_res_name,
//				atom->m_chain_id,
//				atom->m_res_seq,
//				atom->m_coord[0],
//				atom->m_coord[1],
//				atom->m_coord[2]);
			m_atoms.push_back(atom);
		}	

/////////////////////////////////////////////////////////////////////////////////////
// The ENDMDL records are paired with MODEL records to group individual structures
// found in a coordinate entry. 

// COLUMNS         DATA TYPE        FIELD           DEFINITION
// ------------------------------------------------------------------
//  1 -  6         Record name      "ENDMDL"
/////////////////////////////////////////////////////////////////////////////////////

		else if (strncmp(buffer, "ENDMDL", 6) == 0) {
			break;
		}

/////////////////////////////////////////////////////////////////////////////////////
// The END record marks the end of the PDB file. 
//
// COLUMNS       DATA TYPE      FIELD     DEFINITION
// -------------------------------------------------------
//  1 -  6       Record name    "END   "
/////////////////////////////////////////////////////////////////////////////////////

		else if (strncmp(buffer, "END   ", 6) == 0) {
			break;
		}
	}
	fclose(fp);

	Logger::debug("\tDatabase Code: %s", m_id_code);
	Logger::debug("\tClassification: %s", m_classification);
	Logger::debug("\tNumber of Coordinate Records: %d", m_num_coord);
}

void PDB::readPocket(const char *filename)
{	
	FILE *fp;
	char buffer[100];
	int i;

	if (filename != NULL) setFilename(filename);

	if ((fp = fopen(m_filename, "r")) == NULL) {
		Logger::error("Can not open the file: %s\n", m_filename);
		exit(1);
	}

	clearData();
	Logger::debug("Read Pocket file: %s", m_filename);

	i = 0;
	while (!feof(fp)) {
		fgets(buffer, 100, fp);

/////////////////////////////////////////////////////////////////////////////////////
// The ATOM records present the atomic coordinates for standard residues. They also
// present the occupancy and temperature factor for each atom. Heterogen coordinates
// use the HETATM record type. The element symbol is always present on each ATOM
// record; segment identifier and charge are optional. 
//
// COLUMNS        DATA TYPE       FIELD         DEFINITION
// ---------------------------------------------------------------------------------
//  1 -  6        Record name     "ATOM  "
//  7 - 11        Integer         serial        Atom serial number.
// 13 - 16        Atom            name          Atom name.
// 17             Character       altLoc        Alternate location indicator.
// 18 - 20        Residue name    resName       Residue name.
// 22             Character       chainID       Chain identifier.
// 23 - 26        Integer         resSeq        Residue sequence number.
// 27             AChar           iCode         Code for insertion of residues.
// 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
//                                              Angstroms.
// 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
//                                              Angstroms.
// 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
//                                              Angstroms.
// 55 - 60        Real(6.2)       occupancy     Occupancy.
// 61 - 66        Real(6.2)       tempFactor    Temperature factor.
//
// 67 - 70        Integer         pocketID      Pocket identifier.
//
// 73 - 76        LString(4)      segID         Segment identifier, left-justified.
// 77 - 78        LString(2)      element       Element symbol, right-justified.
// 79 - 80        LString(2)      charge        Charge on the atom.
/////////////////////////////////////////////////////////////////////////////////////

		if (strncmp(buffer, "ATOM  ", 6) == 0) {
			PDBAtom atom;
			sscanf(buffer, "%*6c%5d%*c%4c%*c%3c%*c%c%4d%*c%*3c%8lf%8lf%8lf%*6lf%*6lf%4d\n",
				&atom.m_serial,
				atom.m_name,
				atom.m_res_name,
				&atom.m_chain_id,
				&atom.m_res_seq,
				&atom.m_coord[0],
				&atom.m_coord[1],
				&atom.m_coord[2],
				&atom.m_pocket_id);
// 			Logger::debug("%5d, %4s, %3s, %c, %4d, %8.3f, %8.3f, %8.3f, %4d\n",
// 				atom.m_serial,
// 				atom.m_name,
// 				atom.m_res_name,
// 				atom.m_chain_id,
// 				atom.m_res_seq,
// 				atom.m_coord[0],
// 				atom.m_coord[1],
// 				atom.m_coord[2],
// 				atom.m_pocket_id);
			m_atoms.push_back(atom);
		}
	}
	fclose(fp);
}

void PDB::clearData()
{
	int i;

	for (i=0; i<5; i++) m_id_code[i] = 0;
	for (i=0; i<10; i++) m_dep_date[i] = 0;
	for (i=0; i<41; i++) m_classification[i] = 0;
	m_components.clear();
	m_num_coord = 0;
	m_atoms.clear();
}
