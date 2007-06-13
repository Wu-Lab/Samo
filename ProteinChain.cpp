
#include <math.h>
#include <limits.h>

#include "ProteinChain.h"


////////////////////////////////
//
// class ProteinChain

ProteinChain::ProteinChain()
{
	int i;
	m_pdb = NULL;
	m_raw_name[0] = 0;
	m_chain_id = 0;
	m_pocket_id = 0;
	m_range[0] = m_range[1] = 0;
	m_is_backbone = true;
	m_length = 0;
	for (i=0; i<6; i++) m_name[i] = 0;
	for (i=0; i<5; i++) m_id_code[i] = 0;
	for (i=0; i<10; i++) m_dep_date[i] = 0;
	for (i=0; i<41; i++) m_classification[i] = 0;
	m_atoms = NULL;
}

ProteinChain::~ProteinChain()
{
	clearData();
}

void ProteinChain::setPDB(PDB *pdb)
{
	m_pdb = pdb;
}

void ProteinChain::setRawName(const char *rname)
{
	strcpy(m_raw_name, rname);
}

void ProteinChain::setChainID(char cid)
{
	if (cid >= 'a' && cid <= 'z') cid += 'A' - 'a';
	if (cid == '_') cid = ' ';
	m_chain_id = cid;
}

void ProteinChain::setPocketID(int pid)
{
	m_pocket_id = pid;
}

void ProteinChain::setRange(int start, int end)
{
	m_range[0] = start;
	m_range[1] = end;
}

void ProteinChain::setBackbone(bool enable)
{
	m_is_backbone = enable;
}

void ProteinChain::getChain(PDB *pdb, char cid)
{
	PDBAtom *atom;
	int i;

	if (pdb != NULL) setPDB(pdb);
	if (cid != 0) setChainID(cid);

	clearData();

	strcpy(m_id_code, m_pdb->id_code());
	strcpy(m_dep_date, m_pdb->dep_date());
	strcpy(m_classification, m_pdb->classification());
	if (m_chain_id == 0) m_chain_id = m_pdb->getChainID(1);
	if (m_chain_id < -1) {
		Logger::error("Can not find valid chain in PDB file %s!", m_pdb->filename());
		exit(1);
	}
	strcpy(m_name, m_id_code);
	if (m_chain_id != ' ') strncat(m_name, &m_chain_id, 1);
	else strcat(m_name, "_");

	atom = m_pdb->atom_list().getFirst();
	while (atom != NULL) {
		if (_filterChain(atom)) m_length++;
		atom = m_pdb->atom_list().getNext();
	}
	if (m_length == 0)
	{
		Logger::warning("Empty chain! PDB file: %s, Chain ID: %c!", m_pdb->filename(), m_chain_id);
	}
	m_atoms = new PDBAtom [m_length];

	i = 0;
	atom = m_pdb->atom_list().getFirst();
	while (atom != NULL) {
		if (_filterChain(atom)) m_atoms[i++] = *atom;
		atom = m_pdb->atom_list().getNext();
	}

	Logger::debug("Length of the protein chain %s:%c is %d", id_code(), chain_id(), length());
}

void ProteinChain::getAllChains(PDB *pdb)
{
	setChainID(-1);
	getChain(pdb);
}

void ProteinChain::getPocketChain(PDB *pdb, int pid)
{
	PDBAtom *atom, *last;
	int res_seq, i, n;
	double coord[3];

	if (pdb != NULL) setPDB(pdb);
	if (pid != 0) setPocketID(pid);

	if (m_pocket_id == 0) m_pocket_id = 1;
	if (m_pocket_id < -1) {
		Logger::error("Can not find valid pocket in POC file %s!", m_pdb->filename());
		exit(1);
	}
	if (m_chain_id == 0) m_chain_id = m_pdb->getChainID(1);
	if (m_chain_id < -1) {
		Logger::error("Can not find valid chain in POC file %s!", m_pdb->filename());
		exit(1);
	}
	strcpy(m_name, "POC");

	res_seq = INT_MAX;
	atom = m_pdb->atom_list().getFirst();
	while (atom != NULL) {
		if (_filterPocket(atom) && (atom->res_seq() != res_seq)) {
			m_length++;
			res_seq = atom->res_seq();
		}
		atom = m_pdb->atom_list().getNext();
	}
	if (m_length == 0)
	{
		Logger::warning("Empty pocket chain! Pocket file: %s, Pocket ID: %d!", m_pdb->filename(), m_pocket_id);
	}
	m_atoms = new PDBAtom [m_length];

	i = 0;
	res_seq = INT_MAX;
	last = NULL;
	atom = m_pdb->atom_list().getFirst();
	while (atom != NULL) {
		if (_filterPocket(atom)) {
			if (atom->res_seq() != res_seq) {
				if (last != NULL) {
					m_atoms[i] = *last;
					m_atoms[i][0] = coord[0] / n;
					m_atoms[i][1] = coord[1] / n;
					m_atoms[i][2] = coord[2] / n;
					i++;
				}
				res_seq = atom->res_seq();
				coord[0] = (*atom)[0];
				coord[1] = (*atom)[1];
				coord[2] = (*atom)[2];
				last = atom;
				n=1;
			}
			else {
				coord[0] += (*atom)[0];
				coord[1] += (*atom)[1];
				coord[2] += (*atom)[2];
				last = atom;
				n++;
			}
		}
		atom = m_pdb->atom_list().getNext();
	}
	if (last != NULL) {
		m_atoms[i] = *last;
		m_atoms[i][0] = coord[0] / n;
		m_atoms[i][1] = coord[1] / n;
		m_atoms[i][2] = coord[2] / n;
	}
}

void ProteinChain::clearData()
{
	m_length = 0;
	m_name[0] = 0;
	m_id_code[0] = 0;
	m_dep_date[0] = 0;
	m_classification[0] = 0;
	delete[] m_atoms;
	m_atoms = NULL;
}

void ProteinChain::writeChainFile(const char *filename)
{
	FILE *fp;
	int i;

	if ((fp = fopen(filename, "w")) == NULL) {
		Logger::error("Can not open the file: %s\n", filename);
		exit(1);
	}

	fprintf(fp, "%d\n", m_length);

	for (i=0; i<m_length; i++) {
		fprintf(fp, "%8.3f %8.3f %8.3f\n", m_atoms[i][0], m_atoms[i][1], m_atoms[i][2]);
	}

	fclose(fp);
}

void ProteinChain::writeChainCode(FILE *fp)
{
	int i;

	for (i=0; i<m_length; i++) {
		fprintf(fp, "%c", m_atoms[i].resCode());
	}
}

void ProteinChain::writePDBModel(FILE *fp, int model, bool fullchain, double translation[3], double rotation[3][3])
{
	ProteinChain *chain;
	if (m_pocket_id == 0 && fullchain) {
		chain = new ProteinChain;
		chain->setBackbone(false);
		chain->getChain(m_pdb, m_chain_id);
	}
	else {
		chain = this;
	}
	chain->_writePDBModel(fp, model, translation, rotation);
	if (m_pocket_id == 0 && fullchain) delete chain;
}

double **ProteinChain::getMatrix()
{
	double **matrix = Matrix<double>::alloc(m_length, 3);
	int i, j;
	for (i=0; i<m_length; i++) {
		for (j=0; j<3; j++) {
			matrix[i][j] = m_atoms[i][j];
		}
	}
	return matrix;
}

double **ProteinChain::getMatrix(double translation[3], double rotation[3][3])
{
	double **matrix = Matrix<double>::alloc(m_length, 3);
	int i, j, k;
	for (i=0; i<m_length; i++) {
		for (j=0; j<3; j++) {
			matrix[i][j] = translation[j];
		}
	}
	for (i=0; i<m_length; i++) {
		for (j=0; j<3; j++) {
			for (k=0; k<3; k++) {
				matrix[i][j] += rotation[j][k] * m_atoms[i][k];
			}
		}
	}
	return matrix;
}

double ProteinChain::getRMSD(ProteinChain *chain, double translation[3], double rotation[3][3], int *alignment)
{
	double **matrix = getMatrix(translation, rotation);
	double rmsd, dist;
	int i, j, n;
	rmsd = 0;
	n = 0;
	for (i=0; i<m_length; i++) {
		if (alignment[i] >= 0 && alignment[i] < chain->length()) {
			for (j=0; j<3; j++) {
				dist = matrix[i][j] - (*chain)[alignment[i]][j];
				rmsd += dist * dist;
			}
			n++;
		}
	}
	if (n > 0) {
		rmsd = sqrt(rmsd / n);
	}
	Matrix<double>::free(matrix);
	return rmsd;
}

void ProteinChain::_assign(ProteinChain &chain)
{
	int i;
	clearData();
	strcpy(m_name, chain.m_name);
	strcpy(m_id_code, chain.m_id_code);
	strcpy(m_dep_date, chain.m_dep_date);
	strcpy(m_classification, chain.m_classification);
	m_length = chain.m_length;
	m_atoms = new PDBAtom [m_length];
	for (i=0; i<m_length; i++) {
		m_atoms[i] = chain.m_atoms[i];
	}
	m_range[0] = chain.m_range[0];
	m_range[1] = chain.m_range[1];
	m_is_backbone = chain.m_is_backbone;
}

void ProteinChain::_writePDBModel(FILE *fp, int model, double translation[3], double rotation[3][3])
{
	double **matrix;
	int i;
	if (model > 0) fprintf(fp, "MODEL     %4d%66c\n", model, ' ');
	if (translation != NULL && rotation != NULL) {
		matrix = getMatrix(translation, rotation);
	}
	else {
		matrix = getMatrix();
	}
	for (i=0; i<m_length; i++) {
		fprintf(fp, "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%26c\n",
			m_atoms[i].serial(),
			m_atoms[i].name(),
			m_atoms[i].res_name(),
			m_atoms[i].chain_id(),
			m_atoms[i].res_seq(),
			matrix[i][0],
			matrix[i][1],
			matrix[i][2], ' ');
		if ((i == m_length-1) || (m_atoms[i].chain_id() != m_atoms[i+1].chain_id())) {
			fprintf(fp, "TER   %5d      %3s %c%4d%54c\n",
				m_atoms[i].serial()+1,
				m_atoms[i].res_name(),
				m_atoms[i].chain_id(),
				m_atoms[i].res_seq(), ' ');
		}
	}
	Matrix<double>::free(matrix);
	if (model > 0) fprintf(fp, "ENDMDL%74c\n", ' ');
}

inline bool ProteinChain::_filterChain(PDBAtom *atom)
{
	return ((atom->isCa() || !m_is_backbone)
			&& (atom->chain_id() == m_chain_id || m_chain_id == -1)
			&& (atom->res_seq() >= m_range[0] || m_range[0] <= 0)
			&& (atom->res_seq() <= m_range[1] || m_range[1] <= 0));
}

inline bool ProteinChain::_filterPocket(PDBAtom *atom)
{
	return ((atom->pocket_id() == m_pocket_id || m_pocket_id == -1)
			&& (atom->chain_id() == m_chain_id || m_chain_id == -1));
}