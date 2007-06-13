
#include <stdlib.h>

#include "Samo.h"


const char *Samo::m_version = "0.03.000";
const char *Samo::m_year = "2006";


void Samo::copyright()
{
	Logger::info("Samo - Protein Structure Alignment tool [Version %s]", m_version);
	Logger::info("Copyright (C) ZHANGroup@BIOINFOAMSS.ORG 2005-%s. All rights reserved.\n", m_year);
}

void Samo::usage()
{
	char buf[10000];
	Logger::info("Usage 1: Samo [options ...] [{code|file}[:id[:start[:end]]] ...]\n");
	Logger::info("  <code>  is the 4-character PDB assigned identifier");
	Logger::info("  <file>  is the PDB filename");
	Logger::info("  <id>    is the chain identifier in the PDB file");
	Logger::info("  <start> is the start residue sequence number");
	Logger::info("  <end>   is the end residue sequence number");
	Logger::info("");
	Logger::info("Usage 2: Samo -pocket [options ...] [file[:pid[:cid]] ...]\n");
	Logger::info("  <file>  is the Pocket filename");
	Logger::info("  <pid>   is the pocket identifier in the Pocket file");
	Logger::info("  <cid>   is the chain identifier in the Pocket file");
	Logger::info("");
	Logger::info(m_args.printDef(buf));
	Logger::info("Examples:");
	Logger::info("\tSamo 1dhf:a 8dfr");
	Logger::info("\tSamo pdb1dhf.ent:a pdb8dfr.ent");
	Logger::info("\tSamo 1dhf:a:20:120 8dfr:_:50:150");
}

Samo::Samo(int argc, char *argv[])
{
	defineOptions();
	parseOptions(argc, argv);

	m_chain_num = 0;
	m_pdbs = NULL;
	m_chains = NULL;
}

Samo::~Samo()
{
	clear();
}

void Samo::defineOptions()
{
	ArgOption *option;

	///////////////////////////////////////////////////////
	//
	// method options

	option = m_args.addOption("localization", "L");
	option->addHelpInfo("SAMO Localization version.");

	option = m_args.addOption("branch_and_bound", "b");
	option->addHelpInfo("Use branch and bound method instead of iterative method.");

	///////////////////////////////////////////////////////
	//
	// parameter options

	option = m_args.addOption("lambda", "l", true, "6.0");
	option->addHelpInfo("Parameter to balance two objective. Smaller value for smaller RMSD.");

	option = m_args.addOption("sequential_order", "so");
	option->addHelpInfo("Require alignment with sequential order.");

	///////////////////////////////////////////////////////
	//
	// heuristic options

	option = m_args.addOption("heuristic_start", "hs", true, "2");
	option->addHelpInfo("Set heuristic levels for finding initial solutions.");

	///////////////////////////////////////////////////////
	//
	// annealing options

	option = m_args.addOption("annealing_enable", "ae");
	option->addHelpInfo("Enable annealing technique.");

	option = m_args.addOption("annealing_initial", "ai", true, "60.0");
	option->addHelpInfo("Initial value for annealing.");

	option = m_args.addOption("annealing_rate", "ar", true, "0.4");
	option->addHelpInfo("Cooling coefficient for annealing.");

	///////////////////////////////////////////////////////
	//
	// special functions

	option = m_args.addOption("extract_only", "x");
	option->addHelpInfo("Extract chains data from the PDB files.");

	option = m_args.addOption("evaluate", "e", true);
	option->addHelpInfo("Evaluate a given alignment.");

	option = m_args.addOption("improve", "i", true);
	option->addHelpInfo("Improve a given alignment.");

	///////////////////////////////////////////////////////
	//
	// input options

	option = m_args.addOption("pocket", "p");
	option->addHelpInfo("Align two protein pockets instead of protein chains.");

	///////////////////////////////////////////////////////
	//
	// output options

	option = m_args.addOption("output_solfile", "os", true);
	option->addHelpInfo("Output alignment result to a solution file.");

	option = m_args.addOption("output_pdbfile", "op", true);
	option->addHelpInfo("Output alignment result to a PDB file.");

	///////////////////////////////////////////////////////
	//
	// help options

	option = m_args.addOption("debug_level", "d", true, "4", 5, "1", "2", "3", "4", "5");
	option->addHelpInfo("1: Error; 2: Warning; 3:Information; 4:Verbose; 5:Debug.");

	option = m_args.addOption("no_logo", "n");
	option->addHelpInfo("Suppress copyright message.");

	option = m_args.addOption("help", "h");
	option->addHelpInfo("Print this message.");
}

void Samo::parseOptions(int argc, char *argv[])
{
	char buf[10000];
	m_args.parseArguments(argc, argv);
	if (!m_args.isDefined("no_logo")) copyright();
	if (m_args.getNonOptionArgNum() < 1 || m_args.isDefined("help")) {
		usage();
		exit(0);
	}
	Logger::setLogLevel(m_args.getOption("debug_level")->getValueAsInt());
	Logger::debug(m_args.printVal(buf));
}

void Samo::run()
{
	int i;
	char filename[80];

	parseChainInfo();
	if (m_args.isDefined("extract_only")) {
		for (i=0; i<m_chain_num; i++) {
			strcpy(filename, m_chains[i].name());
			strcat(filename, ".d");
			m_chains[i].writeChainFile(filename);
			Logger::info("Length of the protein chain %s is %d", m_chains[i].name(), m_chains[i].length());
		}
	}
	else if (m_args.isDefined("evaluate")) {
		if (m_chain_num != 2) {
			Logger::error("2 protein chains are required for evaluation!");
			exit(1);
		}
		PairAlign palign(&m_chains[0], &m_chains[1]);
		palign.setLambda(m_args.getOption("lambda")->getValueAsDouble());
		palign.setSequentialOrder(m_args.isDefined("sequential_order"));
		palign.evaluate(m_args.getOption("evaluate")->getValue());
		palign.postProcess();
		if (m_args.isDefined("output_solfile")) {
			palign.writeSolutionFile(m_args.getOption("output_solfile")->getValue());
		}
		if (m_args.isDefined("output_pdbfile")) {
			palign.writePDBFile(m_args.getOption("output_pdbfile")->getValue());
		}
	}
	else if (m_args.isDefined("improve")) {
		if (m_chain_num != 2) {
			Logger::error("2 protein chains are required for improvement!");
			exit(1);
		}
		PairAlign palign(&m_chains[0], &m_chains[1]);
		palign.setLambda(m_args.getOption("lambda")->getValueAsDouble());
		palign.setSequentialOrder(m_args.isDefined("sequential_order"));
		palign.improve(m_args.getOption("improve")->getValue());
		palign.postProcess();
		if (m_args.isDefined("output_solfile")) {
			palign.writeSolutionFile(m_args.getOption("output_solfile")->getValue());
		}
		if (m_args.isDefined("output_pdbfile")) {
			palign.writePDBFile(m_args.getOption("output_pdbfile")->getValue());
		}
	}
	else if (m_chain_num <= 1) {
		Logger::error("At least 2 protein chains are required for alignment!");
		exit(1);
	}
	else if (m_chain_num == 2) {
		Logger::beginTimer(1, "Pairwise alignment");
		PairAlign palign(&m_chains[0], &m_chains[1]);
		palign.setLocalization(m_args.isDefined("localization"));
		palign.setBranchAndBound(m_args.isDefined("branch_and_bound"));
		palign.setLambda(m_args.getOption("lambda")->getValueAsDouble());
		palign.setSequentialOrder(m_args.isDefined("sequential_order"));
		palign.setHeuristicStart(m_args.getOption("heuristic_start")->getValueAsInt());
		palign.setAnnealing(m_args.isDefined("annealing_enable"));
		palign.setAnnealingInitial(m_args.getOption("annealing_initial")->getValueAsDouble());
		palign.setAnnealingRate(m_args.getOption("annealing_rate")->getValueAsDouble());
		palign.align();
		palign.postProcess();
		Logger::endTimer(1);
		if (m_args.isDefined("output_solfile")) {
			palign.writeSolutionFile(m_args.getOption("output_solfile")->getValue());
		}
		if (m_args.isDefined("output_pdbfile")) {
			palign.writePDBFile(m_args.getOption("output_pdbfile")->getValue());
		}
	}
	else {
		Logger::beginTimer(1, "Multiple alignment");
		MultiAlign malign(m_chain_num);
//		malign.setBranchAndBound(m_args.isDefined("branch_and_bound"));
		malign.setLambda(m_args.getOption("lambda")->getValueAsDouble());
		malign.setSequentialOrder(m_args.isDefined("sequential_order"));
		malign.setHeuristicStart(m_args.getOption("heuristic_start")->getValueAsInt());
		malign.setAnnealing(m_args.isDefined("annealing_enable"));
		malign.setAnnealingInitial(m_args.getOption("annealing_initial")->getValueAsDouble());
		malign.setAnnealingRate(m_args.getOption("annealing_rate")->getValueAsDouble());
		for (i=0; i<m_chain_num; i++) {
			malign.setChain(i, &m_chains[i]);
		}
		malign.align();
		Logger::endTimer(1);
		if (m_args.isDefined("output_pdbfile")) {
			malign.writePDBFile(m_args.getOption("output_pdbfile")->getValue());
		}
	}
}

void Samo::clear()
{
	delete[] m_pdbs;
	m_pdbs = NULL;
	delete[] m_chains;
	m_chains = NULL;
}

void Samo::parseChainInfo()
{
	int i, start, end;
	char filename[1024], buffer[1024], argument[1024], *token;

	m_chain_num = m_args.getNonOptionArgNum();
	m_pdbs = new PDB [m_chain_num];
	m_chains = new ProteinChain [m_chain_num];

	if (!m_args.isDefined("pocket")) {
		for (i=0; i<m_chain_num; i++) {
			strcpy(argument, m_args.getNonOptionArgument(i));
			m_chains[i].setRawName(argument);
			token = strtok(argument, ":");
			strcpy(filename, token);
			if (strchr(filename, '.') == NULL) {
				strcat(filename, ".ent");
				if (strncmp(filename, "pdb", 3) != 0) {
					strcpy(buffer, "pdb");
					strcat(buffer, filename);
					strcpy(filename, buffer);
				}
			}
			m_pdbs[i].readFile(filename);
			m_chains[i].setPDB(&m_pdbs[i]);
			token = strtok(NULL, ":");
			if (token != NULL) {
				parseChainID(i, token);
				token = strtok(NULL, ":");
				if (token != NULL) {
					start = atoi(token);
					token = strtok(NULL, ":");
					if (token != NULL) {
						end = atoi(token);
					}
					else {
						end = 0;
					}
					m_chains[i].setRange(start, end);
				}
			}
			m_chains[i].getChain();
			if (m_chains[i].length() == 0) exit(1);
		}
	}
	else {
		for (i=0; i<m_chain_num; i++) {
			strcpy(argument, m_args.getNonOptionArgument(i));
			m_chains[i].setRawName(argument);
			token = strtok(argument, ":");
			strcpy(filename, token);
			if (strchr(filename, '.') == NULL) {
				strcat(filename, ".poc");
//				if (strncmp(filename, "pdb", 3) != 0) {
//					strcpy(buffer, "pdb");
//					strcat(buffer, filename);
//					strcpy(filename, buffer);
//				}
			}
			m_pdbs[i].readPocket(filename);
			m_chains[i].setPDB(&m_pdbs[i]);
			token = strtok(NULL, ":");
			if (token != NULL) {
				parsePocketID(i, token);
				token = strtok(NULL, ":");
				if (token != NULL) {
					parseChainID(i, token);
				}
			}
			m_chains[i].getPocketChain();
			if (m_chains[i].length() == 0) exit(1);
		}
	}
}

void Samo::parsePocketID(int i, char *token)
{
	if (token[0] == '#') {
		m_chains[i].setPocketID(m_pdbs[i].getPocketID(atoi(token+1)));
	}
	else {
		m_chains[i].setPocketID(atoi(token));
	}
}

void Samo::parseChainID(int i, char *token)
{
	if (token[0] == '#') {
		m_chains[i].setChainID(m_pdbs[i].getChainID(atoi(token+1)));
	}
	else {
		m_chains[i].setChainID(token[0]);
	}
}