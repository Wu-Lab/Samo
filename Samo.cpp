
#include <iostream>
#include <fstream>

#include "Samo.h"
#include "Options.h"


const char *Samo::m_version = "0.04.000";
const char *Samo::m_year = "2007";


Samo::Samo(int argc, char *argv[])
{
	po::options_description generics("Generic options");
	generics.add_options()
		("help,h", "Show help message")
		("config,c", po::value<string>()->default_value("Samo.cfg"), "Configuration file")
		;

	po::options_description configs("Configurations");
	configs.add_options()
		("nologo", "Suppress logo and copyright information")
		("debug,d", po::value<int>()->default_value(4), "Set debug level - 1: Error; 2: Warning; 3:Information; 4:Verbose; 5:Debug")
		("pocket", "Align two protein pockets instead of protein chains")
		("output-solution", po::value<string>(), "Output alignment result to a solution file")
		("output-pdb", po::value<string>(), "Output alignment result to a PDB file")
		;

	po::options_description parameters("Parameters");
	parameters.add_options()
		("branch-and-bound,b", "Use branch and bound method instead of iterative method")
		("lambda,l", po::value<double>()->default_value(6.0), "Parameter to balance two objective, smaller value for smaller RMSD")
		("sequential-order", "Require alignment with sequential order")
		("heuristic-start", po::value<int>()->default_value(2), "Set heuristic levels for finding initial solutions")
		("annealing", "Enable annealing technique")
		("annealing-initial", po::value<double>()->default_value(60.0), "Initial value for annealing")
		("annealing-rate", po::value<double>()->default_value(0.4), "Cooling coefficient for annealing")
		;

	po::options_description utilities("Utility options");
	utilities.add_options()
		("extract", "Extract chains data from the PDB files")
		("evaluate", po::value<string>(), "Evaluate a given alignment")
		("improve", po::value<string>(), "Improve a given alignment")
		;

	po::options_description hidden;
	hidden.add_options()
		("filename", po::value<vector<string> >(&m_filenames), "input files")
		;

	m_options.add(generics).add(configs).add(parameters).add(utilities).add(hidden);
	m_visible_options.add(generics).add(configs).add(parameters).add(utilities);

	po::positional_options_description p;
	p.add("filename", -1);

	po::options_description cmdline_options, file_options;
	cmdline_options.add(generics).add(configs).add(parameters).add(utilities).add(hidden);
	file_options.add(configs).add(parameters).add(utilities);

	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), m_args);

	copyright();
	ifstream ifs(m_args["config"].as<string>().c_str());
	po::store(po::parse_config_file(ifs, file_options), m_args);

	po::notify(m_args);

//	conflicting_options(m_args, "convert", "compare");

	parseOptions();
	parseFileNames();
}

Samo::~Samo()
{
}

void Samo::copyright()
{
	Logger::info("Samo - Protein Structure Alignment tool [Version %s]", m_version);
	Logger::info("Copyright (C) ZHANGroup@APORC.ORG 2005-%s. All rights reserved.\n", m_year);
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
//	Logger::info(m_args.printDef(buf));
	Logger::info("Examples:");
	Logger::info("\tSamo 1dhf:a 8dfr");
	Logger::info("\tSamo pdb1dhf.ent:a pdb8dfr.ent");
	Logger::info("\tSamo 1dhf:a:20:120 8dfr:_:50:150");
}

void Samo::printOptions()
{
	cout << "Used options (user specified):" << endl;
	print_options(m_args, cout, DisplayOption::specified);
	cout << "Used options (default):" << endl;
	print_options(m_args, cout, DisplayOption::defaulted);
	cout << endl;
}

void Samo::parseOptions()
{
	Logger::setLogLevel(m_args["debug"].as<int>());
	if (!m_args.count("nologo")) {
		copyright();
	}
	if (m_args.count("help") || !m_filenames.size()) {
		usage();
		exit(0);
	}
	if (Logger::isDebug()) {
		printOptions();
	}
}

void Samo::run()
{
	int i;
	char filename[80];

	if (m_args.count("extract")) {
		for (i=0; i<m_chain_num; i++) {
			strcpy(filename, m_chains[i].name());
			strcat(filename, ".d");
			m_chains[i].writeChainFile(filename);
			Logger::info("Length of the protein chain %s is %d", m_chains[i].name(), m_chains[i].length());
		}
	}
	else if (m_args.count("evaluate")) {
		if (m_chain_num != 2) {
			Logger::error("2 protein chains are required for evaluation!");
			exit(1);
		}
		PairAlign palign(&m_chains[0], &m_chains[1]);
		palign.setLambda(m_args["lambda"].as<double>());
		palign.setSequentialOrder(m_args.count("sequential-order"));
		palign.evaluate(m_args["evaluate"].as<string>());
		palign.postProcess();
		output(palign);
	}
	else if (m_args.count("improve")) {
		if (m_chain_num != 2) {
			Logger::error("2 protein chains are required for improvement!");
			exit(1);
		}
		PairAlign palign(&m_chains[0], &m_chains[1]);
		palign.setLambda(m_args["lambda"].as<double>());
		palign.setSequentialOrder(m_args.count("sequential-order"));
		palign.improve(m_args["improve"].as<string>());
		palign.postProcess();
		output(palign);
	}
	else if (m_chain_num <= 1) {
		Logger::error("At least 2 protein chains are required for alignment!");
		exit(1);
	}
	else if (m_chain_num == 2) {
		Logger::beginTimer(1, "Pairwise alignment");
		PairAlign palign(&m_chains[0], &m_chains[1]);
		palign.setBranchAndBound(m_args.count("branch-and-bound"));
		palign.setLambda(m_args["lambda"].as<double>());
		palign.setSequentialOrder(m_args.count("sequential-order"));
		palign.setHeuristicStart(m_args["heuristic-start"].as<int>());
		palign.setAnnealing(m_args.count("annealing"));
		palign.setAnnealingInitial(m_args["annealing-initial"].as<double>());
		palign.setAnnealingRate(m_args["annealing-rate"].as<double>());
		palign.align();
		palign.postProcess();
		Logger::endTimer(1);
		output(palign);
	}
	else {
		Logger::beginTimer(1, "Multiple alignment");
		MultiAlign malign(m_chain_num);
//		malign.setBranchAndBound(m_args.count("branch-and-bound"));
		malign.setLambda(m_args["lambda"].as<double>());
		malign.setSequentialOrder(m_args.count("sequential-order"));
		malign.setHeuristicStart(m_args["heuristic-start"].as<int>());
		malign.setAnnealing(m_args.count("annealing"));
		malign.setAnnealingInitial(m_args["annealing-initial"].as<double>());
		malign.setAnnealingRate(m_args["annealing-rate"].as<double>());
		for (i=0; i<m_chain_num; i++) {
			malign.setChain(i, &m_chains[i]);
		}
		malign.align();
		Logger::endTimer(1);
		output(malign);
	}
}

void Samo::parseFileNames()
{
	int i, start, end;
	char filename[1024], buffer[1024], argument[1024], *token;

	m_chain_num = m_filenames.size();
	m_pdbs.resize(m_chain_num);
	m_chains.resize(m_chain_num);

	if (!m_args.count("pocket")) {
		for (i=0; i<m_chain_num; i++) {
			strcpy(argument, m_filenames[i].c_str());
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
			strcpy(argument, m_filenames[i].c_str());
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