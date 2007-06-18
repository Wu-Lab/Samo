
#include <iostream>
#include <fstream>

#include "Samo.h"
#include "Options.h"

#include "MemLeak.h"


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
		("pocket,p", "Align two protein pockets instead of protein chains")
		("output-solution", po::value<string>(), "Output alignment result to a solution file")
		("output-pdb", po::value<string>(), "Output alignment result to a PDB file")
		;

	po::options_description parameters("Parameters");
	parameters.add_options()
		("branch-and-bound,b", po::value<bool>(&m_params.branch_and_bound), "Use branch and bound method instead of iterative method")
		("lambda,l", po::value<double>(&m_params.lambda)->default_value(6.0), "Parameter to balance two objective, smaller value for smaller RMSD")
		("sequential-order", po::value<bool>(&m_params.sequential_order), "Require alignment with sequential order")
		("heuristic-start", po::value<int>(&m_params.heuristic_start)->default_value(2), "Set heuristic levels for finding initial solutions")
		("annealing", po::value<bool>(&m_params.annealing), "Enable annealing technique")
		("annealing-initial", po::value<double>(&m_params.annealing_initial)->default_value(60.0), "Initial value for annealing")
		("annealing-rate", po::value<double>(&m_params.annealing_rate)->default_value(0.4), "Cooling coefficient for annealing")
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
	if (Logger::log_level() >= Logger::log_level_info()) {
		cout << m_visible_options << "\n\n";
	}
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
		palign.setParams(m_params);
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
		palign.setParams(m_params);
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
		palign.setParams(m_params);
		palign.align();
		palign.postProcess();
		Logger::endTimer(1);
		output(palign);
	}
	else {
		Logger::beginTimer(1, "Multiple alignment");
		MultiAlign malign(m_chain_num);
		malign.setParams(m_params);
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
	string filename;
	vector<string> tokens;

	m_chain_num = m_filenames.size();
	m_pdbs.resize(m_chain_num);
	m_chains.resize(m_chain_num);

	if (!m_args.count("pocket")) {
		for (i=0; i<m_chain_num; i++) {
			m_chains[i].setRawName(m_filenames[i]);
			string_tokenize(tokens, m_filenames[i], ":");
			filename = tokens[0];
			if (filename.find('.') == string::npos) {
				filename += ".ent";
				if (filename.compare(0, 3, "pdb") != 0) {
					filename = "pdb" + filename;
				}
			}
			m_pdbs[i].readFile(filename);
			m_chains[i].setPDB(&m_pdbs[i]);
			if (tokens.size() > 1) {
				parseChainID(i, tokens[1]);
				if (tokens.size() > 2) {
					start = str2int(tokens[2]);
					if (tokens.size() > 3) end = str2int(tokens[3]);
					else end = 0;
					m_chains[i].setRange(start, end);
				}
			}
			m_chains[i].getChain();
			if (m_chains[i].length() == 0) exit(1);
		}
	}
	else {
		for (i=0; i<m_chain_num; i++) {
			m_chains[i].setRawName(m_filenames[i]);
			string_tokenize(tokens, m_filenames[i], ":");
			filename = tokens[0];
			if (filename.find('.') == string::npos) {
				filename += ".poc";
// 				if (filename.compare(0, 3, "pdb") != 0) {
// 					filename = "pdb" + filename;
// 				}
			}
			m_pdbs[i].readPocket(filename);
			m_chains[i].setPDB(&m_pdbs[i]);
			if (tokens.size() > 1) {
				parsePocketID(i, tokens[1]);
				if (tokens.size() > 2) parseChainID(i, tokens[2]);
			}
			m_chains[i].getPocketChain();
			if (m_chains[i].length() == 0) exit(1);
		}
	}
}

void Samo::parsePocketID(int i, const string &token)
{
	if (token[0] == '#') {
		m_chains[i].setPocketID(m_pdbs[i].getPocketID(str2int(token.substr(1))));
	}
	else {
		m_chains[i].setPocketID(str2int(token));
	}
}

void Samo::parseChainID(int i, const string &token)
{
	if (token[0] == '#') {
		m_chains[i].setChainID(m_pdbs[i].getChainID(str2int(token.substr(1))));
	}
	else {
		m_chains[i].setChainID(token[0]);
	}
}