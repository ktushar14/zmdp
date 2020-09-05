/********** tell emacs we use -*- c++ -*- style comments *******************
 $Revision: 1.5 $  $Author: trey $  $Date: 2007-04-08 22:48:04 $

 @file    FastParser.cc
 @brief   No brief

 Copyright (c) 2002-2005, Trey Smith. All rights reserved.

 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License.  You may
 obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 implied.  See the License for the specific language governing
 permissions and limitations under the License.

 ***************************************************************************/

/***************************************************************************
 * INCLUDES
 ***************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <string>

#include "FastParser.h"
#include "MatrixUtils.h"
#include "slaMatrixUtils.h"
#include "sla_cassandra.h"
#include "zmdpCommonDefs.h"

#define POMDP_READ_ERROR_EPS (1e-10)

using namespace std;
using namespace MatrixUtils;

namespace zmdp {

/***************************************************************************
 * STATIC HELPER FUNCTIONS
 ***************************************************************************/

static void trimTrailingWhiteSpace(char* s)
{
	int n = strlen(s);
	int i;
	for (i = n - 1; i >= 0; i--) {
		if (!isspace(s[i])) break;
	}
	s[i + 1] = '\0';
}

/***************************************************************************
 * POMDP FUNCTIONS
 ***************************************************************************/

void FastParser::readGenericDiscreteMDPFromFile(CassandraModel& mdp, const std::string& _fileName)
{
	mdp.fileName = _fileName;
	readModelFromFile(mdp, /* expectPomdp = */ false);
}

void FastParser::readPomdpFromFile(CassandraModel& pomdp, const std::string& _fileName)
{
	pomdp.fileName = _fileName;
	readModelFromFile(pomdp, /* expectPomdp = */ true);
}

void FastParser::readModelFromFile(CassandraModel& problem, bool expectPomdp)
{
	ifstream in;

	timeval startTime, endTime;
	if (zmdpDebugLevelG >= 1) {
		cout << "reading problem (in fast mode) from " << problem.fileName << endl;
		gettimeofday(&startTime, 0);
	}

	in.open(problem.fileName.c_str());
	if (!in) {
		cerr << "ERROR: couldn't open " << problem.fileName << " for reading: " << strerror(errno)
			 << endl;
		exit(EXIT_FAILURE);
	}

	readModelFromStream(problem, in, expectPomdp);

	in.close();

    std::cout << "  closed stream!\n";

	if (zmdpDebugLevelG >= 1) {
		gettimeofday(&endTime, 0);
		double numSeconds =
			(endTime.tv_sec - startTime.tv_sec) + 1e-6 * (endTime.tv_usec - startTime.tv_usec);
		cout << "[file reading took " << numSeconds << " seconds]" << endl;
	}
}

bool StrMatchesStartPrefix(std::string& line)
{
    std::string prefix = "start:";
    std::cout << "\tIn StrMatchesStartPrefix\n";
    std::cout << "\tline:   " << line << "\n";
    std::cout << "\tprefix: " << prefix << "\n";
    return 0 == strncmp(line.c_str(), (prefix.c_str()), strlen(prefix.c_str()));
}

void FastParser::readModelFromStream(CassandraModel& p, std::istream& in, bool expectPomdp)
{
	char buf[1 << 20];
	int lineNumber;
	char sbuf[512];
	bool inPreamble = true;

	kmatrix Rx;
	std::vector<kmatrix> Tx, Ox;

	bool discountSet = false;
	bool valuesSet = false;
	bool numStatesSet = false;
	bool numActionsSet = false;
	bool numObservationsSet = false;
	bool startSet = false;

	const char* rFormat = (expectPomdp ? "R: %d : %d : * : * %lf" : "R: %d : %d : * %lf");

#define PM_PREFIX_MATCHES(X) (0 == strncmp(buf, (X), strlen(X)))

	lineNumber = 1;
	std::string line;
	bool nextFieldIsStart = false;

	while (!in.eof()) {
		std::getline(in, line);
		int numCharsToUnget = (int)line.size() + 1;

		if (!nextFieldIsStart) {
			while (numCharsToUnget--) {
				in.unget();
			}
			in.getline(buf, sizeof(buf));
			if (line.empty()) {
				// std::cout << "  1 std::getline read empty string!" << std::endl;
			}
			// std::cout << "  1 line: " << line << std::endl;
            // std::cout << "  1 buf:  " << buf << std::endl;
		} else {
			if (line.empty()) {
                continue;
            }
            // std::cout << "Using std::getline result\n";
            // memset(buf, 0, sizeof(buf));
            strcpy(buf, "start:");
			// std::cout << "  2 buf: " << buf << std::endl;
            // std::cout << "  2 line (sz): " << line.size() << std::endl;
		}

		if (in.fail() && !in.eof()) {
            if (line.empty()) {
				cerr << "I think this is the end of the file because line is empty(); breaking out "
						"of the while loop\n";
				break;
            } else {
                cerr << "ERROR: " << p.fileName << ": line " << lineNumber
    				 << ": line too long for buffer"
    				 << " (max length " << sizeof(buf) << ")" << endl;
                std::cout << "  line: " << line << std::endl;
                std::cout << "  in.eof(): " << in.eof() << std::endl;
    			exit(EXIT_FAILURE);
            }
		}

		if ('#' == buf[0]) {
            lineNumber++;
            continue;
        }
		trimTrailingWhiteSpace(buf);
		if ('\0' == buf[0]) {
            lineNumber++;
            continue;
        }

        if (inPreamble) {
			if (PM_PREFIX_MATCHES("discount:")) {
                // std::cout << "checking discount\n";
                // std::cout << "\tbuf: " << buf << "\n";
				if (1 != sscanf(buf, "discount: %lf", &p.discount)) {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": syntax error in 'discount' statement" << endl;
					exit(EXIT_FAILURE);
				}
				discountSet = true;
			} else if (PM_PREFIX_MATCHES("values:")) {
                // std::cout << "checking values\n";
                // std::cout << "\tbuf: " << buf << "\n";
				if (1 != sscanf(buf, "values: %s", sbuf)) {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": syntax error in 'values' statement" << endl;
					exit(EXIT_FAILURE);
				}
				if (0 != strcmp(sbuf, "reward")) {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": expected 'values: reward', other types not supported by fast parser"
						 << endl;
					exit(EXIT_FAILURE);
				}
				valuesSet = true;
			} else if (PM_PREFIX_MATCHES("actions:")) {
                // std::cout << "checking actions\n";
                // std::cout << "\tbuf: " << buf << "\n";
				if (1 != sscanf(buf, "actions: %d", &p.numActions)) {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": syntax error in 'actions' statement" << endl;
					exit(EXIT_FAILURE);
				}
				numActionsSet = true;
			} else if (PM_PREFIX_MATCHES("observations:")) {
                // std::cout << "checking observations\n";
                // std::cout << "\tbuf: " << buf << "\n";
				if (expectPomdp) {
					if (1 != sscanf(buf, "observations: %d", &p.numObservations)) {
						cerr << "ERROR: " << p.fileName << ": line " << lineNumber
							 << ": syntax error in 'observations' statement" << endl;
						exit(EXIT_FAILURE);
					}
					numObservationsSet = true;
				} else {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": got unexpected 'observations' statement in MDP" << endl;
					exit(EXIT_FAILURE);
				}
			} else if (PM_PREFIX_MATCHES("states:")) {
                // std::cout << "checking states\n";
                // std::cout << "\tbuf: " << buf << "\n";
				if (1 != sscanf(buf, "states: %d", &p.numStates)) {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": syntax error in 'states' statement" << endl;
					exit(EXIT_FAILURE);
				}
				numStatesSet = true;
				nextFieldIsStart = true;
            } else if (PM_PREFIX_MATCHES("start:")) {
			// } else if (StrMatchesStartPrefix(line)) {
                // std::cout << "CHECKING START\n";
                // std::cout << "\tline(sz): " << line.size() << "\n";
                // std::cout << "\tbuf:  " << buf << "\n";
				if (!numStatesSet) {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": got 'start' statement before 'states' statement" << endl;
					exit(EXIT_FAILURE);
				}
				// readStartVector(p, buf, expectPomdp);
				readStartVector(p, const_cast<char*>(line.c_str()), expectPomdp);
				startSet = true;
				nextFieldIsStart = false;
                std::cout << "  Start set!\n";
			} else {
				// the statement is not one that is expected in the preamble,
				// check that we are ready to transition to parsing the body

#define FP_CHECK_SET(VAR, NAME)                                                                    \
	if (!(VAR)) {                                                                                  \
		cerr << "ERROR: " << p.fileName << ": line " << lineNumber << ": at end of preamble, no '" \
			 << (NAME) << "' statement found" << endl;                                             \
	}

				FP_CHECK_SET(discountSet, "discount");
				FP_CHECK_SET(valuesSet, "values");
				FP_CHECK_SET(numStatesSet, "states");
				FP_CHECK_SET(numActionsSet, "actions");
				FP_CHECK_SET(startSet, "start");
				if (expectPomdp) {
					FP_CHECK_SET(numObservationsSet, "observations");
				} else {
					p.numObservations = -1;
				}

				// initialize data structures
				Rx.resize(p.numStates, p.numActions);
				Tx.resize(p.numActions);
				if (expectPomdp) {
					Ox.resize(p.numActions);
				}
				FOR(a, p.numActions)
				{
					Tx[a].resize(p.numStates, p.numStates);
					if (expectPomdp) {
						Ox[a].resize(p.numStates, p.numObservations);
					}
				}

				// henceforth expect body statements instead of preamble statements
				inPreamble = false;
			}
		}

		if (!inPreamble) {
			if (PM_PREFIX_MATCHES("R:")) {
                // std::cout << "    " << buf << "lineNumber: " << lineNumber << "\n";
				int s, a;
				double reward;
				if (3 != sscanf(buf, rFormat, &a, &s, &reward)) {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": syntax error in R statement\n"
						 << "  (expected format is '" << rFormat << "')" << endl;
					exit(EXIT_FAILURE);
				}
				kmatrix_set_entry(Rx, s, a, reward);
			} else if (PM_PREFIX_MATCHES("T:")) {
                // std::cout << "    " << buf << "lineNumber: " << lineNumber << "\n";
				int s, a, sp;
				double prob;
				if (4 != sscanf(buf, "T: %d : %d : %d %lf", &a, &s, &sp, &prob)) {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": syntax error in T statement" << endl;
					exit(EXIT_FAILURE);
				}
				kmatrix_set_entry(Tx[a], s, sp, prob);
			} else if (PM_PREFIX_MATCHES("O:")) {
                // std::cout << "    " << buf << "lineNumber: " << lineNumber << "\n";
				if (expectPomdp) {
					int s, a, o;
					double prob;
					if (4 != sscanf(buf, "O: %d : %d : %d %lf", &a, &s, &o, &prob)) {
						cerr << "ERROR: " << p.fileName << ": line " << lineNumber
							 << ": syntax error in O statement" << endl;
						exit(EXIT_FAILURE);
					}
					kmatrix_set_entry(Ox[a], s, o, prob);
				} else {
					cerr << "ERROR: " << p.fileName << ": line " << lineNumber
						 << ": got unexpected 'O' statement in MDP" << endl;
					exit(EXIT_FAILURE);
				}
			} else {
				cerr << "ERROR: " << p.fileName << ": line " << lineNumber
					 << ": got unexpected statement type while parsing body" << endl;
				exit(EXIT_FAILURE);
			}
		}

		lineNumber++;
	}

	std::cout << "  It's here\n";

    // post-process
	copy(p.R, Rx);
	Rx.clear();

	p.T.resize(p.numActions);
	p.Ttr.resize(p.numActions);
	if (expectPomdp) {
		p.O.resize(p.numActions);
	}
	FOR(a, p.numActions)
	{
		copy(p.T[a], Tx[a]);
		kmatrix_transpose_in_place(Tx[a]);
		copy(p.Ttr[a], Tx[a]);

#if 1
		// extra error checking
		cvector checkTmp;
		FOR(s, p.numStates)
		{
			copy_from_column(checkTmp, p.Ttr[a], s);
			if (fabs(sum(checkTmp) - 1.0) > POMDP_READ_ERROR_EPS) {
				fprintf(stderr,
						"ERROR: %s: outgoing transition probabilities do not sum to 1 for:\n"
						"  state %d, action %d, transition sum = %.10lf\n",
						p.fileName.c_str(),
						(int)s,
						(int)a,
						sum(checkTmp));
				exit(EXIT_FAILURE);
			}
		}
#endif

		Tx[a].clear();

		if (expectPomdp) {
			copy(p.O[a], Ox[a]);

#if 1
			cmatrix checkObs;

			// extra error checking
			kmatrix_transpose_in_place(Ox[a]);
			copy(checkObs, Ox[a]);

			FOR(s, p.numStates)
			{
				copy_from_column(checkTmp, checkObs, s);
				if (fabs(sum(checkTmp) - 1.0) > POMDP_READ_ERROR_EPS) {
					fprintf(stderr,
							"ERROR: %s: observation probabilities do not sum to 1 for:\n"
							"  state %d, action %d, observation sum = %.10lf\n",
							p.fileName.c_str(),
							(int)s,
							(int)a,
							sum(checkTmp));
					exit(EXIT_FAILURE);
				}
			}
#endif

			Ox[a].clear();
		}
	}

	p.checkForTerminalStates();

#if 1
	// extra error checking
	if (expectPomdp) {
		if (fabs(sum(p.initialBelief) - 1.0) > POMDP_READ_ERROR_EPS) {
			fprintf(stderr,
					"ERROR: %s: initial belief entries do not sum to 1:\n"
					"  entry sum = %.10lf\n",
					p.fileName.c_str(),
					sum(p.initialBelief));
			exit(EXIT_FAILURE);
		}
	}
#endif

	if (zmdpDebugLevelG >= 1) {
		p.debugDensity();
	}
}

void FastParser::readStartVector(CassandraModel& p, char* data, bool expectPomdp)
{
	std::cout << "   FastParser::readStartVector:\n";
	// std::cout << "   data:\n";
	// std::cout << "   " << data << "\n";

	if (expectPomdp) {
		// POMDP case

		int i;
		char* tok;
		double val;

		p.initialBelief.resize(p.numStates);

		// consume 'start:' token at the beginning of the statement
		tok = strtok(data, " ");

		for (i = 0; i < p.numStates; i++) {
			if (NULL != tok) {
				tok = strtok(NULL, " ");
			}
			if (NULL == tok) {
				if (0 == i) {
					double startState = atof(tok);
					if (startState != floor(startState)) {
						cout << "ERROR: " << p.fileName
							 << ": POMDP 'start' statement must either contain a single integer "
							 << "specifying a known start state or a list of the initial "
								"probabilities of "
							 << "all states" << endl;
						exit(EXIT_FAILURE);
					}
					p.initialBelief.push_back((int)startState, 1.0);
					return;
				} else {
					cout << "ERROR: " << p.fileName
						 << ": POMDP 'start' statement must either contain a single integer "
						 << "specifying a known start state or a list of the initial probabilities "
							"of "
						 << "all states" << endl;
					exit(EXIT_FAILURE);
				}
			}

			val = atof(tok);
			if (val > SPARSE_EPS) {
				p.initialBelief.push_back(i, val);
			}
		}
	} else {
		// MDP case

		double x, y;
		int ret = sscanf(data, "start: %lf %lf", &x, &y);
		if ((ret != 1) || (x != floor(x))) {
			cout << "ERROR: " << p.fileName
				 << ": MDP 'start' statement must contain a single integer "
				 << "specifying a known start state" << endl;
			exit(EXIT_FAILURE);
		}

		p.initialState.resize(1);
		int startState;
		if (1 != sscanf(data, "start: %d", &startState)) {
			cout << "ERROR: " << p.fileName
				 << ": MDP 'start' statement must contain a single integer "
				 << "specifying a known start state" << endl;
			exit(EXIT_FAILURE);
		}
		p.initialState.push_back(0, startState);
	}
}

}; // namespace zmdp

/***************************************************************************
 * REVISION HISTORY:
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2007/01/15 17:23:34  trey
 * fix problem that was causing zeros to have entries in the sparse representation of initialBelief
 *
 * Revision 1.3  2006/11/09 21:11:33  trey
 * removed obsolete code referencing variable initialBeliefx
 *
 * Revision 1.2  2006/11/09 20:48:56  trey
 * fixed some MDP vs. POMDP issues
 *
 * Revision 1.1  2006/11/08 16:40:50  trey
 * initial check-in
 *
 *
 ***************************************************************************/
