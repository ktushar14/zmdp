/********** tell emacs we use -*- c++ -*- style comments *******************
 $Revision: 1.9 $  $Author: trey $  $Date: 2006-02-01 01:09:37 $
   
 @file    Solver.h
 @brief   No brief

 Copyright (c) 2002-2005, Trey Smith
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 * The software may not be sold or incorporated into a commercial
   product without specific prior written permission.
 * The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 ***************************************************************************/

#ifndef INCSolver_h
#define INCSolver_h

#include "MDP.h"

namespace zmdp {

class Solver {
public:
  virtual ~Solver(void) {}

  // sets up the problem
  virtual void planInit(const MDP* problem) = 0;

  // plan for a fixed amount of time.  if maxTimeSeconds < 0,
  //   the amount of time is chosen by the solver to optimize
  //   time performance.  returns true if minPrecision has been
  //   reached.
  virtual bool planFixedTime(const state_vector& currentState,
			     double maxTimeSeconds,
			     double minPrecision) = 0;

  virtual int chooseAction(const state_vector& currentState) = 0;

  virtual void setBoundsFile(std::ostream* boundsFile) = 0;
  virtual ValueInterval getValueAt(const state_vector& currentState) const = 0;

  // sets the minimum safety value, for a solver that understands safety
  virtual void setMinSafety(double _minSafety) {}
  
};

}; // namespace zmdp

#endif // INCSolver_h

/***************************************************************************
 * REVISION HISTORY:
 * $Log: not supported by cvs2svn $
 * Revision 1.8  2006/01/28 22:01:10  trey
 * switched include PomdpSim.h -> MDP.h
 *
 * Revision 1.7  2006/01/28 03:07:05  trey
 * improved flexibility for use with mdps
 *
 * Revision 1.6  2005/11/28 20:45:47  trey
 * fixed warning about non-virtual destructor
 *
 * Revision 1.5  2005/10/28 03:50:32  trey
 * simplified license
 *
 * Revision 1.4  2005/10/28 02:51:40  trey
 * added copyright headers
 *
 * Revision 1.3  2005/10/21 20:09:11  trey
 * added namespace zmdp
 *
 * Revision 1.2  2005/01/21 18:07:02  trey
 * preparing for transition to sla matrix types
 *
 * Revision 1.1  2004/11/13 23:29:44  trey
 * moved many files from hsvi to common
 *
 * Revision 1.1.1.1  2004/11/09 16:18:56  trey
 * imported hsvi into new repository
 *
 * Revision 1.1  2003/09/23 21:11:51  trey
 * initial check-in
 *
 * Revision 1.1  2001/08/27 17:49:16  trey
 * initial check-in
 *
 *
 ***************************************************************************/
