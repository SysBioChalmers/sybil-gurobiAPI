#  zzz.R
#  Gurobi support for sybil.
#
#  Copyright (C) 2010-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybilGUROBI.
#
#  SybilGUROBI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SybilGUROBI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybilGUROBI.  If not, see <http://www.gnu.org/licenses/>.


.GRBenv <- new.env()

.onAttach <- function(lib, pkg) {

    sybil::addSolver(solver = "sybilGUROBI",
                     method = "gurobi",
                     probType = list(c("lp", "mip", "qp")))
    
#    SYBIL_SETTINGS("SOLVER_CTRL_PARM",
#                   list(OutputFlag = 0,
#                   # i can not get rid of the logfile completely, so I put
#                   # it to temp; an initial logfile will still be written
#                   # into the working directory
#                   LogFile = tempfile(pattern = "gurobi", fileext = ".log"))
#                   )
    
}


