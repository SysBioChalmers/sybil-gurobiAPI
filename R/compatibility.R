#  compatibility.R
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


GRB_STATUS <- list(
    "LOADED"          = 1,
    "OPTIMAL"         = 2,
    "INFEASIBLE"      = 3,
    "INF_OR_UNBD"     = 4,
    "UNBOUNDED"       = 5,
    "CUTOFF"          = 6,
    "ITERATION_LIMIT" = 7,
    "NODE_LIMIT"      = 8,
    "TIME_LIMIT"      = 9,
    "SOLUTION_LIMIT"  = 10,
    "INTERRUPTED"     = 11,
    "NUMERIC"         = 12,
    "SUBOPTIMAL"      = 13
)

checkSolutionStatus <- function(stat) {
    out <- which(stat != 2)
    return(out)
}

getReturnString <- function(code) {
    out <- ifelse(code == 0, "ok", "not ok")
    return(out)
}

getStatusString <- function(code) {
    out <- names(GRB_STATUS)[code]
    return(out)
}

