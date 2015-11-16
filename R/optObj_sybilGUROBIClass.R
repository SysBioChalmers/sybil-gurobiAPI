#  optObj_sybilGUROBIClass.R
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


#------------------------------------------------------------------------------#
#               definition of the class optObj_sybilGUROBI                     #
#------------------------------------------------------------------------------#

setClass(Class = "optObj_sybilGUROBI",
         representation(grb = "character"),
         contains = "optObj")


#------------------------------------------------------------------------------#
#                                  methods                                     #
#------------------------------------------------------------------------------#

setMethod("delProb", signature(lp = "optObj_sybilGUROBI"),

    function(lp, ...) {

        .GRBenv[[lp@grb]] <- NULL

    }
)


#------------------------------------------------------------------------------#

setMethod("initProb", signature(lp = "optObj_sybilGUROBI"),

    function(lp, to = FALSE, nrows = 0, ncols = 0) {

        prob        <- vector(mode = "list", length = 8)
        names(prob) <- c("A", "obj", "sense", "rhs", "lb", "ub",
                         "modelsense", "modelname")

        repeat{
            pn <- paste(sample(letters, 7), collapse = "")
            if ( (is.null(.GRBenv)) || (! pn %in% ls(.GRBenv)) ) {
                break
            }
        }
        
        lp@grb <- pn

        if (isTRUE(to)) {
            toflag <- 1
        }
        else {
            toflag <- 0
        }

        .GRBenv[[lp@grb]] <- list(lp   = prob,
                                  parm = list(OutputFlag = toflag),
                                  sol  = vector(mode = "list", length = 0))

        .GRBenv[[lp@grb]][["lp"]][["A"]]           <- Matrix(0,
                                                             nrow = nrows,
                                                             ncol = ncols,
                                                             sparse = TRUE)
        .GRBenv[[lp@grb]][["lp"]][["obj"]]         <- numeric(ncols)
        .GRBenv[[lp@grb]][["lp"]][["sense"]]       <- rep("=", nrows)
        .GRBenv[[lp@grb]][["lp"]][["rhs"]]         <- numeric(nrows)
        .GRBenv[[lp@grb]][["lp"]][["lb"]]          <- numeric(ncols)
        .GRBenv[[lp@grb]][["lp"]][["ub"]]          <- numeric(ncols)
        .GRBenv[[lp@grb]][["lp"]][["modelsense"]]  <- "min"
        .GRBenv[[lp@grb]][["lp"]][["modelname"]]   <- lp@grb

        return(lp)
    }
)


#------------------------------------------------------------------------------#

setMethod("backupProb", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        repeat{
            pn <- paste(sample(letters, 7), collapse = "")
            if ( (is.null(.GRBenv)) || (! pn %in% ls(.GRBenv)) ) {
                break
            }
        }

        out <- new("optObj_sybilGUROBI", lp@solver, lp@method, lp@probType)

        out@grb <- pn

        .GRBenv[[out@grb]] <- .GRBenv[[lp@grb]]

        return(out)
    }
)


#------------------------------------------------------------------------------#


setMethod("setSolverParm", signature(lp = "optObj_sybilGUROBI"),

    function(lp, solverParm) {

        if ( ! ((is.data.frame(solverParm)) || (is.list(solverParm))) ) {
            stop(sQuote(solverParm), " must be list or data.frame")
        }

        if (any(is.na(solverParm))) {
            stop(sQuote(solverParm), " contains NA values")
        }

        pn <- names(solverParm)
        for (i in seq(along = solverParm)) {
            .GRBenv[[lp@grb]][["parm"]][[pn[i]]] <- solverParm[[pn[i]]]
        }
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolverParm", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        out <- .GRBenv[[lp@grb]][["parm"]]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_sybilGUROBI",
                                 lpdir = "character"),

    function(lp, lpdir) {

        .GRBenv[[lp@grb]][["lp"]][["modelsense"]] <- ifelse(lpdir == "max",
                                                               "max", "min")

    }
)


#------------------------------------------------------------------------------#

setMethod("getObjDir", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        out <- .GRBenv[[lp@grb]][["lp"]][["modelsense"]]

        return(out)

    }
)


#------------------------------------------------------------------------------#

setMethod("getNumRows", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        out <- nrow(.GRBenv[[lp@grb]][["lp"]][["A"]])

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumCols", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        out <- ncol(.GRBenv[[lp@grb]][["lp"]][["A"]])

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsToProb", signature(lp = "optObj_sybilGUROBI"),

    # i: vector containing the new row indices (must be ascending)
    # cind: list, containing the column indices of the new nz elements
    # nzval: list, containing the new nz elements
    #
    # i, type, lb, ub, cind and nzval must have the same length
    #
    # type can be one of the following:
    # "F" = free variable                -INF <  x <  INF
    # "L" = variable with lower bound      lb <= x <  INF
    # "U" = variable with upper bound    -INF <  x <= ub
    # "D" = double-bounded variable        lb <= x <= ub
    # "E" = fixed variable                 lb  = x  = ub

    function(lp, i, type, lb, ub, cind, nzval, rnames = NULL) {

        stopifnot(length(lb) == length(ub))
        ubc      <- type == "U"
        clb      <- lb
        clb[ubc] <- ub[ubc]
        
        nc  <- ncol(.GRBenv[[lp@grb]][["lp"]][["A"]])
        nr  <- nrow(.GRBenv[[lp@grb]][["lp"]][["A"]])
        mat <- Matrix(0, nrow = length(i), ncol = nc)
        
        .GRBenv[[lp@grb]][["lp"]][["A"]] <- rBind(
                                       .GRBenv[[lp@grb]][["lp"]][["A"]], mat)

        .GRBenv[[lp@grb]][["lp"]][["rhs"]] <- append(
                      .GRBenv[[lp@grb]][["lp"]][["rhs"]], numeric(length(i)))

        .GRBenv[[lp@grb]][["lp"]][["sense"]] <- append(
                   .GRBenv[[lp@grb]][["lp"]][["sense"]], rep("=", length(i)))

        
        for (k in seq(along = i)) {
            ltype <- switch(EXPR = type[k],
                            "L" = { ">=" },
                            "U" = { "<=" },
                            "E" = { "=" },
                                  { "=" }
            )

            .GRBenv[[lp@grb]][["lp"]][["A"]][(nr+k), cind[[k]]] <- nzval[[k]]
    
            .GRBenv[[lp@grb]][["lp"]][["rhs"]][(nr+k)] <- clb[k]
    
            .GRBenv[[lp@grb]][["lp"]][["sense"]][(nr+k)] <- ltype
        }


    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBnds", signature(lp = "optObj_sybilGUROBI"),

    function(lp, j, lb, ub) {

        .GRBenv[[lp@grb]][["lp"]][["lb"]][j] <- lb
        .GRBenv[[lp@grb]][["lp"]][["ub"]][j] <- ub
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBndsObjCoefs", signature(lp = "optObj_sybilGUROBI"),

    function(lp, j, lb, ub, obj_coef) {

        .GRBenv[[lp@grb]][["lp"]][["lb"]][j]  <- lb
        .GRBenv[[lp@grb]][["lp"]][["ub"]][j]  <- ub
        .GRBenv[[lp@grb]][["lp"]][["obj"]][j] <- obj_coef

    }
)


#------------------------------------------------------------------------------#

setMethod("getColsLowBnds", signature(lp = "optObj_sybilGUROBI"),

    function(lp, j) {

        out <- .GRBenv[[lp@grb]][["lp"]][["lb"]][j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsUppBnds", signature(lp = "optObj_sybilGUROBI"),

    function(lp, j) {

        out <- .GRBenv[[lp@grb]][["lp"]][["ub"]][j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeRowsBnds", signature(lp = "optObj_sybilGUROBI"),

    function(lp, i, lb, ub) {

        ubc      <- .GRBenv[[lp@grb]][["lp"]][["sense"]][i] == "<="
        clb      <- lb
        clb[ubc] <- ub[ubc]

        .GRBenv[[lp@grb]][["lp"]][["rhs"]][i] <- clb
    }
)


#------------------------------------------------------------------------------#

setMethod("setRhsZero", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        nr  <- nrow(.GRBenv[[lp@grb]][["lp"]][["A"]])
        .GRBenv[[lp@grb]][["lp"]][["rhs"]]   <- rep(0, nr)
        .GRBenv[[lp@grb]][["lp"]][["sense"]] <- rep("=", nr)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeObjCoefs", signature(lp = "optObj_sybilGUROBI"),

    function(lp, j, obj_coef) {

        .GRBenv[[lp@grb]][["lp"]][["obj"]][j] <- obj_coef
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjCoefs", signature(lp = "optObj_sybilGUROBI"),

    function(lp, j) {

        out <- .GRBenv[[lp@grb]][["lp"]][["obj"]][j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeMatrixRow", signature(lp = "optObj_sybilGUROBI"),

    function(lp, i, j, val) {

        .GRBenv[[lp@grb]][["lp"]][["A"]][i, j] <- val

    }
)


#------------------------------------------------------------------------------#

setMethod("loadLPprob", signature(lp = "optObj_sybilGUROBI"),

    function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype,
             lpdir = "max", rub = NULL, ctype = NULL,
             cnames = NULL, rnames = NULL, pname = NULL) {

        stopifnot(is(mat, "Matrix"))

        crtype <- sapply(rtype,
                         function(x) switch(EXPR = x,
                                            "L" = { ">=" },
                                            "U" = { "<=" },
                                            "E" = { "="  },
                                                  { "="  }))

        # optimization direction
        setObjDir(lp, lpdir = lpdir)
        
        # constraint matrix
        .GRBenv[[lp@grb]][["lp"]][["A"]] <- as(mat, "CsparseMatrix")

        # column (variable) bounds and objective function
        .GRBenv[[lp@grb]][["lp"]][["lb"]]  <- lb
        .GRBenv[[lp@grb]][["lp"]][["ub"]]  <- ub
        .GRBenv[[lp@grb]][["lp"]][["obj"]] <- obj

        # variable type
        .GRBenv[[lp@grb]][["lp"]][["vtypes"]] <- ctype

        # right hand side
        if (is.null(rub)) {
            crlb      <- rlb
        }
        else {
            ubc       <- rtype == "U"
            crlb      <- rlb
            crlb[ubc] <- rub[ubc]
        }
        
        .GRBenv[[lp@grb]][["lp"]][["sense"]] <- crtype
        .GRBenv[[lp@grb]][["lp"]][["rhs"]]   <- crlb

        # problem name
        if (!is.null(pname)) {
            .GRBenv[[lp@grb]][["lp"]][["modelname"]] <- pname
        }

    }
)


#------------------------------------------------------------------------------#

setMethod("loadQobj", signature(lp = "optObj_sybilGUROBI", mat = "Matrix"),

    function(lp, mat) {

        .GRBenv[[lp@grb]][["lp"]][["Q"]] <- as(mat, "CsparseMatrix")

    }
)


#------------------------------------------------------------------------------#

setMethod("loadQobj", signature(lp = "optObj_sybilGUROBI", mat = "numeric"),

    function(lp, mat) {

        dmat <- Matrix::Diagonal(length(mat), mat)

        .GRBenv[[lp@grb]][["lp"]][["Q"]] <- as(dmat, "CsparseMatrix")

    }
)


#------------------------------------------------------------------------------#

setMethod("scaleProb", signature(lp = "optObj_sybilGUROBI"),

    function(lp, opt) {

        .GRBenv[[lp@grb]][["parm"]][["ScaleFlag"]] <- opt

    }
)


#------------------------------------------------------------------------------#

setMethod("solveLp", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        if (length(.GRBenv[[lp@grb]][["parm"]]) > 0) {
            gp <- .GRBenv[[lp@grb]][["parm"]]
        }
        else {
            gp <- NULL
        }

        .GRBenv[[lp@grb]][["sol"]] <- gurobi::gurobi(.GRBenv[[lp@grb]][["lp"]],gp)

        .GRBenv[[lp@grb]][["sol"]][["status"]] <- GRB_STATUS[[.GRBenv[[lp@grb]][["sol"]][["status"]]]]
                
        if (is.null(.GRBenv[[lp@grb]][["sol"]][["objval"]])) {
            out <- 1
            .GRBenv[[lp@grb]][["sol"]][["objval"]] <- 0
        }
        else {
            out <- 0
        }
        
#        if (!is.null(.GRBenv[[lp@grb]][["lp"]][["advStartSet"]])) {
            if (!is.null(.GRBenv[[lp@grb]][["sol"]][["vbasis"]])) {
                .GRBenv[[lp@grb]][["lp"]][["vbasis"]] <- .GRBenv[[lp@grb]][["sol"]][["vbasis"]]
            }
            if (!is.null(.GRBenv[[lp@grb]][["sol"]][["cbasis"]])) {
                .GRBenv[[lp@grb]][["lp"]][["cbasis"]] <- .GRBenv[[lp@grb]][["sol"]][["cbasis"]]
            }
#            .GRBenv[[lp@grb]][["lp"]][["advStartSet"]] <- TRUE
#        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjVal", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        out <- .GRBenv[[lp@grb]][["sol"]][["objval"]]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRedCosts", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        out <- .GRBenv[[lp@grb]][["sol"]][["rc"]]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolStat", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        out <- .GRBenv[[lp@grb]][["sol"]][["status"]]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getFluxDist", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        fld <- .GRBenv[[lp@grb]][["sol"]][["x"]]
        if (is.null(fld)) {
            out <- rep(0, getNumCols(lp))
        }
        else {
            out <- fld
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColPrim", signature(lp = "optObj_sybilGUROBI"),

    function(lp, j) {

        out <- .GRBenv[[lp@grb]][["sol"]][["x"]][j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumNnz", signature(lp = "optObj_sybilGUROBI"),

    function(lp) {

        out <- nnzero(.GRBenv[[lp@grb]][["lp"]][["A"]])

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("writeProb", signature(lp = "optObj_sybilGUROBI"),

    function(lp, fname, ff = "lp") {

        if (length(.GRBenv[[lp@grb]][["parm"]]) > 0) {
            gp <- .GRBenv[[lp@grb]][["parm"]]
        }
        else {
            gp <- vector(mode = "list", length = 0)

        }

        if (grepl(".", fname, fixed = TRUE)) {
            fn <- fname
        }
        else {
            fn <- paste(fname, ff, sep = ".")
        }

        gp[["ResultFile"]] <- fn
        gurobi::gurobi(.GRBenv[[lp@grb]][["lp"]], gp)


        return(TRUE)


        #out <- save(.GRBenv[[lp@grb]][["lp"]], file = fname)

        #return(out)
    }
)


#------------------------------------------------------------------------------#


























