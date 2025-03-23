## machinery for adding parameter transformations

## BMB: I'm a little puzzled about at what level (spec/calibrator)
## to do the operations. It is easy to do the transformations at
## the level of the specification, as long as the parameters to be
## transformed are *dynamical* parameters (i.e., relevant to the model
## spec, not e.g. dispersion parameters or other auxiliary parameters
## defined the purpose of model calibration)
##
## on the other hand, arguments aren't there in specs, we need the calibrator
## for that. For now I will write a separate function to modify the
## calibrator (after building the calibrator from the modified spec), but
## presumably these can be done together ??

## This machinery relies on mp_tmb_coef to do back-transformation

## FIXME: allow transformation table to be more transparent/augmented?
## add 'identity' to trans table?

#' data frame of transitions
#' name: label
#' link: link function (to be looked up via get())
#' invlink: inverse link function, as character value with %s in place of argument
trans_table <- data.frame(name = c("log", "logit"),
                          link = c("log", "qlogis"),
                          invlink = c("exp(%s)", "1/(1+exp(-%s))"))

mk_par_names <- function(trans_vec)  sprintf("%s_%s", trans_vec, names(trans_vec))

#' @param spec macpan2 model specification
#' @param trans_vec named vector
#' @examples
#' sir <- mp_tmb_library("starter_models", "sir", package = "macpan2")
#' mp_trans_pars(sir, c(beta = "log", gamma = "log"))
#' try(mp_trans_pars(sir, c(beta = "junk")))
mp_trans_pars <- function(spec, trans_vec) {
    trans_ind <- match(trans_vec, trans_table$name)
    ## FIXME: error for unmatched values
    if (any(is.na(trans_ind))) {
        stop("unknown transformation(s) requested: ",
             paste(trans_vec[is.na(trans_ind)], collapse = ", "))
    }
    trans_funs <- trans_table$link[trans_ind]
    trans_funs <- lapply(trans_funs, get, mode = "function")
    trans_pars <- mk_par_names(trans_vec)
    sstr <- trans_table$invlink[trans_ind]
    trans_ifuns <- sprintf(sstr, trans_pars)
    trans_exprs <- Map(function(p, ifun) reformulate(response = p, ifun),
                       names(trans_vec),
                       trans_ifuns)
    
    trans_defs <- Map(function(tr, p) tr(mp_default_list(spec)[[p]]),
                      trans_funs,
                      names(trans_vec))
    names(trans_defs) <- mk_par_names(trans_vec)
    trans_spec <- mp_tmb_insert(spec, phase = "before"
                              , at = 1
                              , expressions = trans_exprs,
                              , default = trans_defs
                                )
    return(trans_spec)
}

mp_trans_args <- function(cal, trans_vec) {
    browser()
    trans_pars <- sprintf("%s_%s", trans_vec, names(trans_vec))
    new_pars <- cal$cal_args$par
    new_pars[match(names(trans_vec), new_pars)] <- trans_pars
    cal$cal_args$par <- new_pars
    return(cal)
}

