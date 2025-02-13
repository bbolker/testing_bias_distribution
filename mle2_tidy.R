## improvement on tidy method in broom package (allows different conf methods)
tidy.mle2 <- function (x, conf.int = FALSE, conf.level = 0.95,
          conf.method = "quad", ...) 
{
    ## check_ellipses("exponentiate", "tidy", "mle2", ...)
    co <- bbmle::coef(bbmle::summary(x))
    ret <- broom:::as_tidy_tibble(co, new_names = c("estimate", 
        "std.error", "statistic", "p.value"))
    if (conf.int) {
        ci <- bbmle::confint(x, level = conf.level, method = conf.method)
        if (is.null(dim(ci))) {
            ci <- matrix(ci, nrow = 1)
        }
        colnames(ci) <- c("conf.low", "conf.high")
        ci <- as_tibble(ci)
        ret <- dplyr::bind_cols(ret, ci)
    }
    as_tibble(ret)
}
