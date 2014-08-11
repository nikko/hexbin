"xy.hexplot" <-
function (formula, data = parent.frame(), allow.multiple = FALSE, 
    outer = FALSE, auto.key = FALSE, aspect = "fill", layout = NULL,
    panel = "panel.hexplot", 
    prepanel = NULL, scales = list(), strip = TRUE,
    xbins=25, style =c("colorscale", "centroids", "lattice",
                 "nested.lattice", "nested.centroids"),
    minarea = 0.05, maxarea = 0.8,mincnt = 1, maxcnt = NULL,
    trans = NULL, colorcut, density = NULL, border = FALSE, 
    pen = NULL, colramp=BTY, xbnds = range(x), ybnds = range(y)
    xlab, xlim, ylab, ylim, ..., subset = TRUE) 
{
    dots <- list(...)
    ##groups <- eval(substitute(groups), data, parent.frame())
    subset <- eval(substitute(subset), data, parent.frame())
    form <- latticeParseFormula(formula, data, subset = subset, 
        groups = groups, multiple = allow.multiple, outer = outer, 
        subscripts = TRUE)
    groups <- form$groups
    if (!is.function(panel)) 
        panel <- eval(panel)
    if (!is.function(strip)) 
        strip <- eval(strip)
    if ("subscripts" %in% names(formals(panel))) 
        subscripts <- TRUE
    if (subscripts) 
        subscr <- form$subscr
    prepanel <- if (is.function(prepanel)) 
        prepanel
    else if (is.character(prepanel)) 
        get(prepanel)
    else eval(prepanel)
    cond <- form$condition
    number.of.cond <- length(cond)
    y <- form$left
    x <- form$right
    if (number.of.cond == 0) {
        strip <- FALSE
        cond <- list(as.factor(rep(1, length(x))))
        layout <- c(1, 1, 1)
        number.of.cond <- 1
    }
    if (missing(xlab)) 
        xlab <- form$right.name
    if (missing(ylab)) 
        ylab <- form$left.name
    if (!(is.numeric(x) && is.numeric(y))) 
        warning("x and y are not both numeric")
    foo <- do.call("trellis.skeleton", c(list(aspect = aspect, 
        strip = strip, panel = panel, xlab = xlab, ylab = ylab), 
        dots))
    dots <- foo$dots
    foo <- foo$foo
    foo$call <- match.call()
    foo$fontsize.normal <- 10
    foo$fontsize.small <- 8
    if (is.list(foo$xlab) && !is.characterOrExpression(foo$xlab$label)) 
        foo$xlab$label <- form$right.name
    if (is.list(foo$ylab) && !is.characterOrExpression(foo$ylab$label)) 
        foo$ylab$label <- form$left.name
    if (is.character(scales)) 
        scales <- list(relation = scales)
    foo <- c(foo, do.call("construct.scales", scales))
    have.xlim <- !missing(xlim)
    if (!is.null(foo$x.scales$limit)) {
        have.xlim <- TRUE
        xlim <- foo$x.scales$limit
    }
    have.ylim <- !missing(ylim)
    if (!is.null(foo$y.scales$limit)) {
        have.ylim <- TRUE
        ylim <- foo$y.scales$limit
    }
    have.xlog <- !is.logical(foo$x.scales$log) || foo$x.scales$log
    have.ylog <- !is.logical(foo$y.scales$log) || foo$y.scales$log
    if (have.xlog) {
        xlog <- foo$x.scales$log
        xbase <- if (is.logical(xlog)) 
            10
        else if (is.numeric(xlog)) 
            xlog
        else if (xlog == "e") 
            exp(1)
        x <- log(x, xbase)
        if (have.xlim) 
            xlim <- log(xlim, xbase)
    }
    if (have.ylog) {
        ylog <- foo$y.scales$log
        ybase <- if (is.logical(ylog)) 
            10
        else if (is.numeric(ylog)) 
            ylog
        else if (ylog == "e") 
            exp(1)
        y <- log(y, ybase)
        if (have.ylim) 
            ylim <- log(ylim, ybase)
    }
    cond.max.level <- unlist(lapply(cond, nlevels))
    id.na <- is.na(x) | is.na(y)
    for (var in cond) id.na <- id.na | is.na(var)
    if (!any(!id.na)) 
        stop("nothing to draw")
    foo$condlevels <- lapply(cond, levels)
    foo$panel.args.common <- dots
    if (subscripts) 
        foo$panel.args.common$groups <- groups
    layout <- compute.layout(layout, cond.max.level, skip = foo$skip)
    plots.per.page <- max(layout[1] * layout[2], layout[2])
    number.of.pages <- layout[3]
    foo$skip <- rep(foo$skip, length = plots.per.page)
    foo$layout <- layout
    nplots <- plots.per.page * number.of.pages
    foo$panel.args <- as.list(1:nplots)
    cond.current.level <- rep(1, number.of.cond)
    panel.number <- 1
    for (page.number in 1:number.of.pages) if (!any(cond.max.level - 
        cond.current.level < 0)) 
        for (plot in 1:plots.per.page) {
            if (foo$skip[plot]) 
                foo$panel.args[[panel.number]] <- FALSE
            else if (!any(cond.max.level - cond.current.level < 
                0)) {
                id <- !id.na
                for (i in 1:number.of.cond) {
                  var <- cond[[i]]
                  id <- id & if (is.shingle(var)) 
                    ((var >= levels(var)[[cond.current.level[i]]][1]) & 
                      (var <= levels(var)[[cond.current.level[i]]][2]))
                  else (as.numeric(var) == cond.current.level[i])
                }
                foo$panel.args[[panel.number]] <- list(x = x[id], 
                  y = y[id])
                if (subscripts) 
                  foo$panel.args[[panel.number]]$subscripts <- subscr[id]
                cond.current.level <- cupdate(cond.current.level, 
                  cond.max.level)
            }
            panel.number <- panel.number + 1
        }
    foo <- c(foo, limits.and.aspect(prepanel.default.xyplot, 
        prepanel = prepanel, have.xlim = have.xlim, xlim = xlim, 
        have.ylim = have.ylim, ylim = ylim, x.relation = foo$x.scales$relation, 
        y.relation = foo$y.scales$relation, panel.args.common = foo$panel.args.common, 
        panel.args = foo$panel.args, aspect = aspect, nplots = nplots, 
        x.axs = foo$x.scales$axs, y.axs = foo$y.scales$axs))
    if (is.null(foo$key) && !is.null(groups) && (is.list(auto.key) || 
        (is.logical(auto.key) && auto.key))) 
        foo$key <- do.call("simpleKey", c(list(levels(as.factor(groups))), 
            if (is.list(auto.key)) auto.key else list()))
    class(foo) <- "trellis"
    foo
}
