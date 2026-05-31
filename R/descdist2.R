#' @title Descdist but with colors
#' @description
#' Adds bg/data color to descdist plot
#'
#' @export
descdist2 = function (data, discrete = FALSE, boot = NULL, method = "unbiased",
graph = TRUE, print = TRUE, obs.col = "red", obs.pch = 16,
boot.col = "orange", data.col = "lightgrey")
{
  if (missing(data) || !is.vector(data, mode = "numeric"))
    stop("data must be a numeric vector")
  if (length(data) < 4)
    stop("data must be a numeric vector containing at least four values")
  moment <- function(data, k) {
    m1 <- mean(data)
    return(sum((data - m1)^k)/length(data))
  }
  if (method == "unbiased") {
    skewness <- function(data) {
      sd <- sqrt(moment(data, 2))
      n <- length(data)
      gamma1 <- moment(data, 3)/sd^3
      unbiased.skewness <- sqrt(n * (n - 1)) * gamma1/(n -
                                                         2)
      return(unbiased.skewness)
    }
    kurtosis <- function(data) {
      n <- length(data)
      var <- moment(data, 2)
      gamma2 <- moment(data, 4)/var^2
      unbiased.kurtosis <- (n - 1)/((n - 2) * (n - 3)) *
        ((n + 1) * gamma2 - 3 * (n - 1)) + 3
      return(unbiased.kurtosis)
    }
    standdev <- function(data) {
      sd(data)
    }
  }
  else if (method == "sample") {
    skewness <- function(data) {
      sd <- sqrt(moment(data, 2))
      return(moment(data, 3)/sd^3)
    }
    kurtosis <- function(data) {
      var <- moment(data, 2)
      return(moment(data, 4)/var^2)
    }
    standdev <- function(data) {
      sqrt(moment(data, 2))
    }
  }
  else stop("The only possible value for the argument method are 'unbiased' or 'sample'")
  res <- list(min = min(data), max = max(data), median = median(data),
              mean = mean(data), sd = standdev(data), skewness = skewness(data),
              kurtosis = kurtosis(data), method = method)
  skewdata <- res$skewness
  kurtdata <- res$kurtosis
  if (graph) {
    if (!is.null(boot)) {
      if (!is.numeric(boot) || boot < 10) {
        stop("boot must be NULL or a integer above 10")
      }
      n <- length(data)
      databoot <- matrix(sample(data, size = n * boot,
                                replace = TRUE), nrow = n, ncol = boot)
      s2boot <- sapply(1:boot, function(iter) skewness(databoot[,
                                                                iter])^2)
      kurtboot <- sapply(1:boot, function(iter) kurtosis(databoot[,
                                                                  iter]))
      kurtmax <- max(10, ceiling(max(kurtboot)))
      xmax <- max(4, ceiling(max(s2boot)))
    }
    else {
      kurtmax <- max(10, ceiling(kurtdata))
      xmax <- max(4, ceiling(skewdata^2))
    }
    ymax <- kurtmax - 1
    plot(skewdata^2, kurtmax - kurtdata, pch = "",
         xlim = c(0, xmax), ylim = c(0, ymax), yaxt = "n",
         xlab = "square of skewness", ylab = "kurtosis",
         main = "Cullen and Frey graph")
    yax <- as.character(kurtmax - 0:ymax)
    axis(side = 2, at = 0:ymax, labels = yax)
    if (!discrete) {
      p <- exp(-100)
      lq <- seq(-100, 100, 0.1)
      q <- exp(lq)
      s2a <- (4 * (q - p)^2 * (p + q + 1))/((p + q + 2)^2 *
                                              p * q)
      ya <- kurtmax - (3 * (p + q + 1) * (p * q * (p +
                                                     q - 6) + 2 * (p + q)^2)/(p * q * (p + q + 2) *
                                                                                (p + q + 3)))
      p <- exp(100)
      lq <- seq(-100, 100, 0.1)
      q <- exp(lq)
      s2b <- (4 * (q - p)^2 * (p + q + 1))/((p + q + 2)^2 *
                                              p * q)
      yb <- kurtmax - (3 * (p + q + 1) * (p * q * (p +
                                                     q - 6) + 2 * (p + q)^2)/(p * q * (p + q + 2) *
                                                                                (p + q + 3)))
      s2 <- c(s2a, s2b)
      y <- c(ya, yb)
      polygon(s2, y, col = data.col, border = data.col)
      lshape <- seq(-100, 100, 0.1)
      shape <- exp(lshape)
      s2 <- 4/shape
      y <- kurtmax - (3 + 6/shape)
      lines(s2[s2 <= xmax], y[s2 <= xmax], lty = 2)
      lshape <- seq(-100, 100, 0.1)
      shape <- exp(lshape)
      es2 <- exp(shape^2)
      s2 <- (es2 + 2)^2 * (es2 - 1)
      y <- kurtmax - (es2^4 + 2 * es2^3 + 3 * es2^2 -
                        3)
      lines(s2[s2 <= xmax], y[s2 <= xmax], lty = 3)
      legend2(xmax * 0.7, ymax * 1.15, legend = "Theoretical", # 1.1
              bty = "n", cex = 0.8, text.font = 2)
      legend2(xmax * 0.75, ymax * 1.09, pch = 8, legend = "normal",
              bty = "n", cex = 0.8)
      legend2(xmax * 0.75, ymax * 1.05, pch = 2, legend = "uniform",
              bty = "n", cex = 0.8)
      legend2(xmax * 0.75, ymax * 1.01, pch = 7, legend = "exponential",
              bty = "n", cex = 0.8)
      legend2(xmax * 0.75, ymax * 0.97, pch = 3, legend = "logistic",
              bty = "n", cex = 0.8)
      legend2(xmax * 0.75, ymax * 0.93, fill = data.col,
              legend = "beta", bty = "n", cex = 0.8)
      legend2(xmax * 0.75, ymax * 0.89, lty = 3, legend = "lognormal",
              bty = "n", cex = 0.8)
      legend2(xmax * 0.75, ymax * 0.85, lty = 2, legend = "gamma",
              bty = "n", cex = 0.8)
      legend2(xmax * 0.7, ymax * 0.845, legend = c("(Weibull is close to \n gamma and lognormal)"),
              bty = "n", cex = 0.6)
      legend2(xmax * 0.7, ymax * 0.74, legend = "Empirical",
              bty = "n", cex = 0.8, text.font = 2)
      legend2(xmax * 0.75, ymax * 0.70, pch = obs.pch, legend = "data",
              bty = "n", cex = 0.8, pt.cex = 1.2, col = obs.col)
      if (!is.null(boot)) {
        legend(xmax * 0.75, ymax * 0.67, pch = 1, legend = "bootstrap",
               bty = "n", cex = 0.8, col = boot.col)
      }
    }
    else {
      p <- exp(-10)
      lr <- seq(-100, 100, 0.1)
      r <- exp(lr)
      s2a <- (2 - p)^2/(r * (1 - p))
      ya <- kurtmax - (3 + 6/r + p^2/(r * (1 - p)))
      p <- 1 - exp(-10)
      lr <- seq(100, -100, -0.1)
      r <- exp(lr)
      s2b <- (2 - p)^2/(r * (1 - p))
      yb <- kurtmax - (3 + 6/r + p^2/(r * (1 - p)))
      s2 <- c(s2a, s2b)
      y <- c(ya, yb)
      polygon(s2, y, col = data.col, border = data.col)
      legend2(xmax * 0.73, ymax * 1.15, legend = "Theoretical",
              bty = "n", cex = 0.8, text.font = 2)
      legend2(xmax * 0.78, ymax * 1.09, pch = 8, legend = "normal",
              bty = "n", cex = 0.8)
      legend2(xmax * 0.78, ymax * 1.135, fill = data.col,
              legend = "negative \nbinomial", bty = "n", cex = 0.8)
      legend2(xmax * 0.78, ymax * 0.98, lty = 2, legend = "Poisson",
              bty = "n", cex = 0.8)
      legend2(xmax * 0.73, ymax * 0.92, legend = "Empirical",
              bty = "n", cex = 0.8, text.font = 2)
      legend2(xmax * 0.78, ymax * 0.88, pch = obs.pch, legend = "data",
              bty = "n", cex = 0.8, pt.cex = 1.2, col = obs.col)
      if (!is.null(boot)) {
        legend2(xmax * 0.78, ymax * 0.87, pch = 1, legend = "bootstrap",
                bty = "n", cex = 0.8, col = boot.col)
      }
      llambda <- seq(-100, 100, 0.1)
      lambda <- exp(llambda)
      s2 <- 1/lambda
      y <- kurtmax - (3 + 1/lambda)
      lines(s2[s2 <= xmax], y[s2 <= xmax], lty = 2)
    }
    if (!is.null(boot)) {
      points(s2boot, kurtmax - kurtboot, pch = 1, col = boot.col,
             cex = 0.5)
    }
    points(skewness(data)^2, kurtmax - kurtosis(data), pch = obs.pch,
           cex = 2, col = obs.col)
    points(0, kurtmax - 3, pch = 8, cex = 1.5, lwd = 2)
    if (!discrete) {
      points(0, kurtmax - 9/5, pch = 2, cex = 1.5, lwd = 2)
      points(2^2, kurtmax - 9, pch = 7, cex = 1.5, lwd = 2)
      points(0, kurtmax - 4.2, pch = 3, cex = 1.5, lwd = 2)
    }
  }
  if (!print)
    invisible(structure(res, class = "descdist"))
  else structure(res, class = "descdist")
}

# legend ------------------------------------------------------------------

legend2 = function (x, y = NULL, legend, fill = NULL, col = par("col"),
                    border = "black", lty, lwd, pch, angle = 45, density = NULL,
                    bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"),
                    box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex,
                    pt.lwd = lwd, xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
                    adj = c(0, 0.5), text.width = NULL, text.col = par("col"),
                    text.font = NULL, merge = do.lines && has.pch, trace = FALSE,
                    plot = TRUE, ncol = 1, horiz = FALSE, title = NULL, inset = 0,
                    xpd, title.col = text.col[1], title.adj = 0.5, title.cex = cex[1],
                    title.font = text.font[1], seg.len = 2)
{
  if (missing(legend) && !missing(y) && (is.character(y) ||
                                         is.expression(y))) {
    legend <- y
    y <- NULL
  }
  mfill <- !missing(fill) || !missing(density)
  if (!missing(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd = op))
    par(xpd = xpd)
  }
  if (is.null(text.font))
    text.font <- par("font")
  title <- as.graphicsAnnot(title)
  if (length(title) > 1)
    stop("invalid 'title'")
  legend <- as.graphicsAnnot(legend)
  if (any(vapply(legend, is.language, NA)))
    legend <- as.expression(legend)
  n.leg <- length(legend)
  if (n.leg == 0)
    stop("'legend' is of length 0")
  auto <- if (is.character(x))
    match.arg(x, c("bottomright", "bottom", "bottomleft",
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y, setLab = FALSE)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2)
      stop("invalid coordinate lengths")
  }
  else nx <- 0
  reverse.xaxis <- par("xaxp")[1] > par("xaxp")[2]
  reverse.yaxis <- par("yaxp")[1] > par("yaxp")[2]
  xlog <- par("xlog")
  ylog <- par("ylog")
  cex <- rep(cex, length.out = n.leg)
  x.intersp <- rep(x.intersp, length.out = n.leg)
  seg.len <- rep(seg.len, length.out = n.leg)
  rect2 <- function(left, top, dx, dy, density = NULL, angle,
                    ...) {
    r <- left + dx
    if (xlog) {
      left <- 10^left
      r <- 10^r
    }
    b <- top - dy
    if (ylog) {
      top <- 10^top
      b <- 10^b
    }
    rect(left, top, r, b, angle = angle, density = density,
         ...)
  }
  segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx
    if (xlog) {
      x1 <- 10^x1
      x2 <- 10^x2
    }
    y2 <- y1 + dy
    if (ylog) {
      y1 <- 10^y1
      y2 <- 10^y2
    }
    segments(x1, y1, x2, y2, ...)
  }
  points2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    points(x, y, ...)
  }
  text2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    text(x, y, ...)
  }
  colwise <- function(x, n, ncol, n.legpercol, fun, reverse = FALSE) {
    xmat <- matrix(c(rep(x, length.out = n), rep(0L, n.legpercol *
                                                   ncol - n)), ncol = ncol)
    res <- apply(xmat, 2, fun)
    res[res == 0L] <- max(res)
    if (reverse)
      -res
    else res
  }
  rowwise <- function(x, n, ncol, n.legpercol, fun, reverse = FALSE) {
    xmat <- matrix(c(rep(x, length.out = n), rep(0L, n.legpercol *
                                                   ncol - n)), ncol = ncol)
    res <- apply(xmat, 1, fun)
    if (reverse)
      -res
    else res
  }
  if (trace) {
    catn <- function(...) do.call(cat, c(lapply(list(...),
                                                formatC), "\n"))
    fv <- function(...) paste(vapply(lapply(list(...), formatC),
                                     paste, collapse = ",", ""), collapse = ", ")
  }
  n.legpercol <- if (horiz) {
    if (ncol != 1)
      warning(gettextf("horizontal specification overrides: Number of columns := %d",
                       n.leg), domain = NA)
    ncol <- n.leg
    1
  }
  else ceiling(n.leg/ncol)
  Cex <- cex * par("cex")
  if (is.null(text.width))
    text.width <- max(abs(mapply(strwidth, legend, cex = cex,
                                 font = text.font, MoreArgs = list(units = "user"))))
  else if ((length(text.width) > 1L && any(is.na(text.width))) ||
           (all(!is.na(text.width)) && (!is.numeric(text.width) ||
                                        any(text.width < 0))))
    stop("'text.width' must be numeric, >= 0, or a scalar NA")
  if (auto.text.width <- all(is.na(text.width))) {
    text.width <- abs(mapply(strwidth, legend, cex = cex,
                             font = text.font, MoreArgs = list(units = "user")))
    ncol <- ceiling(n.leg/n.legpercol)
  }
  xyc <- xyinch(par("cin"), warn.log = FALSE)
  xc <- Cex * xyc[1L]
  yc <- Cex * xyc[2L]
  if (any(xc < 0))
    text.width <- -text.width
  xchar <- xc
  xextra <- 0
  y.intersp <- rep(y.intersp, length.out = n.legpercol)
  yextra <- rowwise(yc, n = n.leg, ncol = ncol, n.legpercol = n.legpercol,
                    fun = function(x) max(abs(x)), reverse = reverse.yaxis) *
    (y.intersp - 1)
  ymax <- sign(yc[1]) * max(abs(yc)) * max(1, mapply(strheight,
                                                     legend, cex = cex, font = text.font, MoreArgs = list(units = "user"))/yc)
  ychar <- yextra + ymax
  ymaxtitle <- title.cex * par("cex") * xyc[2L] * max(1, strheight(title,
                                                                   cex = title.cex, font = title.font, units = "user")/(title.cex *
                                                                                                                          par("cex") * xyc[2L]))
  ychartitle <- yextra[1] + ymaxtitle
  if (trace)
    catn("  xchar=", fv(xchar), "; (yextra, ychar)=", fv(yextra,
                                                         ychar))
  if (mfill) {
    xbox <- xc * 0.25
    ybox <- yc * 0.2
    dx.fill <- max(xbox)
  }
  do.lines <- (!missing(lty) && (is.character(lty) || any(lty >
                                                            0))) || !missing(lwd)
  has.pch <- !missing(pch) && length(pch) > 0
  if (do.lines) {
    x.off <- if (merge)
      -0.7
    else 0
  }
  else if (merge)
    warning("'merge = TRUE' has no effect when no line segments are drawn")
  if (has.pch) {
    if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L],
                                                      type = "c") > 1) {
      if (length(pch) > 1)
        warning("not using pch[2..] since pch[1L] has multiple chars")
      np <- nchar(pch[1L], type = "c")
      pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
    }
    if (!is.character(pch))
      pch <- as.integer(pch)
  }
  if (is.na(auto)) {
    if (xlog)
      x <- log10(x)
    if (ylog)
      y <- log10(y)
  }
  if (nx == 2) {
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top <- y[2L]
    w <- diff(x)
    h <- diff(y)
    w0 <- w/ncol
    x <- mean(x)
    y <- mean(y)
    if (missing(xjust))
      xjust <- 0.5
    if (missing(yjust))
      yjust <- 0.5
  }
  else {
    yc <- rowwise(yc, n.leg, ncol, n.legpercol, fun = function(x) max(abs(x)),
                  reverse = reverse.yaxis)
    h <- sum(ychar) + yc[length(yc)] + (!is.null(title)) *
      ychartitle
    xch1 <- colwise(xchar, n.leg, ncol, n.legpercol, fun = function(x) max(abs(x)),
                    reverse = reverse.xaxis)
    x.interspCol <- colwise(x.intersp, n.leg, ncol, n.legpercol,
                            fun = max)
    seg.lenCol <- colwise(seg.len, n.leg, ncol, n.legpercol,
                          fun = max)
    text.width <- colwise(text.width, n = if (auto.text.width)
      n.leg
      else ncol, ncol, n.legpercol = if (auto.text.width)
        n.legpercol
      else 1, fun = function(x) max(abs(x)), reverse = reverse.xaxis)
    w0 <- text.width + (x.interspCol + 1) * xch1
    if (mfill)
      w0 <- w0 + dx.fill
    if (do.lines)
      w0 <- w0 + (seg.lenCol + x.off) * xch1
    w <- sum(w0) + 0.5 * xch1[ncol]
    if (!is.null(title) && (abs(tw <- strwidth(title, units = "user",
                                               cex = title.cex, font = title.font) + 0.5 * title.cex *
                                par("cex") * xyc[1L])) > abs(w)) {
      xextra <- (tw - w)/2
      w <- tw
    }
    if (is.na(auto)) {
      left <- x - xjust * w
      top <- y + (1 - yjust) * h
    }
    else {
      usr <- par("usr")
      inset <- rep_len(inset, 2)
      insetx <- inset[1L] * (usr[2L] - usr[1L])
      left <- switch(auto, bottomright = , topright = ,
                     right = usr[2L] - w - insetx, bottomleft = ,
                     left = , topleft = usr[1L] + insetx, bottom = ,
                     top = , center = (usr[1L] + usr[2L] - w)/2)
      insety <- inset[2L] * (usr[4L] - usr[3L])
      top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] +
                      h + insety, topleft = , top = , topright = usr[4L] -
                      insety, left = , right = , center = (usr[3L] +
                                                             usr[4L] + h)/2)
    }
  }
  if (plot && bty != "n") {
    if (trace)
      catn("  rect2(", left, ",", top, ", w=", w, ", h=",
           h, ", ...)", sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL,
          lwd = box.lwd, lty = box.lty, border = box.col)
  }
  xt <- left + xc + xextra + rep(c(0, cumsum(w0))[1L:ncol],
                                 each = n.legpercol, length.out = n.leg)
  topspace <- 0.5 * ymax + (!is.null(title)) * ychartitle
  yt <- top - topspace - cumsum((c(0, ychar)/2 + c(ychar,
                                                   0)/2)[1L:n.legpercol])
  yt <- rep(yt, length.out = n.leg)
  if (mfill) {
    if (plot) {
      if (!is.null(fill))
        fill <- rep_len(fill, n.leg)
      rect2(left = xt * 0.982, top = yt + ybox/2, dx = xbox, dy = ybox,
            col = fill, density = density, angle = angle,
            border = border)
    }
    xt <- xt + dx.fill * 0.02
  }
  if (plot && (has.pch || do.lines))
    col <- rep_len(col, n.leg)
  if (missing(lwd) || is.null(lwd))
    lwd <- par("lwd")
  if (do.lines) {
    if (missing(lty) || is.null(lty))
      lty <- 1
    lty <- rep_len(lty, n.leg)
    lwd <- rep_len(lwd, n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0) &
      !is.na(lwd)
    if (trace)
      catn("  segments2(", xt[ok.l] + x.off * xchar[ok.l],
           ",", yt[ok.l], ", dx=", (seg.len * xchar)[ok.l],
           ", dy=0, ...)")
    if (plot)
      segments2((xt[ok.l] + x.off * xchar[ok.l]) * 0.975, yt[ok.l],
                dx = ((seg.len * xchar) * 0.2)[ok.l], dy = 0, lty = lty[ok.l],
                lwd = lwd[ok.l], col = col[ok.l])
    xt <- xt + (seg.len + x.off) * xchar * 0.002
  }
  if (has.pch) {
    pch <- rep_len(pch, n.leg)
    pt.bg <- rep_len(pt.bg, n.leg)
    pt.cex <- rep_len(pt.cex, n.leg)
    pt.lwd <- rep_len(pt.lwd, n.leg)
    ok <- !is.na(pch)
    if (!is.character(pch)) {
      ok <- ok & (pch >= 0 | pch <= -32)
    }
    else {
      ok <- ok & nzchar(pch)
    }
    x1 <- (if (merge && do.lines)
      xt - (seg.len/2) * xchar
      else xt)[ok]
    y1 <- yt[ok]
    if (trace)
      catn("  points2(", x1, ",", y1, ", pch=", pch[ok],
           ", ...)")
    if (plot)
      points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok],
              bg = pt.bg[ok], lwd = pt.lwd[ok])
  }
  xt <- xt + x.intersp * xc * 0.4
  if (plot) {
    if (!is.null(title))
      text2(left + w * title.adj, top - ymaxtitle, labels = title,
            adj = c(title.adj, 0), cex = title.cex, col = title.col,
            font = title.font)
    text2(xt, yt, labels = legend, adj = adj, cex = cex,
          col = text.col, font = text.font)
  }
  invisible(list(rect = list(w = w, h = h, left = left, top = top),
                 text = list(x = xt, y = yt)))
}
