#' Extract the top up and down regulated MR per sample
#'
#' @param vipres numeric matrix of NES as returned by the viper function
#' @param nn integer of number of top MR to consider for each regulon
#' @param direction character, which tail should be returned, defaults to both up and down
#' @param make_regulon logical, whether to set class of return object to "regulon"
#' @return list object with top dysregulated MR per sample
#' @export

vpres2regulon <- function(vpres, nn = 25, direction = c('both', 'up', 'down'), make_regulon = FALSE) {

    dir <- match.arg(direction)

    switch(dir,
           both = {idx <- c(1:nn, (nrow(vpres)-(nn-1)):nrow(vpres))
           tfm <- rep(c(1,-1), each = nn) },
           up = {idx <- 1:nn
           tfm <- rep(1, nn)},
           down = {idx <- (nrow(vpres)-(nn-1)):nrow(vpres)
           tfm <- rep(-1, nn)}
    )

    tmp <- apply(vpres, 2, function(i) {
        names(tfm) <- rownames(vpres[order(i, decreasing = TRUE), ])[idx]
        list(tfmode = tfm, likelihood = unname(abs(tfm)))
    })

    if(make_regulon) class(tmp) <- "regulon"
    return(tmp)
}



