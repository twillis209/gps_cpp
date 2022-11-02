library(argparse)
library(evd)
library(fitdistrplus)

parser <- ArgumentParser(description = 'Computes p-value for GPS statistic using permuted data')
parser$add_argument('-g', '--gps', type = 'double', help = 'GPS statistic computed using non-permuted data', required = T)
parser$add_argument('-p', '--perm_file', type = 'character', help = 'Path to file containing realisations of the GPS value generated from permuted data', required = T)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

perms <- scan(args$perm_file)

fgev.fit <- tryCatch(
   fgev(perms),
  error = function(c) {
    msg <- conditionMessage(c)
    if(msg == "observed information matrix is singular; use std.err = FALSE"){
      fgev(perms, std.err = F)
    } else {
      stop(msg)
      }
    }
)

fgev.fitdist <- fitdist(perms, 'gev', start = list(loc = fgev.fit$estimate[['loc']], scale = fgev.fit$estimate[['scale']], shape = fgev.fit$estimate[['shape']]))

loc <- fgev.fitdist$estimate[['loc']]
loc.sd <- fgev.fitdist$sd[['loc']]
scale <- fgev.fitdist$estimate[['scale']]
scale.sd <- fgev.fitdist$sd[['scale']]
shape <- fgev.fitdist$estimate[['shape']]
shape.sd <- fgev.fitdist$sd[['shape']]

pvalue <- pgev(args$gps, loc = loc, scale = scale, shape = shape, lower.tail = F)

write.table(data.frame(gps = args$gps, n = length(perms), loc = loc, loc.sd = loc.sd, scale = scale, scale.sd = scale.sd, shape = shape, shape.sd = shape.sd, pval = pvalue), sep = '\t', file = args$output_file, quote = F, row.names = F)
