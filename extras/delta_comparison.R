conf = 0.1
c = 2/3
lambda = 2/3
alpha = 0.01
R_c_lam = (1 - lambda)/(1 - lambda + c)
m_conf = ceiling(-log(1/conf, base = 1 - R_c_lam) - 1e-13)
i_0 = 10000

deltas_RESET = rep(0, i_0 + 1)

for (j in i_0:(2*i_0)) {
  d = 0
  while (stats::pbinom((d + 1), floor((j - (d + 1)) * alpha + 1e-13) + 1 + (d + 1), R_c_lam) <= conf) {
    d = d + 1
  }
  deltas_RESET[j - i_0 + 1] = d
}

################################################################

conf = 0.1
c = 1/2
lambda = 1/2
alpha = 0.01
R_c_lam = (1 - lambda)/(1 - lambda + c)
m_conf = ceiling(-log(1/conf, base = 1 - R_c_lam) - 1e-13)
i_0 = 10000

deltas_FDP_SD = rep(0, i_0 + 1)

for (j in i_0:(2*i_0)) {
  d = 0
  while (stats::pbinom((d + 1), floor((j - (d + 1)) * alpha + 1e-13) + 1 + (d + 1), R_c_lam) <= conf) {
    d = d + 1
  }
  deltas_FDP_SD[j - i_0 + 1] = d
}

df = data.frame(x = rep(i_0:(2*i_0), 2), delta = c(2*deltas_RESET, deltas_FDP_SD), Method = c(rep('RESET', i_0 + 1), rep('FDP-SD', i_0 + 1)))
df$ratio = (2*deltas_RESET)/(deltas_FDP_SD)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

p1 = ggplot(df) + geom_line(aes(x = x, y = delta, group = Method, colour = Method)) + geom_line(aes(x = x, y = ratio*(300), colour = 'RESET:FDP-SD')) +
  scale_y_continuous(expression(delta), limits = c(0, 300), sec.axis = sec_axis(~.*(1/300), name = "Ratio of RESET/FDP-SD")) + xlab('Index')+ scale_colour_manual(values = c(gg_color_hue(2)[1], gg_color_hue(2)[2], "black"))
####################################################
#Doing it again but with i_0 smaller
####################################################


c = 2/3
lambda = 2/3
alpha = 0.01
R_c_lam = (1 - lambda)/(1 - lambda + c)
m_conf = ceiling(-log(1/conf, base = 1 - R_c_lam) - 1e-13)
i_0 = 1000
conf = 0.1

deltas_RESET = rep(0, i_0 + 1)

for (j in i_0:(2*i_0)) {
  d = 0
  while (stats::pbinom((d + 1), floor((j - (d + 1)) * alpha + 1e-13) + 1 + (d + 1), R_c_lam) <= conf) {
    d = d + 1
  }
  deltas_RESET[j - i_0 + 1] = d
}

################################################################

c = 1/2
lambda = 1/2
alpha = 0.01
R_c_lam = (1 - lambda)/(1 - lambda + c)
m_conf = ceiling(-log(1/conf, base = 1 - R_c_lam) - 1e-13)
i_0 = 1000
conf = 0.1

deltas_FDP_SD = rep(0, i_0 + 1)

for (j in i_0:(2*i_0)) {
  d = 0
  while (stats::pbinom((d + 1), floor((j - (d + 1)) * alpha + 1e-13) + 1 + (d + 1), R_c_lam) <= conf) {
    d = d + 1
  }
  deltas_FDP_SD[j - i_0 + 1] = d
}

df = data.frame(x = rep(i_0:(2*i_0), 2), delta = c(2*deltas_RESET, deltas_FDP_SD), Method = c(rep('RESET', i_0 + 1), rep('FDP-SD', i_0 + 1)))
df$ratio = (2*deltas_RESET)/(deltas_FDP_SD) #a single at the start

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ggplot(df) + geom_line(aes(x = x, y = delta, group = Method, colour = Method)) + geom_line(aes(x = x, y = ratio*(30), colour = 'RESET:FDP-SD')) +
  scale_y_continuous(expression(delta), limits = c(0, 30), sec.axis = sec_axis(~.*(1/30), name = "Ratio of RESET/FDP-SD")) + xlab('Index') + scale_colour_manual(values = c(gg_color_hue(2)[1], gg_color_hue(2)[2], "black"))