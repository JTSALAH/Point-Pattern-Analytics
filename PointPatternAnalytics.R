# Point Pattern Analytics

# 0. Load in Packages
require(spatstat)
require(spatstat.data)
require(spdep)
require(sf)
require(here)

# 1. Prepare Data
# 1.1: Load in Data
pts <- read.csv(here("data", "cactus.csv"))
boundary <- read.csv(here("data", "cactus_boundaries.csv"),header=T)

# 1.2: Create a spatstat object  with our pts data
ppp.window <- owin(xrange=c(boundary$Xmin, boundary$Xmax),
                   yrange=c(boundary$Ymin, boundary$Ymax))
ppp <- ppp(pts$East, pts$North, window=ppp.window)


# 2. Plot raw data and density
par(mfrow = c(1,2), oma=c(0,0,0,1))
plot(ppp, main = "Points")
plot(density(ppp,1), main = "Density")  #the 1 alters the smoothing parameter

# 3. Take a look at the data's ppp summary
summary(ppp)

# 4. Create plotting template
ppp_plot = function(fun_name, none, iso, trans) {
  par(mfrow = c(1,4))
  plot(none, main = paste(fun_name, "none"),legend=F)         
  plot(none, . - r~r, main = paste(fun_name, "none"), legend=F)  
  plot(iso, . - r~r, main = paste(fun_name, "iso"), legend=F)
  plot(trans, . - r~r, main = paste(fun_name, "trans"), legend=F)
}

# 5. Ripleys K with various corrections
# 5.0: All lines
Kall <- Kest(ppp)

# 5.1: 1:1 expectation (no correction)
Knone <- Kest(ppp, correction="none")

# 5.2: Isotropic edge correction
Kiso <- Kest(ppp, correction="isotropic")

# 5.3: Translate (toroidal) edge correction
Ktrans <- Kest(ppp, correction="trans")

# 5.4: Plot!
ppp_plot("K", Knone, Kiso, Ktrans)

# 6. L with various corrections
# 6.1: 1:1 expectation (no correction)
Lnone <- Lest(ppp, correction="none")

# 6.2: Isotropic edge correction
Liso <- Lest(ppp, correction="isotropic")

# 6.3: Translate (toroidal) edge correction
Ltrans <- Lest(ppp, correction="trans")

# 6.4: Plot!
ppp_plot("L", Lnone, Liso, Ltrans)

# 7. Pair correlation function, g, with various corrections
# 7.1: 1:1 expectation (no correction)
Pnone <- pcf(ppp, correction="none")

# 7.2: Isotropic edge correction
Piso <- pcf(ppp, correction="isotropic")

# 7.3: Translate (toroidal) edge correction
Ptrans <- pcf(ppp, correction="trans")

# 7.4: Plot!
par(mfrow = c(1,3))
plot(Pnone, main = "Pnone",legend=F, ylim=c(0,3))         
plot(Piso, main = "Piso", legend=F, ylim=c(0,3))
plot(Ptrans, main = "Ptrans", legend=F, ylim=c(0,3))

# 8. Basic G function
# Note: G & F don't have an isometric or trans correction, but they have similar corrections.
# 8.1: 1:1 expectation (no correction)
Gnone <- Gest(ppp, correction="none")

# 8.2: Reduced sample or border correction
Grs <- Gest(ppp, correction="rs")

# 8.3: Best (determines best correction for dataset)
Gbest = Gest(ppp, correction="best")

# 8.4: Plot!
par(mfrow = c(1,3))
plot(Gnone, main = "Gnone",legend=F)         
plot(Grs, main = "Grs", legend=F)
plot(Gbest, main = "Gbest", legend=F)

# 9. Basic F function
# 9.1: 1:1 expectation (no correction)
Fnone <- Fest(ppp, correction="none")

# 9.2: Reduced sample or border correction
Frs <- Fest(ppp, correction="rs")

# 9.3: Best (determines best correction for dataset)
Fbest = Fest(ppp, correction="best")

# 9.4: Plot!
par(mfrow = c(1,3))
plot(Fnone, main = "Fnone",legend=F)         
plot(Frs, main = "Frs", legend=F)
plot(Fbest, main = "Fbest", legend=F)

# ---- Point Pattern Process Envelopes  ----

# 1. Create a Lest simulated envelope of global and pointwise confidence under CSR
# 1.1: Create a global & pointwise (non-global) Envelope
Lcsr   <- envelope(ppp, Lest, nsim=99, rank=1, correction="trans", global=F)
Lcsr.g <- envelope(ppp, Lest, nsim=99, rank=1, correction="trans", global=T)

# 1.2: Plot point-wise envelope
plot(Lcsr, . - r~r, shade=c("hi", "lo"), legend=F)

# 1.3: Plot global envelope
plot(Lcsr.g, . - r~r, shade=c("hi", "lo"), legend=F)

# 2. Create a pcf simulated envelope of pointwise confidence under CSR
# 2.1: Create a pair correlation function, g, with trans correction
Ptrans <- pcf(ppp, correction="trans")

# 2.2: Create a fine envelope
Penv <- envelope(ppp,pcf, nsim=99, rank=1, stoyan=0.15, correction="trans", global=F)#stoyan = bandwidth; set to default

# 2.3: Create a coarse envelope
Penv.coarse <- envelope(ppp, pcf, nsim=99, rank=1, stoyan=0.3, correction="trans", global=F)

# 2.4: Plot!
# 2.4.1: Plot no-envelope Ptrans
plot(Ptrans, legend=FALSE, ylim = c(0,3))

# 2.4.2: Plot our fine envelope
plot(Penv, shade=c("hi", "lo"), legend=FALSE, ylim = c(0,3))

# 2.4.3: Plot our coarse envelope
plot(Penv.coarse, shade=c("hi", "lo"), legend=F, ylim = c(0,3))

# 3. Create a Gest simulated envelope of pointwise confidence under CSR
# 3.1: Create a G estimation with trans correction
Gtrans <- Gest(ppp, correction="rs")

# 3.2: Create a pointwise Gest envelope
Genv <- envelope(ppp, Gest, nsim=99, rank=1, correction="rs", global=F)

# 3.3: Create a nearest neighbor distance variable for our plot
nn.dist <- nndist(ppp)
max(nn.dist)

# 3.4: Plot our trans G
plot(Gtrans, legend=F)

# 3.5: Plot G with our pointwise envelope & nearest neighbor distances
plot(Genv, shade=c("hi", "lo"), legend=F)
plot(ecdf(nn.dist), add=T)

# 4. Mark-Correlation Analysis
data(spruces)

# 4.1: Create an envelope for spruces
MCFenv <- envelope(spruces, markcorr, nsim=99, correction="iso", global=F)

# 4.2: Plot envelope
plot(MCFenv,  shade=c("hi", "lo"), legend=F)
