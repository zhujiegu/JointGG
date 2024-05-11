library(mvQuad)
nw <- createNIGrid(dim=2, type="GHe", level=4)

nw <- createNIGrid(dim=2, type="GHe", level=9,ndConstruction = "sparse")

# no rescaling
plot(nw, main="initial grid",, pch=8)

rescale(nw, m = mu, C = sigma, dec.type = 1)
plot(nw, main="no dec.", xlim=c(-6,6), ylim=c(-6,6), pch=8)


rescale(nw, m = mu, C = params$G, dec.type = 2)
plot(nw, main="Cholesky dec.", pch=8)
