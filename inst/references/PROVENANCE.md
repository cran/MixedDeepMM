# Software origins and attribution

MixedDeepMM includes continuous Gaussian, Student-t and skew-t deep mixture
engines alongside project-developed continuous/mixed-response interfaces,
Monte Carlo estimation, diagnostics, and cluster assessment functions.

The `legacy-*.R` modules originate from the research integration archive
`DStMM-main.zip`. Their deep Gaussian model ancestry includes the open-source
`deepgmm` package (Cinzia Viroli and Geoffrey J. McLachlan; CRAN maintainer
Suren Rathnayake), publicly distributed under **GPL (>= 3)**:
https://CRAN.R-project.org/package=deepgmm.

The full MixedDeepMM distribution adopts **GPL (>= 3)**. This covers the
package's distribution terms while preserving the original authors'
attribution. Derived source components may be modified and redistributed
under these terms provided their original licensing conditions are followed.

Other files (`mixed-*.R`, response-type infrastructure, related utilities)
were developed and integrated for the MixedDeepMM research project.
For the statistical scope and differences relative to published M1DGMM and
M2DGMM estimation algorithms, see `inst/MODEL_SCOPE.md`.

Reference: Viroli, C. and McLachlan, G. J. (2019), Deep Gaussian Mixture
Models, Statistics and Computing 29:43-51,
https://doi.org/10.1007/s11222-017-9793-z.
