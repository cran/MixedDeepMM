# Model scope

For entirely continuous data, the original continuous DGMM, RDMM and DStMM
engines are used. For mixed responses, a conditional-measurement single-head or
two-head MCEM estimator combines Gaussian continuous measurements with link
models for binary, ordinal, nominal and count responses. The two-head model
has component-specific continuous-observation loadings, intercepts and residual
variances plus a shared shallow Gaussian-mixture tail.

Mixed Student-t and skew-t distributions extend this parameterisation. The
mixed model does not replicate all the layers, pruning, or numerical estimation
steps of the original published MDGMM procedure. Reported BIC values use a
represented-coefficient count, and Monte Carlo evaluation adds simulation
variability. Investigate parameter recovery, boundary cases, and sensitivity
to MC integration in scientific use.
