# ------------------------------------------------------------------------------
# Event Study Estimator for Walmart data 

using Plots, PrettyTables, CSV, DataFrames, GLM, StatsBase, LaTeXStrings

YEARS = collect(1977:1999)
T0 = 1985
B = 1000
ests = [(p = 2, outcome = "retail"), (p = 1, outcome = "wholesale")]


for est in ests
  outcome = est.outcome
  p = est.p
  
  lims = outcome == "retail" ? [-0.175, 0.3] : [-0.4, 0.2]
  outcome_var = "log_$(outcome)_emp"

  # Load bootstrap results

  # First row is point estimate
  factor_est = DataFrame(CSV.File(
    "data/factor_est_$(outcome)_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0)_B_$(B)_p_$(p).csv"
  ))

  factor_est = Matrix(factor_est)

  factor_estimate = factor_est[1, :]
  factor_upper_95 = zeros(size(factor_est)[2])
  factor_lower_95 = zeros(size(factor_est)[2])
  factor_upper_90 = zeros(size(factor_est)[2])
  factor_lower_90 = zeros(size(factor_est)[2])
  factor_se_bootstrap = zeros(size(factor_est)[2])
  for i in 1:size(factor_est)[2]
    col = factor_est[:, i]

    factor_upper_95[i] = quantile(col, 0.975)
    factor_lower_95[i] = quantile(col, 0.025)
    factor_upper_90[i] = quantile(col, 0.95)
    factor_lower_90[i] = quantile(col, 0.05)
    factor_se_bootstrap[i] = std(col)
  end


  did2s_est = DataFrame(CSV.File(
    "data/did2s_est_$(outcome)_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0)_B_$(B).csv"
  ))

  did2s_estimate = did2s_est.estimate
  did2s_se_bootstrap = did2s_est.std_error
  did2s_upper = did2s_est.upper
  did2s_lower = did2s_est.lower






  # Plot TE Estimates ------------------------------------------------------------


  # Extract point estimates 
  rel_year = collect(-22:13)

  idx_pre  = (rel_year .<  0)
  idx_post = (rel_year .>= 0)

  # ATT
  # mean(factor_estimate[idx_post])

  # Linear extrapolation
  pre_ests = DataFrame(
    rel_year   = rel_year[(idx_pre) .& (rel_year .>= -15)], 
    did2s_estimate  = did2s_estimate[(idx_pre) .& (rel_year .>= -15)],
    factor_estimate = factor_estimate[(idx_pre) .& (rel_year .>= -15)]
  )
  did2s_extrapolation = lm(
    @formula(did2s_estimate ~ 1 + rel_year), 
    pre_ests
  ) |> coef
  factor_extrapolation = lm(
    @formula(factor_estimate ~ 1 + rel_year), 
    pre_ests
  ) |> coef



  pgfplotsx()
  # gr(size = (600, 300))

  did2s_plot = hline([0], linestyle = :dash, linecolor = :black)
  # title!("TWFE Model Imputation")
  xlabel!("Event Time")
  ylabel!("Coefficient")
  ylims!(lims[1], lims[2])
  Plots.abline!(
    did2s_extrapolation[2], did2s_extrapolation[1],
    color = RGBA(154/255, 36/255, 21/255, 1)
  )
  scatter!(
    rel_year[idx_pre],
    did2s_estimate[idx_pre], 
    # yerror = 1.96 .* did2s_se_bootstrap[idx_pre],
    yerror = (
      did2s_estimate[idx_pre] - did2s_lower[idx_pre], 
      did2s_upper[idx_pre] - did2s_estimate[idx_pre]
    ),
    legend = false,
    markersize = 2, 
    markerstrokecolor = RGBA(154/255, 36/255, 21/255, 1),
    markercolor = RGBA(154/255, 36/255, 21/255, 1)
  )
  scatter!(
    rel_year[idx_post],
    did2s_estimate[idx_post],
    # yerror = 1.96 .* did2s_se_bootstrap[idx_post],
    yerror = (
      did2s_estimate[idx_post] - did2s_lower[idx_post], 
      did2s_upper[idx_post] - did2s_estimate[idx_post]
    ),
    legend = false,
    markersize = 2, 
    markerstrokecolor = RGBA(16/255, 120/255, 149/255, 1),
    markercolor = RGBA(16/255, 120/255, 149/255, 1)
  )


  factor_plot = plot(size = (600, 400) .* 0.75)
  hline!([0], linestyle = :dash, linecolor = :black)
  # title!("Factor Model Imputation")
  xlabel!("Event Time")
  ylabel!("Coefficient")
  ylims!(lims[1], lims[2])
  Plots.abline!(
    factor_extrapolation[2], factor_extrapolation[1],
    color = RGBA(154/255, 36/255, 21/255, 1)
  )
  scatter!(
    rel_year[idx_pre],
    factor_estimate[idx_pre],
    yerror = 1.96 .* factor_se_bootstrap[idx_pre],
    # yerror = (
    #   factor_estimate[idx_pre] - factor_lower_95[idx_pre], 
    #   factor_upper_95[idx_pre] - factor_estimate[idx_pre]
    # ),
    legend = false,
    markersize = 2, 
    markerstrokecolor = RGBA(154/255, 36/255, 21/255, 1),
    markercolor = RGBA(154/255, 36/255, 21/255, 1)
  )
  scatter!(
    rel_year[idx_post],
    factor_estimate[idx_post],
    yerror = 1.96 .* factor_se_bootstrap[idx_post],
    # yerror = (
    #   factor_estimate[idx_post] - factor_lower_95[idx_post], 
    #   factor_upper_95[idx_post] - factor_estimate[idx_post]
    # ),
    legend = false,
    markersize = 2, 
    markerstrokecolor = RGBA(16/255, 120/255, 149/255, 1),
    markercolor = RGBA(16/255, 120/255, 149/255, 1)
  )

  print("plotting YEARS: $(minimum(YEARS)):$(maximum(YEARS)); T0: $(T0); p: $(p)\n")

  did2s_plot = plot(
    did2s_plot, size = (600, 400) .* 0.65
    # did2s_plot, size = (600, 300) .* 1.25
  )
  factor_plot = plot(
    factor_plot, size = (600, 400) .* 0.65
    # factor_plot, size = (600, 300) .* 1.25
  )

  combined_plot = plot(
    did2s_plot, factor_plot, layout = (1,2), legend = false, size = (1200, 400)
    # did2s_plot, factor_plot, layout = (2, 1), legend = false, size = (600, 800)
  ) 


  # Export to .tex
  savefig(did2s_plot, "figures/did2s_$(outcome)_bootstrap_$(B).tex")
  savefig(factor_plot, "figures/factor_$(outcome)_p_$(p)_bootstrap_$(B).tex")

  # Export to .pdf
  savefig(did2s_plot, "figures/did2s_$(outcome)_bootstrap_$(B).pdf")
  savefig(factor_plot, "figures/factor_$(outcome)_p_$(p)_bootstrap_$(B).pdf")

  # Export for slides .tex
  did2s_plot_slide = plot(did2s_plot)
  ylabel!("")
  savefig(
    plot(did2s_plot_slide, size=(475,300)), 
    "figures/did2s_$(outcome)_bootstrap_$(B)_slides.tex"
  )
  factor_plot_slide = plot(factor_plot)
  ylabel!("")
  savefig(
    plot(factor_plot_slide, size=(475,300)), 
    "figures/factor_$(outcome)_p_$(p)_bootstrap_$(B)_slides.tex"
  )

  # savefig(factor_plot, "figures/factor_$(outcome)_p_$(p)_bootstrap_$(B).tex")


  # --------

end



