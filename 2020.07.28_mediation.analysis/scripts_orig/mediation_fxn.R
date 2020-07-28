"
Mediation analysis and plotting

Runs mediation analysis of 1 mediator on 1 independent and dependent variable

REQUIRED
  dat = dataframe containing all variables below
  iv = character string of independent variable name
  dv = character string of dependent variable name
  mediator = character string of mediator variable name

OPTIONAL
  plot = logical if should output a flowchart plot of the mediation results.
         Default is FALSE

EXAMPLE
  mediation.fxn(dat = dat,
                iv = 'FeNO.PreBro', 
                dv = 'dACC', 
                mediator = 'module_P337.1_EOS_02',
                plot = TRUE)
"
  
mediation.fxn <- function(dat, iv, dv, mediator, plot=FALSE){
  require(tidyverse, quietly = TRUE)
  require(mediation, quietly = TRUE)
  require(speedglm, quietly = TRUE)
  options(warn=-1)
  
  ##### Filter to data with all variables #####
  dat.sub <- dat %>% 
    dplyr::select(1, all_of(c(iv, dv, mediator))) %>% 
    drop_na()

  ##### Fit linear models #####
  # Total effect of independent var on dependent var
  fit.totEffect <- lm(reformulate(termlabels=iv, response=dv), 
                      data=dat.sub)
  summ.totEffect <- summary(fit.totEffect)
  # Effect of independent var on mediator
  #Must be signif for mediation analysis
  fit.mediator <- lm(reformulate(termlabels=iv, response=mediator),
                     data=dat.sub)
  summ.mediator <- summary(fit.mediator)
  # Effect mediator on dependent var
  fit.dv <- lm(reformulate(termlabels=c(iv,mediator), response=dv),
               data=dat.sub)
  summ.dv <- summary(fit.dv)
  
  ##### Mediation analysis #####
  mediation <- mediate(model.m = fit.mediator,
                     model.y = fit.dv, 
                     treat = iv, 
                     mediator = mediator,
                     boot = TRUE,
                     use_speed = TRUE)
  mediation.summ <- summary(mediation)
  
  ##### Extract results #####
  pval.all <- as.data.frame(summ.dv$coefficients) %>% 
    rownames_to_column("variable") %>% 
    mutate(model = "DV~IV+mediator")
  pval.all <- as.data.frame(summ.mediator$coefficients) %>% 
    rownames_to_column("variable") %>% 
    mutate(model = "mediator~IV") %>% 
    bind_rows(pval.all)
  pval.all <- as.data.frame(summ.totEffect$coefficients) %>% 
    rownames_to_column("variable") %>% 
    mutate(model = "DV~IV") %>% 
    bind_rows(pval.all) %>% 
    filter(variable != '(Intercept)') %>% 
    rename(Pval ='Pr(>|t|)')
  
  med.summ <- data.frame(
    #ACME = average causal mediation effects
    #ADE = average direct effects
    model="mediation",
    variable = c("ACME","ADE","Total Effect","Prop. Mediated"),
    Estimate = c(mediation.summ$d0, mediation.summ$z0,
                 mediation.summ$tau.coef, mediation.summ$n0),
    
    CI.95.lower =c(mediation.summ$d0.ci[1], mediation.summ$z0.ci[1],
                   mediation.summ$tau.ci[1], mediation.summ$n0.ci[1]),
    
    CI.95_upper = c(mediation.summ$d0.ci[2], mediation.summ$z0.ci[2],
                    mediation.summ$tau.ci[2], mediation.summ$n0.ci[2]),
    
    Pval = c(mediation.summ$d0.p, mediation.summ$z0.p,
                  mediation.summ$tau.p, mediation.summ$n0.p))
    
  results.all <- full_join(pval.all, med.summ, 
                    by = c("variable", "Estimate", "Pval", "model")) %>% 
    dplyr::select(model, everything()) %>% 
    mutate(Pval.symbol = ifelse(Pval <= 0.001, "***",
                                ifelse(Pval <= 0.01, "**",
                                       ifelse(Pval <= 0.05, "*", ""))))
  
  ##### Save results #####
  assign("mediation_results", results.all, envir = .GlobalEnv)
  
  dir.create("results/mediation", showWarnings = FALSE)
  filename <- paste("results/mediation/", paste(dv, iv, mediator, sep="_"),
                 ".csv", sep="")
  write_csv(results.all, filename)

  #### Plot ####
  if(plot){
    print("Creating plot")
    require(DiagrammeR, quietly = TRUE)
    require(DiagrammeRsvg, quietly = TRUE)
    require(rsvg, quietly = TRUE)
    
    ##### Edge labels #####
    #create edge labels with model effect size and p-value symbol
    iv.med.lab <- paste(
      round(results.all[results.all$model == "mediator~IV",
                          'Estimate'], digits=3),
      results.all[results.all$model == "mediator~IV",
                  'Pval.symbol'], 
      sep="")
    
    med.dv.lab <- paste(
      round(results.all[results.all$model == "DV~IV+mediator" &
                        results.all$variable == mediator,
                        'Estimate'], digits=3),
      results.all[results.all$model == "DV~IV+mediator" &
                    results.all$variable == mediator,
                  'Pval.symbol'], 
      sep="")
      
    iv.dv.lab <- paste(
      round(results.all[results.all$model == "DV~IV",
                        'Estimate'], digits=3),
      results.all[results.all$model == "DV~IV",
                  'Pval.symbol'], 
      " (",
      round(results.all[results.all$model == "DV~IV+mediator" &
                          results.all$variable == iv,
                        'Estimate'], digits=3),
      results.all[results.all$model == "DV~IV+mediator" &
                    results.all$variable == iv,
                  'Pval.symbol'], 
      ")",
      sep="")

    ##### Format df for diagram #####
    nodes <- create_node_df(n=3,
                            nodes = c("A","B","C"),
                            type="lower",
                            label = c(iv,mediator,dv),
                            shape = "rectangle",
                            color = "black",
                            style = "filled", fillcolor="white",
                            fixedsize = FALSE,
                            fontcolor = "black")
    edges <- create_edge_df(from = c(1,1,2),
                            to = c(3,2,3),
                            color = "black",
                            label = c(iv.dv.lab, iv.med.lab, med.dv.lab))
    
    graph <- create_graph(nodes_df = nodes,
                          edges_df = edges,
                          attr_theme = "lr")
    #### Save plot ####
    dir.create("figs/mediation", showWarnings = FALSE)
    filename.plot <- paste("figs/mediation/", "mediation_", 
                           paste(dv,iv,mediator, sep="_"), ".pdf", sep="")
    
    title.plot <- results.all %>% 
      filter(model == "mediation") %>% 
      dplyr::select(variable, Estimate, Pval.symbol) %>% 
      mutate(Estimate = round(Estimate, digits=3))
    title.plot <- paste(title.plot$variable, title.plot$Estimate, 
                        title.plot$Pval.symbol, 
                        collapse = ", ", sep=" ")
      
    export_graph(graph=graph, file_name=filename.plot, title=title.plot,
                 height=200, width=800)
  }
}
