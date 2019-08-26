shinyServer(
  function(input, output, session){
    
    LINKAGE_SAMPLE=1000
    
    # updateCollapse(session,id="what_snSMART", open = NULL, close = "what_snSMART", style = NULL)
    # 
    # updateCollapse(session,id="goal_of_app", open = NULL, close = "goal_of_app", style = NULL)
    
    
    formula <- eventReactive(input$goButton,{
      
      max_beta1=1/min(c(max(input$piA,input$piB,input$piC),0.99))
      str_beta1=paste0("The value of ","\u03B21"," should be within 1~",max_beta1)
      
      validate(
        need(input$piA<1 & input$piA>0, "The value of \u03C0A should be within 0.01~0.99"),
        need(input$piB<1 & input$piB>0, "The value of \u03C0B should be within 0.01~0.99"),
        need(input$piC<1 & input$piC>0, "The value of \u03C0C should be within 0.01~0.99"),
        need((input$beta1<=max_beta1 & input$beta1>1)|
               is.na(max_beta1)|
               input$piA>1 | input$piA<0|
               input$piB>1 | input$piB<0|
               input$piC>1 | input$piC<0
               ,
             str_beta1),
        need(input$beta0<1 & input$beta0>0, "The value of \u03B20 should be within 0.01~0.99"),
        need(input$coverage<1 & input$coverage>0, "The value of 1-\u03B1 should be within 0.01~0.99"),
        need(input$power<1 & input$power>0, "The value of 1-\u03BE should be within 0.01~0.99"),
        
        need(input$muA<1 & input$muA>0, "The value of \u03bcA should be within 0.01~0.99"),
        need(input$muB<1 & input$muB>0, "The value of \u03bcB should be within 0.01~0.99"),
        need(input$muC<1 & input$muC>0, "The value of \u03bcC should be within 0.01~0.99"),
        need(input$nA>0, "The value of nA should be greater than 0"),
        need(input$nB>0, "The value of nB should be greater than 0"),
        need(input$nC>0, "The value of nC should be greater than 0")
      )
      
      K=3
      magic_Z=1.5
      pi_A=as.numeric(input$piA)
      pi_B=as.numeric(input$piB)
      pi_C=as.numeric(input$piC)
      
      muA_input=as.numeric(input$muA)
      muB_input=as.numeric(input$muB)
      muC_input=as.numeric(input$muC)
      
      nA_input=as.numeric(input$nA)
      nB_input=as.numeric(input$nB)
      nC_input=as.numeric(input$nC)
      
      LIST_OF_PIS=c(pi_A,pi_B,pi_C)
      LIST_OF_PIS=LIST_OF_PIS[order(LIST_OF_PIS,decreasing = T)]
      
      COVRAGE=as.numeric(input$coverage)
      # COVRAGE_BONFF=1-(1-COVRAGE)/(K-1)
      # COVRAGE=COVRAGE_BONFF
      POW = as.numeric(input$power) #when POW is small, we need to decrease sample size lower limit, need to estimate the lower limit or there will not be a opposite sign
      
      # SAMPLE_SIZE_LLIMIT=1
      # SAMPLE_SIZE_ULIMIT=1000
      
      # CIL_MIN=0.1# can not be too small
      # CIL_MAX=2
      # CIL_STEP=0.01
      
      SS_LOW=1
      SS_HIGH=300
      CONVERGE_TOL=0.001
      CIL_MIN=0.01# can not be too small
      CIL_STEP_I=0.01
      
      ###################
      # check if pis are all the same
      # see if the first two treatment are the same-if not the same then do computation, if the same, do computation with the first and last
      # need all treatment arm for computation
      if(LIST_OF_PIS[1]==LIST_OF_PIS[2]){
        piA=LIST_OF_PIS[1]
        piB=LIST_OF_PIS[K]
        piC=LIST_OF_PIS[setdiff(1:K,c(1,K))]
      }else{
        piA=LIST_OF_PIS[1]
        piB=LIST_OF_PIS[2]
        piC=LIST_OF_PIS[setdiff(1:K,c(1,2))]
      }
      
      # generate truncated beta1 and beta0
      set.seed(1991)
      # generate pareto truncated at 1.5 with location 1 and scale 3, mean 1.5
      beta1_mean=as.numeric(input$beta1)
      pareto_beta=1/(1-1/beta1_mean)
      # beta1_sample=rtrunc(LINKAGE_SAMPLE, spec="pareto",a=1,b=1/max(c(piA,piB,piC)),1,pareto_beta)
      # beta1_sample=rep(1,LINKAGE_SAMPLE)
      beta0_mean=as.numeric(input$beta0)
      # beta0_sample=rbeta(LINKAGE_SAMPLE,beta0_mean*2,2-beta0_mean*2)
      # beta0_sample=rep(1,LINKAGE_SAMPLE)
      # plot(hist(beta0))
      # load functions
      
      CIL_MAX=(piA-piB)*magic_Z
      beta1_sample=round(mean(rtrunc(99999, spec="pareto",a=1,b=1/max(c(piA,piB,piC)),1,pareto_beta)),3)
      beta0_sample=rep(beta0_mean,1)
      
      sample_size_list_pair1=NULL
      
      error_count_pair1=0
      warn_count_pair1=0
      
      sim_count=0
      
      error_round_pair1=NULL
      warn_round_pair1=NULL
      error_mesg_pair1=NULL
      warn_mesg_pair1=NULL
      
      error_round_pair1=NULL
      warn_round_pair1=NULL
      error_mesg_pair1=NULL
      warn_mesg_pair1=NULL
      
      # calculate priors
      
      beta1=beta1_sample
      beta0=beta0_sample
      
      aA = muA_input*nA_input
      bA = nA_input-aA
      # aA = 0.25*2
      # bA = 2-aA
      cA = cFunc(aA, bA, beta1)
      dA = dFunc(aA, bA, beta1)
      eA = eFunc(aA, bA, beta0, K)
      fA = fFunc(aA, bA, beta0, K)
      
      aB = muB_input*nB_input
      bB = nB_input-aB
      # aB = 0.25*2
      # bB = 2-aB
      cB = cFunc(aB, bB, beta1)
      dB = dFunc(aB, bB, beta1)
      eB = eFunc(aB, bB, beta0, K)
      fB = fFunc(aB, bB, beta0, K)
      
      aC = muC_input*nC_input
      bC = nC_input-aC
      # aC = 0.25*2
      # bC = 2-aC
      cC = cFunc(aC, bC, beta1)
      dC = dFunc(aC, bC, beta1)
      eC = eFunc(aC, bC, beta0, K)
      fC = fFunc(aC, bC, beta0, K)
      
      ciZ=qnorm(1-(1-COVRAGE)/2) # critical value of coverage
      
      # i=0
      
      withProgress(message = 'Sample Size Calculating', value = 0, {
        for(CIL_I in seq(CIL_MAX,CIL_MIN,by=-CIL_STEP_I))
        {
          # i=i+1
          print(CIL_I)
          # print(CIL_I)
          # get sample size solved
          # print(CIL_I)
          # CIL_I=0.26
          if(CIL_I/2==(max(c(piA,piB,piC))-min(c(piA,piB,piC)))){
            next
          }
          tryCatch({
            
            fun <- function (x) {
              2*ciZ*sqrt(var1_o1o2_diff(K=K,
                                        piA=piA, piB=piB, piC=piC,
                                        beta1=beta1, beta0=beta0,
                                        aA=aA, bA=bA, cA=cA, dA=dA, eA=eA, fA=fA,
                                        aB=aB, bB=bB, cB=cB, dB=dB, eB=eB, fB=fB,
                                        aC=aC, bC=bC, cC=cC, dC=dC, eC=eC, fC=fC,
                                        n=x)-mean_o1o2_diff(K=K,
                                                            piA=piA, piB=piB, piC=piC,
                                                            beta1=beta1, beta0=beta0,
                                                            aA=aA, bA=bA, cA=cA, dA=dA, eA=eA, fA=fA,
                                                            aB=aB, bB=bB, cB=cB, dB=dB, eB=eB, fB=fB,
                                                            aC=aC, bC=bC, cC=cC, dC=dC, eC=eC, fC=fC,
                                                            n=x)^2)-CIL_I
            }
            
            # for(i in 1:300){
            #   print(fun(i))
            # }
            
            # ciL=0.28
            
            # sample_size_tmp_pair1=ceiling(uniroot(fun, c(SS_LOW, SS_HIGH),tol = CONVERGE_TOL)$root)
            sample_size_tmp_pair1=ceiling(uniroot(fun, c(SS_LOW, SS_HIGH),tol = CONVERGE_TOL)$root)
            # All <- uniroot.all(sample_size_equation, c(0, 8))
            # calculate powerpower
            # 1-pnorm((ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))+pnorm((-ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))
            mu_o1o2_diff=mean_o1o2_diff(K=K,
                                        piA=piA, piB=piB, piC=piC,
                                        beta1=beta1, beta0=beta0,
                                        aA=aA, bA=bA, cA=cA, dA=dA, eA=eA, fA=fA,
                                        aB=aB, bB=bB, cB=cB, dB=dB, eB=eB, fB=fB,
                                        aC=aC, bC=bC, cC=cC, dC=dC, eC=eC, fC=fC,
                                        n=sample_size_tmp_pair1)
            mu_o1o2_sq_diff=var1_o1o2_diff(K=K,
                                           piA=piA, piB=piB, piC=piC,
                                           beta1=beta1, beta0=beta0,
                                           aA=aA, bA=bA, cA=cA, dA=dA, eA=eA, fA=fA,
                                           aB=aB, bB=bB, cB=cB, dB=dB, eB=eB, fB=fB,
                                           aC=aC, bC=bC, cC=cC, dC=dC, eC=eC, fC=fC,
                                           n=sample_size_tmp_pair1)
            var_o1o2_diff=mu_o1o2_sq_diff-mu_o1o2_diff^2
            # pow_pair1=(1-pnorm((ciL/2-mu_o1o2_diff)/sqrt(var_o1o2_diff)))
            pow_pair1=(1-pnorm((CIL_I/2-mu_o1o2_diff)/sqrt(var_o1o2_diff)))
            
            # total_process=length(seq(CIL_MAX,CIL_MIN,by=-CIL_STEP_I))
            # incProgress(1/total_process, detail = paste0(floor(i/total_process*100)," percent done"))
            # 
            # total_process=POW
            setProgress(pow_pair1/POW, detail = paste0(floor(pow_pair1/POW*100)," percent done"))
            # ?incProgress
            
            if(pow_pair1>POW) break
          },
          error = function(c) {
            # error_round_tmp_pair1=cbind(i,beta1,beta0,CIL_I)
            # error_round_pair1<<-rbind(error_round_pair1,error_round_tmp_pair1)
            # next
          },
          warning = function(c) {
            # warn_round_tmp_pair1=cbind(i,beta1,beta0,CIL_I)
            # warn_round_pair1<<-rbind(warn_round_pair1,warn_round_tmp_pair1)
            # print(i)
            # print(CIL_I)
            # next
          },
          finally = {# posterior_sample_burn<<-window(posterior_sample,start=BURNING, end=MCMC_SAMPLE)
            # posterior_sample_cmb<<-do.call(rbind, posterior_sample_burn)
          }
          )
          print(CIL_I)
          print(pow_pair1)
          print(mu_o1o2_diff)
          print(mu_o1o2_sq_diff)
          print(sample_size_tmp_pair1)
        }
        

      })
      
      sample_size_tmp_pair1

    })
      




 
    
    # output$snSMART_describe <- renderUI(HTML("<ul><li>...text...</li><li>...more text...</li></ul>"))
    
    
    output$myResult=renderUI(
      {
        input$goButton
        result=formula()
        
        parmsettings <- paste0("Design parameters: \n", "Treatment Expectation Settings:", sep = " ")
        # str0 <- isolate(paste0("With given settings, ","the estimated sample size per arm for an snSMART is: ", result, "."))
        
        str0 <- h4(paste0("With given settings, ","the estimated sample size per arm for an snSMART is: ", result, "."))
        
        all_pis=c(as.numeric(input$piA),as.numeric(input$piB),as.numeric(input$piC))
        all_pis=all_pis[order(all_pis,decreasing = T)]
        delta=all_pis[1]-all_pis[2]
        max_response=max(as.numeric(input$piA),as.numeric(input$piB),as.numeric(input$piC))
        type_I_error=1-as.numeric(input$coverage)
        
        str1 <- h4(style='margin-top:-20px;',paste0("This implies that for an snSMART with sample size of ", 
                                                    formula(), " per arm", " (", formula()*3," in total for three agents): "))
        str2 <- h4(paste0("The probability of successfully identifying the best treatment is ",
                          input$power, ", when the difference of response rates between the best and second best treatment is at least ",
                          delta, 
                          " and the response rate of the best treatment is ", max_response, "."))
        # str3 <- paste0("(ii) The probability that we mistakenly claim there is a best treatment when there is no best treatment is ",type_I_error,".")
        
        
        # Text version of results for the output file.
        result_title <- "------------------------\n***     Results:     ***\n------------------------\n\n"
        
        str0t <- paste0("With given settings, ","the estimated sample size per arm for a snSMART is: ", result, ".\n")
        
        str1t <- paste0("This implies that for a snSMART with sample size of ", 
                        formula(), " per arm", " (", formula()*3," in total for three agents): ")
        str2t <- paste0("\nThe probability of successfully identifying the best treatment is ",
                        input$power, " when the difference of response rates between the best \nand second best treatment is at least ",
                        delta, 
                        " and the response rate of the best treatment is ", max_response, ".")
        
        results <<- paste0(result_title, str0t,'\n', str1t, str2t,'\n\n--------------------------------------------------')
        
        # Generate HTML of results
        HTML(paste0(str0,'<br/><br/>',str1, str2,'<br/><br/>'))
      }
    )
    
    
    # When the Refresh button is pressed, reload the entire page, 
    # which resets to the initial parameter settings.
    observeEvent(input$refreshButton, {
      session$reload()
    })
    
    
    
    # When the Download Results button is pressed, generate a text file
    # with all parameter inputs and results of the function
    # Generates and sends the output file back to the button when it is triggered
    output$downloadButton <- downloadHandler(
      filename = function() {
        paste("snSMART_results", ".txt", sep = "")
      },
      content = function(file) {
        
        
        
        
        # formatted text output for snSMART download
        
        title <- "--------------------------------------------------\nParameters and results of snSMART Sample Size App\n--------------------------------------------------\n\n"
        
        # parameters and settings
        
        title_parms <- "------------------------\nParameters and settings:\n------------------------\n\n"
        
        # ------------------------------------------------------------------------------
        # Design settings
        design_title <- "*** Design settings: ***\n\n"
        
        # Treatment Expectation Setting
        
        title_exp  <- "Treatment Expectation Settings: \n"
        
        piAt <- "    Response rate (ranges from 0.01 to 0.99) for treatment A: "
        piBt <- "    Response rate (ranges from 0.01 to 0.99) for treatment B: "
        piCt <- "    Response rate (ranges from 0.01 to 0.99) for treatment C: "
        
        Av <- input$piA
        Bv <- input$piB
        Cv <- input$piC
        
        trt_exp_parms <- paste0(title_exp, " \n", piAt, Av, " \n", piBt, Bv, " \n", piCt, Cv, " \n \n")
        
        # Linkage Parameter Setting
        
        
        title_linkage <- "Linkage Parameter Settings: \n"
        
        beta1t <- "    Linkage parameter (ranges from 1.00 to 1/largest response rate) for first stage responders. 
    (A smaller value leads to more conservative sample size calculation because two stages are less correlated: "
        beta0t <- "    Linkage parameter (ranges from 0.01 to 0.99) for first stage non-responders. 
    A larger value leads to a more conservative sample size calculation because two stages are less correlated: "
        
        beta1v <- input$beta1
        beta0v <- input$beta0
        
        linkage_parms <- paste0(title_linkage, " \n", beta1t, beta1v, " \n", beta0t, beta0v, " \n \n")
        
        
        # Coverage and Power Setting
        
        
        title_conv <- "Coverage and Power Settings: \n"
        
        coveraget <- "    Coverage rate (ranges from 0.01 to 0.99) for the posterior difference of top two treatments: "
        powert <- "    Probability (ranges from 0.01 to 0.99) for identify the best treatment: "
        
        coveragev <- input$coverage
        powerv <- input$power
        
        coverage_parms <- paste0(title_conv, " \n", coveraget, coveragev, " \n", powert, powerv, " \n \n")
        
        
        # ------------------------------------------------------------------------------
        
        
        
        # ------------------------------------------------------------------------------
        # Prior settings
        
        prior_title <- "*** Prior settings: ***\n\n"
        
        # Treatment Prior mean Settings
        
        title_prior_mean <- "Treatment prior mean settings: \n"
        
        muAt <- "    Prior mean (ranges from 0.01 to 0.99) for treatment A: "
        muBt <- "    Prior mean (ranges from 0.01 to 0.99) for treatment B: "
        muCt <- "    Prior mean (ranges from 0.01 to 0.99) for treatment C: "
        
        muAv <- input$muA
        muBv <- input$muB
        muCv <- input$muC
        
        prior_mean_parms <- paste0(prior_title, title_prior_mean, " \n", muAt, muAv, " \n", muBt, muBv, " \n", muCt, muCv, " \n \n")
        
        
        
        # Treatment Prior Sample Size Settings
        
        title_prior_ss <- "Treatment prior sample size settings*: \n  (*Note: The larger value of prior sample size, the higher degree of belief for prior means.)\n"
        
        nAt <- "    Prior sample size (larger than 0) for treatment A: "
        nBt <- "    Prior sample size (larger than 0) for treatment B: "
        nCt <- "    Prior sample size (larger than 0) for treatment C: "
        
        nAv <- input$nA
        nBv <- input$nB
        nCv <- input$nC
        
        prior_ss_parms <- paste0(title_prior_ss, " \n", nAt, nAv, " \n", nBt, nBv, " \n", nCt, nCv, " \n \n")
        
        
        
        # ------------------------------------------------------------------------------
        
        
        # ------------------------------------------------------------------------------
        # Full parameters output
        parms <- paste0(title, title_parms, design_title, trt_exp_parms, linkage_parms, 
                        coverage_parms, prior_mean_parms, prior_ss_parms, "\n\n")
        
        # ------------------------------------------------------------------------------
        
        
        
        
        
        
        
        writeLines(paste0(parms, " \n", results), file)
      }
    )
    
    

    
    # End SHINY server
  }
)