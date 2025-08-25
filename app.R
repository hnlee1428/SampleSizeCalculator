library(shiny)
library(ggplot2)
library(tidyverse)
library(DT)
library(plotly)

ui <- fluidPage(
  #titlePanel("Inference & Sampling Toolkit for Finite Populations"),
  tags$head(
    tags$style(HTML("
      .app-title {
        text-align: center;
        font-size: 28px;
        font-weight: 600;
        margin: 30px 0 10px 0;
      }
    "))
  ),
  tags$head(
    tags$style(HTML("
      .nav-tabs > li > a {
        font-size: 18px !important;
      }
    "))
  ),
  div(
    class = "app-title",
    tagList(
      icon("calculator"),  # Try "clipboard", "chart-bar", or another if you prefer
      "Inference & Sampling Toolkit for Finite Populations"
    )
  ),
  
  uiOutput("dynamic_title"), 
  
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tabs == 'inference' || input.tabs == 'plot' || input.tabs == 'table'",  # Show only on 'inference' or 'plot' tab,
        radioButtons("confidence_choice_inference", label = HTML("Choose significance level (&alpha;)<br>(0.01 recommended)"), choices = c(0.01, 0.05)),
        numericInput("N_inference", "Population size (N):", value = 1000, min = 5, step = 1),
        numericInput("n", "Sample size (n):", value = 40, min = 5, step = 1),
        numericInput("x", "Number of negatives in the sample (x) :", min = 0, value = 0, step = 1)
      ),
      conditionalPanel(
        condition = "input.tabs == 'plot' || input.tabs == 'table'",  # Show only on 'Plot' tab,
        h4("Settings for Future Sampling", style = "background-color: greenyellow;"),
        sliderInput("add_x", "Additional negatives beyond your x :", min = 0, max = 5, value = 1, step = 1),
        numericInput("line_y", "Reference line (y-axis):", value = NA)
      )#,
      # conditionalPanel(
      #   condition = "input.tabs == 'table'",  # Show only on 'Plot' tab,
      #   radioButtons("confidence_choice", label = HTML("Choose significance level (&alpha;)<br><small>(0.01 recommended)</small>"), choices = c(0.01, 0.05)),
      #   numericInput("N_table", "Population size (N):", value = 1000, min = 5, step = 1),
      #   numericInput("q_table", "Maximum acceptable negative proportion (q):", value = 0.1, min = 0, max = 0.5),
      #   sliderInput("m_range", "Range of maximum allowed negatives (m) :", min = 0, max = 10, value = c(0, 2), step = 1)
      # )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Description", 
                 withMathJax(),
                 tags$div(style = "padding: 20px; max-width: 800px;",
                          tags$h4("Our inference here relies on the hypergeometric distribution, which models sampling without replacement from a finite population."),
                          
                          tags$h3("Hypergeometric Distribution"),
                          tags$h4("Notation:"),
                          tags$ul(
                            tags$li("\\( N \\): total population size", style = "font-size: 16px;"),
                            tags$li("\\( K \\): number of negatives in the population", style = "font-size: 16px;"),
                            tags$li("\\( n \\): sample size", style = "font-size: 16px;"),
                            tags$li("\\( x \\): number of negatives in the sample", style = "font-size: 16px;")
                          ),
                          
                          tags$h4("Probability mass function (PMF):"),
                          tags$h4("$$ P(X = x | N, K, n) = \\frac{\\binom{K}{x} \\binom{N-K}{n-x}}{\\binom{N}{n}} $$"),
                          
                          tags$h4("\\( P(X = x | N, K, n) \\) is the probability of observing exactly x negative items in a sample of size n, drawn without replacement from a population of N items, of which K are negative items.")
                 )
        ),
        tabPanel("Inference about K", value = "inference", 
                 div(style = "color: #004080; font-size: 18px; font-weight: bold; text-align: left;",
                     uiOutput("inference_out")),
                 div(style = "color: #4d4d4d; font-size: 18px; font-style: italic; text-align: left;",
                     uiOutput("inference_out2"))),
        tabPanel("Update Sample Size", value = "plot", plotlyOutput("n_K_estimation_plot"),
                 HTML("<br>"),
                 div(style = "color: #4d4d4d; font-style: italic; font-size: 18px",
                        textOutput("plot_description")))#,
        #tabPanel("Update Sample Size (Table)", value = "table", DTOutput("n_K_estimation_table"))
      )
    )
  )
)

server <- function(input, output){
  output$dynamic_title <- renderUI({
    if (input$tabs == "inference") {
      h3("Inference Settings")
    } else if (input$tabs == "plot") {
      h3("Inference Settings")
    } else if (input$tabs == "table") {
      h3("Inference Settings")
    }
  })
  
  n.calculate.f <- function(N = 1000, q = 0.1, allowed_no_negatives = 0){
    
    n_seq <- seq(1, N, by = 1)
    verify_n <- (n_seq-allowed_no_negatives)/n_seq
    n_seq <- n_seq[verify_n>q]
    
    no.of.negatives.in.population <- as.integer(N*q) # K
    no.of.positives.in.population <- N-no.of.negatives.in.population 
    
    res_seq <- c()
    for(k in 1:length(n_seq)){
      sample.size <- n_seq[k]
      # P(Type I error)
      x_q <- 0:allowed_no_negatives
      prob <- 0
      for(i in 1:length(x_q)){
        prob <- prob + dhyper(x_q[i],
                              no.of.negatives.in.population,
                              no.of.positives.in.population,
                              sample.size,
                              log = FALSE)
      }
      res_seq <- c(res_seq, 1-prob)
    }
    res_df <- cbind.data.frame(n = n_seq,
                               confidence_level = res_seq,
                               N = N,
                               q = q,
                               allowed_no_negatives = allowed_no_negatives) %>%
      mutate(one_or_not = if_else(round(confidence_level,3) == 1, 1, 0)) %>%
      group_by(one_or_not) %>%
      mutate(row_id = row_number()) %>%
      ungroup() %>%
      dplyr::filter(!(one_or_not == 1 & row_id > 10)) %>%
      dplyr::select(-c("one_or_not", "row_id"))
    return(res_df)
  }
  
  infer_K_at_alpha.f <- function(alpha, N, sample.size, x_negative){
    K_alpha <- 1
    prob <- 1
    res <- c()
    while((prob > alpha) & K_alpha <= N){
      prob <- 0
      for(xx in 0:x_negative){
        prob <- prob + dhyper(xx,
                              K_alpha,
                              N-K_alpha,
                              sample.size,
                              log = FALSE)
      }
      res <- rbind(res, cbind(significance = prob, 
                              confidence = 1-prob, 
                              lowestreportable_K = K_alpha))
      K_alpha <- K_alpha + 1}
    res <- as.data.frame(res)
    return(res$lowestreportable_K[nrow(res)])
  }
  
  
  inference <- reactive({
    req(input$confidence_choice_inference, input$N_inference, input$n, input$x)
    if(input$N_inference >= input$n & input$n >= input$x){
      infer_K_at_alpha.f(alpha = as.numeric(input$confidence_choice_inference),
                         N = input$N_inference, 
                         sample.size = input$n,
                         x_negative = input$x)
    }else{
      return("Make sure that N \U2265 n \U2265 x.")}
    
    
  })
  
  output$inference_out <- renderUI({
    #inference()
    HTML(paste0(
      '<span style="color: #CC0000;"> <br>K</span>',
      '<span style="color: #4d4d4d;">: the number of negative itmes in the population, you want to infer.</span><br>',
      '<span style="color: #4d4d4d;"> As you observed</span>',
      '<span style="color: #1338BE;"> ', input$x, '</span>',
      '<span style="color: #4d4d4d;"> negative items in the sample, you can say:</span><br>',
      '<span style="color: #004080;"> With </span>',
      '<span style="color: #004080;">', as.integer(100*(1-as.numeric(input$confidence_choice_inference))),'</span>',
      '<span style="color: #004080;">% confidence, </span>',
      '<span style="color: #004080; background-color: yellow;">K is smaller than </span>',
      '<span style="color: #004080; background-color: yellow;">', inference() ,'</span>',
      '<span style="color: #004080; background-color: yellow;">. (i.e., K/N < </span>',
      '<span style="color: #004080; background-color: yellow;">', round(inference()/input$N_inference, 4) ,'</span>',
      '<span style="color: #004080; background-color: yellow;">).</span><br><br>'))
  })
  
  output$inference_out2 <- renderUI({
    #inference()
    HTML(paste0(
     '<span style="color: #000000;">Note: Supporting values of K below the current level (',
      '<span style="color: #000000;">', inference() ,'</span>',
      '<span style="color: #000000;">) requires a larger sample size. While additional negatives may be observed, the overall proportion is typically preserved through random sampling.</span>'
    ))
  })
  update_n_inference <- reactive({
    req(input$confidence_choice_inference, input$N_inference, input$n, input$x)
    max.K.to.show <- inference()
    
    allowed_no_negatives <- input$x+input$add_x
    n_seq <- input$n:as.integer(input$N_inference*0.3)
    
    res <- 0:allowed_no_negatives %>%
      map(~c())
    for(i in 0:allowed_no_negatives){
      for(j in 1:length(n_seq)){
        res[[i+1]] <- c(res[[i+1]],
                        infer_K_at_alpha.f(alpha = as.numeric(input$confidence_choice_inference), 
                                           N = input$N_inference, sample.size = n_seq[j], x_negative = i))
      }}
    res <- res %>%
      map(~cbind.data.frame(K.to.show = .x, n.to.show = n_seq) %>%
            dplyr::filter(K.to.show<=max.K.to.show)) %>%
      map2(0:allowed_no_negatives,
           ~cbind.data.frame(.x, potential_x = rep(.y, nrow(.x)))) %>%
      map_dfr(~.x)
    return(res)
    
  })
  
  output$n_K_estimation_plot <- renderPlotly({
    res_data <- update_n_inference()
    y_line <- if(!is.na(input$line_y)){input$line_y}else{max(max(res_data$K.to.show)/input$N_inference)}

    no.group <- length(unique(res_data$potential_x))
    x_ticks <- pretty(res_data$n.to.show, n = 10)
    y_ticks <- seq(0, max(res_data$K.to.show)/input$N_inference, by = 0.01)
    p <- plot_ly(x = ~res_data$n.to.show, y = ~res_data$K.to.show/input$N_inference,
                 type = 'scatter',
                 mode = 'markers',
                 color = ~factor(res_data$potential_x), colors = RColorBrewer::brewer.pal(max(3, no.group), "Set2")[seq_len(no.group)],
                 symbol = ~factor(res_data$potential_x),
                 symbols = c("circle-open", "cross", "x", "diamond-open", "triangle-up",
                             "star", "cross-open", "circle", "x-open", "square")[1:no.group],
                 marker = list(size = 3),
                 showlegend = TRUE) %>%
      plotly::layout(
        xaxis = list(title = "Sample Size",
                     tickvals = x_ticks),
        yaxis = list(title = "Minimum q to support K/N < q",
                     tickvals = y_ticks),
        legend = list(title = list(text = "Possible observed value of x: "),
                      orientation = "h",     # horizontal layout
                      x = 0.5,               # center horizontally
                      xanchor = "center",    # align x to center of legend box
                      y = 1.1,
                      itemsizing = "constant"),                # position above the plot
        bargap = 0.1,
        shapes = list(list( type = "line",
            x0 = 0, x1 = 1,           # full width (in paper coordinates)
            xref = "paper",           # use plot width instead of data scale
            y0 = y_line, y1 = y_line,       # y position for horizontal line
            line = list(
              dash = "dot",          # makes it dotted
              width = 1,
              color = "dodgeblue4"))))
    
    if(no.group == 1) {
      p <- p %>% plotly::add_trace(
        data = NULL,
        x = numeric(0), y = numeric(0),
        name = unique(res_data$potential_x),
        type = "scatter", mode = "markers",
        marker = list(color = RColorBrewer::brewer.pal(max(3, no.group), "Set2")[1]),
        showlegend = TRUE,
        inherit = FALSE 
      )
    }
    p
  })
  
  output$plot_description <- renderText({
    req(input$add_x)
  
      if(input$x+input$add_x == 0){
      "The graph of q versus n shows the required sample size to support K/N < q, assuming no negatives are observed in a future sample."
      }else{
        paste0("Note: The graph of q versus n shows the required sample size to support K/N < q",
               ", assuming no negatives are observed in a future sample, with allowance for up to ",
               input$x + input$add_x, " negatives.")
      }
    
  })
  
  # output$n_K_estimation_table <- renderDT({
  #   pivot_table<-update_n_inference()
  #   datatable(pivot_table,rownames = FALSE, options = list(pageLength = 10))
  # })
  
  # pmf_data <- reactive({
  #   req(input$N, input$q, input$m)
  #   n.calculate.f(N = input$N, q = input$q, allowed_no_negatives = input$m)
  # })
  # 
  # 
  # output$n_confidence_plot <- renderPlotly({
  #   x_ticks <- pretty(pmf_data()$n, n = 20)
  #   y_ticks <- seq(0, 1, by = 0.1)  # or define your own
  #   plot_ly(x = ~pmf_data()$n, y = ~pmf_data()$confidence_level,
  #           type = 'scatter',
  #           mode = 'markers',
  #           marker = list(color = 'darkorange', size = 8)) %>%
  #     plotly::layout(
  #       xaxis = list(title = "Sample Size",
  #                    tickvals = x_ticks),
  #       yaxis = list(title = "Confidence Level",
  #                    tickvals = y_ticks),
  #       bargap = 0.1)
  # })
  # 
  # seq_vals <- reactive({
  #   seq(input$m_range[1], input$m_range[2], by = 1)
  # })
  # 
  # output$pmf_table <- renderDT({
  #   
  #   # table_data <- seq_vals() %>%  
  #   #   map(~n.calculate.f(N = input$N_table, q = input$q_table, allowed_no_negatives = .x) %>%
  #   #         # mutate(cl_group = if_else(confidence_level<(1-as.numeric(input$confidence_choice)), 0, 1)) %>%
  #   #         mutate(cl_group = if_else(confidence_level<0.99, 0, 0.99)) %>%
  #   #         group_by(cl_group) %>%
  #   #         dplyr::filter(n == min(n),
  #   #                       cl_group>0) %>%
  #   #         ungroup()) %>%
  #   #   map_dfr(~.x)
  #   # 
  #   # 
  #   # pivot_table <- table_data %>%
  #   #   group_by(allowed_no_negatives) %>%
  #   #   summarise(n = min(n), cl_group = unique(cl_group), .groups = "drop") %>%
  #   #   pivot_wider(names_from = allowed_no_negatives,
  #   #               values_from = n,
  #   #               names_prefix = "m =") %>%
  #   #   mutate(cl_group = paste(100*cl_group, "%"))
  #   # colnames(pivot_table)[1] <- c("Confidence level = (1-\u03B1)100%")
  #   
  #   pivot_table<-as.data.frame(update_n_inference())
  #   datatable(pivot_table,rownames = FALSE)
  # })
}

shinyApp(ui, server)
