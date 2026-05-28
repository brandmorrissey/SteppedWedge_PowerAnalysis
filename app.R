required_packages <- c(
  "shiny",
  "dplyr",
  "lubridate",
  "purrr",
  "tidyr",
  "ggplot2",
  "scales",
  "sandwich",
  "lmtest"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {
  
  if (!(pkg %in% installed)) {
    install.packages(pkg)
  }
  
  library(pkg, character.only = TRUE)
}

set.seed(123)

parse_numeric_vector <- function(x) {
  as.numeric(trimws(unlist(strsplit(x, ","))))
}

make_site_schedule <- function(n_sites, active_starts, active_length_months, step_id) {
  sites <- paste0("Site_", seq_len(n_sites))
  
  tibble(
    site = sites,
    step = paste0("Step ", step_id),
    start_active = active_starts
  ) %>%
    mutate(
      end_active = start_active %m+% months(active_length_months),
      site = factor(site, levels = sites)
    )
}

make_design <- function(
    n_sites,
    n_per_site_month,
    trial_start,
    trial_end,
    active_starts,
    active_length_months,
    visit_months,
    step_id
) {
  
  sites <- paste0("Site_", seq_len(n_sites))
  
  trial_months <- seq(
    from = trial_start,
    to = trial_end,
    by = "month"
  )
  
  intervention_start <- make_site_schedule(
    n_sites = n_sites,
    active_starts = active_starts,
    active_length_months = active_length_months,
    step_id = step_id
  )
  
  enrollment_dat <- expand.grid(
    site = sites,
    entry_month = trial_months,
    person = seq_len(n_per_site_month)
  ) %>%
    as_tibble() %>%
    mutate(
      person_id = paste(site, entry_month, person, sep = "_")
    )
  
  followup_dat <- expand.grid(
    person_id = enrollment_dat$person_id,
    visit_month = visit_months
  ) %>%
    as_tibble() %>%
    left_join(enrollment_dat, by = "person_id") %>%
    left_join(intervention_start, by = "site") %>%
    mutate(
      observation_month = entry_month %m+% months(visit_month),
      
      intervention = ifelse(
        observation_month >= start_active &
          observation_month < end_active,
        1,
        0
      ),
      
      site = factor(site),
      month_factor = factor(observation_month),
      visit_factor = factor(visit_month)
    ) %>%
    filter(observation_month <= trial_end)
  
  followup_dat
}

simulate_one <- function(
    n_sites,
    n_per_site_month,
    trial_start,
    trial_end,
    active_starts,
    active_length_months,
    visit_months,
    step_id,
    p_control,
    p_intervention,
    ICC,
    alpha
) {
  
  dat <- make_design(
    n_sites = n_sites,
    n_per_site_month = n_per_site_month,
    trial_start = trial_start,
    trial_end = trial_end,
    active_starts = active_starts,
    active_length_months = active_length_months,
    visit_months = visit_months,
    step_id = step_id
  )
  
  sites <- levels(factor(dat$site))
  
  beta_0 <- qlogis(p_control)
  beta_intervention <- qlogis(p_intervention) - qlogis(p_control)
  
  site_variance <- (ICC * (pi^2 / 3)) / (1 - ICC)
  site_sd <- sqrt(site_variance)
  
  site_effects <- rnorm(length(sites), mean = 0, sd = site_sd)
  
  dat <- dat %>%
    mutate(
      site_re = site_effects[as.numeric(site)],
      linear_pred = beta_0 + beta_intervention * intervention + site_re,
      prob = plogis(linear_pred),
      outcome = rbinom(n(), size = 1, prob = prob)
    )
  
  fit <- tryCatch(
    glm(
      outcome ~ intervention + month_factor + visit_factor + site,
      data = dat,
      family = binomial
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NA)
  
  robust_vcov <- tryCatch(
    sandwich::vcovCL(
      fit,
      cluster = dat$site
    ),
    error = function(e) NULL
  )
  
  if (is.null(robust_vcov)) return(NA)
  
  test <- tryCatch(
    lmtest::coeftest(
      fit,
      vcov = robust_vcov
    ),
    error = function(e) NULL
  )
  
  if (is.null(test)) return(NA)
  
  p_value <- test["intervention", "Pr(>|z|)"]
  
  p_value < alpha
}

estimate_power <- function(
    n_sims,
    n_sites,
    n_per_site_month,
    trial_start,
    trial_end,
    active_starts,
    active_length_months,
    visit_months,
    step_id,
    p_control,
    p_intervention,
    ICC,
    alpha
) {
  
  mean(
    replicate(
      n_sims,
      simulate_one(
        n_sites = n_sites,
        n_per_site_month = n_per_site_month,
        trial_start = trial_start,
        trial_end = trial_end,
        active_starts = active_starts,
        active_length_months = active_length_months,
        visit_months = visit_months,
        step_id = step_id,
        p_control = p_control,
        p_intervention = p_intervention,
        ICC = ICC,
        alpha = alpha
      )
    ),
    na.rm = TRUE
  )
}

make_wedge_plot_data <- function(
    n_sites,
    trial_start,
    trial_end,
    active_starts,
    active_length_months,
    step_id
) {
  
  sites <- paste0("Site_", seq_len(n_sites))
  
  schedule <- make_site_schedule(
    n_sites = n_sites,
    active_starts = active_starts,
    active_length_months = active_length_months,
    step_id = step_id
  )
  
  full_period <- tibble(
    site = factor(sites, levels = rev(sites)),
    step = paste0("Step ", rev(step_id)),
    start = trial_start,
    end = trial_end %m+% months(1),
    period = "Control / usual care"
  )
  
  active_period <- schedule %>%
    transmute(
      site = factor(as.character(site), levels = rev(sites)),
      step = step,
      start = start_active,
      end = pmin(end_active, trial_end %m+% months(1)),
      period = "Active intervention"
    )
  
  bind_rows(full_period, active_period)
}

ui <- fluidPage(
  
  titlePanel("Stepped-Wedge Power Analysis"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4("Study design"),
      
      textInput(
        "sites_per_step",
        "Sites per step",
        value = "2, 2, 1"
      ),
      helpText("Number of sites assigned to each rollout step. Example: 2, 2, 1 means Step 1 has two sites, Step 2 has two sites, and Step 3 has one site."),
      
      textInput(
        "step_offsets",
        "Step start offsets, in months",
        value = "3, 6, 9"
      ),
      helpText("Months after the study start when each rollout step begins. Enter one value per step."),
      
      dateInput(
        "trial_start",
        "Study start month",
        value = as.Date("2026-07-01")
      ),
      helpText("The first month when participant enrollment can occur."),
      
      dateInput(
        "trial_end",
        "Study end month",
        value = as.Date("2029-03-01")
      ),
      helpText("The final month included in the study period."),
      
      numericInput(
        "active_length",
        "Length of active intervention period, in months",
        value = 12,
        min = 1,
        step = 1
      ),
      helpText("How long each site remains in the active intervention period after its step begins."),
      
      textInput(
        "visit_months",
        "Follow-up visit months",
        value = "0, 3, 6"
      ),
      helpText("Repeated measurement points for each participant. Example: 0, 3, 6 means baseline, 3-month, and 6-month follow-up."),
      
      hr(),
      
      h4("Power assumptions"),
      
      numericInput(
        "p_control",
        "Control probability",
        value = 0.60,
        min = 0.01,
        max = 0.99,
        step = 0.01
      ),
      helpText("The expected outcome probability before the intervention or under usual care."),
      
      numericInput(
        "p_intervention",
        "Intervention probability",
        value = 0.80,
        min = 0.01,
        max = 0.99,
        step = 0.01
      ),
      helpText("The expected outcome probability during the active intervention period."),
      
      numericInput(
        "icc",
        "Intraclass correlation coefficient, ICC",
        value = 0.03,
        min = 0,
        max = 0.30,
        step = 0.01
      ),
      helpText("How similar outcomes are within the same site. Higher ICC usually lowers power."),
      
      numericInput(
        "alpha",
        "Significance level",
        value = 0.05,
        min = 0.001,
        max = 0.20,
        step = 0.01
      ),
      helpText("The p-value threshold used to determine whether the intervention effect is statistically significant."),
      
      textInput(
        "sample_sizes",
        "Participants per site per month",
        value = "5, 8, 10, 12, 15, 20"
      ),
      helpText("Recruitment rates to test. Enter values separated by commas."),
      
      numericInput(
        "n_sims",
        "Number of simulations",
        value = 100,
        min = 10,
        step = 10
      ),
      helpText("More simulations give more stable power estimates but take longer to run."),
      
      actionButton("run", "Run Power Analysis")
    ),
    
    mainPanel(
      h3("Stepped-Wedge Design"),
      plotOutput("wedge_plot", height = "350px"),
      
      h3("Estimated Power"),
      plotOutput("power_plot", height = "350px"),
      
      h3("Power Results"),
      tableOutput("power_table")
    )
  )
)

server <- function(input, output, session) {
  
  design_inputs <- reactive({
    
    sites_per_step <- parse_numeric_vector(input$sites_per_step)
    step_offsets <- parse_numeric_vector(input$step_offsets)
    visit_months <- parse_numeric_vector(input$visit_months)
    
    validate(
      need(all(!is.na(sites_per_step)), "Sites per step must be numeric."),
      need(all(sites_per_step > 0), "Each step must include at least one site."),
      need(all(sites_per_step == floor(sites_per_step)), "Sites per step must be whole numbers."),
      need(all(!is.na(step_offsets)), "Step start offsets must be numeric."),
      need(length(sites_per_step) == length(step_offsets),
           "Enter one step start offset for each step."),
      need(all(!is.na(visit_months)), "Visit months must be numeric."),
      need(input$trial_end > input$trial_start,
           "Study end date must be after study start date."),
      need(input$p_intervention > input$p_control,
           "Intervention probability should be greater than control probability.")
    )
    
    n_sites <- sum(sites_per_step)
    
    step_id <- rep(seq_along(sites_per_step), times = sites_per_step)
    
    active_starts <- input$trial_start %m+% months(
      rep(step_offsets, times = sites_per_step)
    )
    
    list(
      n_sites = n_sites,
      step_id = step_id,
      sites_per_step = sites_per_step,
      trial_start = input$trial_start,
      trial_end = input$trial_end,
      active_starts = active_starts,
      active_length_months = input$active_length,
      visit_months = visit_months
    )
  })
  
  output$wedge_plot <- renderPlot({
    
    d <- design_inputs()
    
    wedge_dat <- make_wedge_plot_data(
      n_sites = d$n_sites,
      trial_start = d$trial_start,
      trial_end = d$trial_end,
      active_starts = d$active_starts,
      active_length_months = d$active_length_months,
      step_id = d$step_id
    )
    
    ggplot(wedge_dat) +
      geom_rect(
        aes(
          xmin = start,
          xmax = end,
          ymin = as.numeric(site) - 0.35,
          ymax = as.numeric(site) + 0.35,
          fill = period
        )
      ) +
      scale_y_continuous(
        breaks = seq_along(levels(wedge_dat$site)),
        labels = levels(wedge_dat$site)
      ) +
      scale_x_date(
        date_breaks = "1 year",
        date_labels = "%Y",
        expand = expansion(mult = c(0.01, 0.03))
      ) +
      labs(
        x = "Year",
        y = "Site",
        fill = NULL
      ) +
      theme_minimal()
  })
  
  power_results <- eventReactive(input$run, {
    
    d <- design_inputs()
    sample_sizes <- parse_numeric_vector(input$sample_sizes)
    
    validate(
      need(all(!is.na(sample_sizes)), "Sample sizes must be numeric."),
      need(all(sample_sizes > 0), "Sample sizes must be greater than zero.")
    )
    
    trial_month_count <- length(seq(d$trial_start, d$trial_end, by = "month"))
    
    withProgress(message = "Running simulations...", value = 0, {
      
      tibble(
        n_per_site_month = sample_sizes
      ) %>%
        mutate(
          total_sites = d$n_sites,
          total_enrolled = n_per_site_month * d$n_sites * trial_month_count,
          total_observations = total_enrolled * length(d$visit_months),
          
          power = map_dbl(
            n_per_site_month,
            function(n) {
              
              incProgress(1 / length(sample_sizes))
              
              estimate_power(
                n_sims = input$n_sims,
                n_sites = d$n_sites,
                n_per_site_month = n,
                trial_start = d$trial_start,
                trial_end = d$trial_end,
                active_starts = d$active_starts,
                active_length_months = d$active_length_months,
                visit_months = d$visit_months,
                step_id = d$step_id,
                p_control = input$p_control,
                p_intervention = input$p_intervention,
                ICC = input$icc,
                alpha = input$alpha
              )
            }
          )
        )
    })
  })
  
  output$power_table <- renderTable({
    power_results() %>%
      mutate(
        power = round(power, 3)
      )
  })
  
  output$power_plot <- renderPlot({
    
    ggplot(
      power_results(),
      aes(x = n_per_site_month, y = power)
    ) +
      geom_line() +
      geom_point(size = 3) +
      geom_hline(yintercept = 0.80, linetype = "dashed") +
      scale_y_continuous(
        limits = c(0, 1),
        labels = percent_format(accuracy = 1)
      ) +
      labs(
        x = "Participants per Site per Month",
        y = "Estimated Power"
      ) +
      theme_minimal()
  })
}

shinyApp(ui, server)
