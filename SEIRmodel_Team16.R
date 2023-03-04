library(shiny)
library(deSolve)
library(dplyr)

#e0 = 1400 means that some exposed people traveled to the city, thus the spread of the disease started. Therefore they are not counted towards total population in week0, which is 2040000.
best_beta = 0.0002
best_gamma = 0.95
best_theta = 0.005
best_e0 = 1400

COVIDData <- read_excel("../data/COVIDdata.xlsm")
actual_new_infection = COVIDData$`Newly Positive (Confirmed Infections)`
actual_new_infection <- actual_new_infection[!is.na(actual_new_infection)]
actual_new_infection
SEIRmodel <- function(times, state, parameters) {
  with (as.list(c(state, parameters)), {
    dS <- - beta * I * S
    dE <- beta * I * S - theta * E
    dI <- theta * E - gamma * I
    dR <- gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}

predict <- function(beta = 0.000002565, gamma = 0.8, theta = 0.3, S0 = 1886876, E0 = 5000, I0 = 1033, R0 = 147091, 
                    tmax = 21) {
  times <- seq(0, tmax, 1)
  
  parameters = c(beta = beta, gamma = gamma, theta = theta)
  
  iniState <- c(S = S0, E = E0, I = I0, R = R0)
  
  yhat <- as.data.frame(ode(iniState, times, SEIRmodel, parameters))
  yhat$`pred_CumulativeInfections` <- yhat$I + yhat$R
  
  yhat <- yhat %>%
    arrange(time) %>%
    mutate(`New Infected` = `pred_CumulativeInfections` - lag(`pred_CumulativeInfections`, 
                                                              default =first(`pred_CumulativeInfections`)))
  
  #  yhat$`New Infected` <- theta * yhat$E4
  
  return(yhat)
}

calc_RMSE <- function(pred_new_infection, actual_new_infection = actual_new_infection){
  sqrt(mean((pred_new_infection - actual_new_infection)^2))
}

best_result <- predict(beta = best_beta, gamma = best_gamma, theta = best_theta, S0 = 2040000, E0 = best_e0, I0 = 0, R0 = 0,tmax = 21)
pred_new_infection <- best_result$`New Infected`[2:22]
best_rmse = calc_RMSE(pred_new_infection = pred_new_infection, actual_new_infection = actual_new_infection)
print(best_result)
print(best_result$`New Infected`)

plotSEIR <- function(beta = best_beta, gamma = best_gamma, theta = best_theta, S0 = 2040000, E0 = best_e0, I0 = 0, R0 = 0, 
                    tmax = 21) {
  times <- seq(0, tmax, 1)
  
  parameters = c(beta = beta, gamma = gamma, theta = theta)
  
  iniState <- c(S = S0, E = E0, I = I0, R = R0)
  
  solution <- ode(iniState, times, SEIRmodel, parameters)
  
  par(mar = c(5, 5, 0, 0), oma = c(0, 0, 1, 1), mgp = c(2.5, 1, 0), xpd = T)
  plot(NA, xlab = "Time", ylab = "Number of individuals", ylim = c(0, 2040000), 
       xlim = c(0, 60), cex.lab = 2)
  lines(x = solution[, "time"], y = solution[, "S"], col = "black", lwd = 3)
  lines(x = solution[, "time"], y = solution[, "E"], col = "yellow", lwd = 3)
  lines(x = solution[, "time"], y = solution[, "I"], col = "red", lwd = 3)
  lines(x = solution[, "time"], y = solution[, "R"], col = "blue", lwd = 3)
  
  legend(x = 'top', y = 1100, legend = c("S", "E", "I", "R"), 
         col = c("black", "yellow", "red", "blue"), lwd = 3, horiz = T, bty = "n", cex = 2) 
  
  }

SEIRmodel <- function(times, state, parameters) {
  with (as.list(c(state, parameters)), {
    dS <- - beta * I * S
    dE <- beta * I * S - theta * E
    dI <- theta * E - gamma * I
    dR <- gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}

app <- shinyApp(
  ui = fluidPage(
    titlePanel("SEIR Model"),
    
    sidebarLayout(
      sidebarPanel(      
        sliderInput("transmission", "Transmission Rate:",
                    min = 0, max = 0.001,
                    value = best_beta, step = 0.0001),
        
        sliderInput("incubation", "Incubation Rate:",
                    min = 0, max = 1,
                    value = best_theta, step = 0.001),
        
        sliderInput("recovery", "Recovery Rate:",
                    min = 0, max = 1,
                    value = best_gamma, step = 0.01),

        sliderInput("S0", "Initial Susceptible Individuals:",
                    min = 1, max = 2040000,
                    value = 2040000, step = 1000),
        
        sliderInput("E0", "Initial Exposed Individuals:",
                    min = 0, max = 2040000,
                    value = best_e0, step = 100),
        
        sliderInput("I0", "Initial Infected Individuals:",
                    min = 0, max = 2040000,
                    value = 0, step = 100),
        
        sliderInput("tmax", "Time:",
                    min = 0, max = 30,
                    value = 21, step = 1)
        
      ),
      
      mainPanel(
        plotOutput("plot")
      ),
      
      "left"
    )
  ),
  
  server = function(input, output) {
    output$plot <- renderPlot({
      plotSEIR(beta = input$transmission, theta = input$incubation, gamma = input$recovery, E0 = input$E0, S0 = input$S0,I0 = input$I0, R0 = 0, tmax = input$tmax)
    })
  }
)

runApp(app)
cat('The final best RMSE = ', best_rmse)