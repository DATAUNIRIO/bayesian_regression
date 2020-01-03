#https://www.r-bloggers.com/creating-and-using-a-simple-bayesian-linear-model-in-brms-and-r/
# https://github.com/tw0handt0uch/summerofblogdown/blob/master/content/post/2019-12-01-creating-and-using-a-simple-bayesian-linear-model-in-brms-and-r.Rmd
# https://rileyking.netlify.com/post/creating-and-using-a-simple-bayesian-linear-model-in-brms-and-r/

# Creating and Using a Simple, Bayesian Linear Model (in brms and R) 
  
library(tidyverse)
library(styler)
library(ggExtra)
library(knitr)
library(brms)
library(cowplot)
library(gridExtra)
library(skimr)
library(DiagrammeR)
library(rayshader)
library(av)
library(rgl)

setwd('C:/Users/Steven/Documents/GitHub/bayesian_regression/Using a Simple Bayesian Linear Model/')
ablation_dta_tbl <- read.csv(file = "abl_data_2.csv")
ablation_dta_tbl <- ablation_dta_tbl %>% select(temp, time)
ablation_dta_tbl %>% skim()

# Analise Exploratoria de Dados
# Diagrama de Dispersao

scatter_1_fig <- ablation_dta_tbl %>% ggplot(aes(x = time, y = temp)) +
  geom_point(
    colour = "#2c3e50",
    fill = "#2c3e50",
    size = 2,
    alpha = 0.4
  ) +
  labs(
    x = "Ablation Time (seconds)",
    y = "Tissue Temperature (deg C)",
    title = "Ablation Time vs. Tissue Temperature",
    subtitle = "Simulated Catheter RF Ablation"
  )

scatter_hist_1_fig <- ggMarginal(scatter_1_fig,
                                 type = "histogram",
                                 color = "white",
                                 alpha = 0.7,
                                 fill = "#2c3e50",
                                 xparams = list(binwidth = 1),
                                 yparams = list(binwidth = 2.5)
)

# ggExtra needs these explit calls to display in Markdown docs *shrug*
grid::grid.newpage()
grid::grid.draw(scatter_hist_1_fig)




grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle, fillcolor = yellow]        
      tab1 [label = 'Step 1: Propose a distribution for the response variable \n Choose a maximum entropy distribution given the constraints you understand']
      tab2 [label = 'Step 2: Parameterize the mean \n The mean of the response distribution will vary linearly across the range of predictor values']
      tab3 [label = 'Step 3: Set priors \n Simulate what the model knows before seeing the data.  Use domain knowledge as constraints.']
      tab4 [label = 'Step 4: Define the model \n Create the model using the observed data, the likelihood function, and the priors']
      tab5 [label = 'Step 5: Draw from the posterior \n Plot plausible lines using parameters visited by the Markov chains']
      tab6 [label = 'Step 6: Push the parameters back through the model \n Simulate real data from plausible combinations of mean and sigma']
      # edge definitions with the node IDs
      tab1 -> tab2 -> tab3 -> tab4 -> tab5 -> tab6;
      }
      ")

#----------------------------------------------------------------------------------------
# Step 1: Set priors 
#----------------------------------------------------------------------------------------

#We know some things about these data and we should use it to help regularize to model through the priors.
#Temperature is a continuous variable so we want a continuous distribution. We also know from the nature of the treatment that there isn’t really any physical mechanism within the device that would be expected to cool down the tissue below normal body temperature. Since only heating is expected, the slope should be positive or zero.
# McElreath emphasizes simulating from the priors to visualize “what the model knows before it sees the data”. Here are some priors to consider. Let’s evaluate.

# Set seed for repeatability
set.seed(1999)

# number of sims
n <- 150

# random draws from the specified prior distributions
# lognormal distribution is used to constrain slopes to positive values
a <- rnorm(n, 75, 15)

b <- rnorm(n, 0, 1)
b_ <- rlnorm(n, 0, 0.8)

# calc mean of time and temp for later use
mean_temp <- mean(ablation_dta_tbl$temp)
mean_time <- mean(ablation_dta_tbl$time)

# dummy tibble to feed ggplot()
empty_tbl <- tibble(x = 0)

# y = b(x - mean(var_1)) + a is equivalent to:
# y = bx + (a - b * mean(var_1))

# in this fig we use the uninformed prior that generates some unrealistic values
prior_fig_1 <- empty_tbl %>% ggplot() +
  geom_abline(
    intercept = a - b * mean_time,
    slope = b,
    color = "#2c3e50",
    alpha = 0.3,
    size = 1
  ) +
  ylim(c(0, 150)) +
  xlim(c(0, 150)) +
  labs(
    x = "time (sec)",
    y = "Temp (C)",
    title = "Prior Predictive Simulations",
    subtitle = "Uninformed Prior"
  )

# in this fig we confine the slopes to broad ranges informed by what we know about the domain
prior_fig_2 <- empty_tbl %>% ggplot() +
  geom_abline(
    intercept = a - b_ * mean_time,
    slope = b_,
    color = "#2c3e50",
    alpha = 0.3,
    size = 1
  ) +
  ylim(c(0, 150)) +
  xlim(c(0, 150)) +
  labs(
    x = "time (sec)",
    y = "Temp (C)",
    title = "Prior Predictive Simulations",
    subtitle = "Mildly Informed Prior"
  )

plot_grid(prior_fig_1, prior_fig_2)


empty_tbl %>% ggplot() +
  geom_abline(
    intercept = a - b_ * mean_time,
    slope = b_,
    color = "#2c3e50",
    alpha = 0.3,
    size = 1
  ) +
  ylim(c(37, 100)) +
  xlim(c(15, 40)) +
  labs(
    x = "time (sec)",
    y = "Temp (C)",
    title = "Prior Predictive Simulations",
    subtitle = "Mildly Informed Prior, Original Data Range"
  )


#----------------------------------------------------------------------------------------
# Step 2: Define the model
#----------------------------------------------------------------------------------------

#Here I use the brm() function in brms to build what I’m creatively calling: “model_1”. This one uses the un-centered data for time. This function uses Markov Chain Monte Carlo to survey the parameter space. After the warm up cycles, the relative amount of time the chains spend at each parameter value is a good approximation of the true posterior distribution. I’m using a lot of warm up cycles because I’ve heard chains for the uniform priors on sigma can take a long time to converge. This model still takes a bit of time to chug through the parameter space on my modest laptop.

model_1 <-
 brm(
   data = ablation_dta_tbl, family = gaussian,
   temp ~ 1 + time,
   prior = c(
     prior(normal(75, 15), class = Intercept),
     prior(lognormal(0, .8), class = b),
     prior(uniform(0, 30), class = sigma)
   ),
   iter = 41000, warmup = 40000, chains = 4, cores = 4,
   seed = 4
 )


#----------------------------------------------------------------------------------------
# Step 3: Draw from the posterior
#----------------------------------------------------------------------------------------
# The fruits of all my labor! The posterior holds credible combinations for sigma and the slope and intercept (which together describe the mean of the outcome variable we care about). Let’s take a look.

post_samplesM1_tbl <-
  posterior_samples(model_1) %>%
  select(-lp__) %>%
  round(digits = 3)

post_samplesM1_tbl %>%
  head(10) %>%
  kable(align = rep("c", 3))

plot(model_1)


#Posterior_summary() can grab the model results in table form.

mod_1_summary_tbl <-
  posterior_summary(model_1) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>%
  mutate_at(.vars = c(2:5), funs(as.numeric(.)))

mod_1_summary_tbl %>%
  kable(align = rep("c", 5))

#----------------------------------------------------------------------------------------
#Now let’s see what changes if the time data is centered. Everything is the same here in model_2 except the time_c data which is transformed by subtracting the mean from each value.
#----------------------------------------------------------------------------------------

ablation_dta_tbl <- ablation_dta_tbl %>% mutate(time_c = time - mean(time))

model_2 <-
 brm(
   data = ablation_dta_tbl, family = gaussian,
   temp ~ 1 + time_c,
   prior = c(
     prior(normal(75, 15), class = Intercept),
     prior(lognormal(0, .8), class = b),
     prior(uniform(0, 30), class = sigma)
   ),
   iter = 41000, warmup = 40000, chains = 4, cores = 4,
   seed = 4
 )

#Plotting model_2 to compare with the output of model_1 above.

plot_mod_2_fig <- plot(model_2)


post_samplesM2_tbl <-
  posterior_samples(model_2) %>%
  select(-lp__) %>%
  round(digits = 3)

post_samplesM2_tbl %>%
  head(10) %>%
  kable(align = rep("c", 3))

#----------------------------------------------------------------------------------------
# Visualize the original data (centered and un-centered versions) along with plausible values for regression line of the mean:
#----------------------------------------------------------------------------------------

mean_regressionM1_fig <-
  ablation_dta_tbl %>%
  ggplot(aes(x = time, y = temp)) +
  geom_point(
    colour = "#481567FF",
    size = 2,
    alpha = 0.6
  ) +
  geom_abline(aes(intercept = b_Intercept, slope = b_time),
              data = post_samplesM1_tbl,
              alpha = 0.1, color = "gray50"
  ) +
  geom_abline(
    slope = mean(post_samplesM1_tbl$b_time),
    intercept = mean(post_samplesM1_tbl$b_Intercept),
    color = "blue", size = 1
  ) +
  labs(
    title = "Regression Line Representing Mean of Slope",
    subtitle = "Data is As-Observed (No Centering of Predictor)",
    x = "Time (s)",
    y = "Temperature (C)"
  )

mean_regressionM2_fig <-
  ablation_dta_tbl %>%
  ggplot(aes(x = time_c, y = temp)) +
  geom_point(
    color = "#55C667FF",
    size = 2,
    alpha = 0.6
  ) +
  geom_abline(aes(intercept = b_Intercept, slope = b_time_c),
              data = post_samplesM2_tbl,
              alpha = 0.1, color = "gray50"
  ) +
  geom_abline(
    slope = mean(post_samplesM2_tbl$b_time_c),
    intercept = mean(post_samplesM2_tbl$b_Intercept),
    color = "blue", size = 1
  ) +
  labs(
    title = "Regression Line Representing Mean of Slope",
    subtitle = "Predictor Data (Time) is Centered",
    x = "Time (Difference from Mean Time in seconds)",
    y = "Temperature (C)"
  )


combined_mean_fig <-
  ablation_dta_tbl %>%
  ggplot(aes(x = time, y = temp)) +
  geom_point(
    colour = "#481567FF",
    size = 2,
    alpha = 0.6
  ) +
  geom_point(
    data = ablation_dta_tbl, aes(x = time_c, y = temp),
    colour = "#55C667FF",
    size = 2,
    alpha = 0.6
  ) +
  geom_abline(aes(intercept = b_Intercept, slope = b_time),
              data = post_samplesM1_tbl,
              alpha = 0.1, color = "gray50"
  ) +
  geom_abline(
    slope = mean(post_samplesM1_tbl$b_time),
    intercept = mean(post_samplesM1_tbl$b_Intercept),
    color = "blue", size = 1
  ) +
  geom_abline(aes(intercept = b_Intercept, slope = b_time_c),
              data = post_samplesM2_tbl,
              alpha = 0.1, color = "gray50"
  ) +
  geom_abline(
    slope = mean(post_samplesM2_tbl$b_time_c),
    intercept = mean(post_samplesM2_tbl$b_Intercept),
    color = "blue", size = 1
  ) +
  labs(
    title = "Regression Line Representing Mean of Slope",
    subtitle = "Centered and Un-Centered Predictor Data",
    x = "Time (s)",
    y = "Temperature (C)"
  )

combined_predicts_fig <- combined_mean_fig + 
  ylim(c(56,90)) +
  labs(title = "Points Represent Observed Data (Green is Centered)",
       subtitle = "Regression Line Represents Rate of Change of Mean (Grey Bands are Uncertainty)")
combined_predicts_fig



p_spaceM1_fig <- 
  post_samplesM1_tbl[1:1000, ] %>%
  ggplot(aes(x = b_time, y = b_Intercept, color = sigma)) +
  geom_point(alpha = 0.5) +
  geom_density2d(color = "gray30") +
  scale_color_viridis_c() +
  labs(
    title = "Parameter Space - Model 1 (Un-Centered)",
    subtitle = "Intercept Represents the Expected Temp at Time = 0"
  )

p_spaceM1_fig

p_spaceM2_fig <- 
  post_samplesM2_tbl[1:1000, ] %>%
  ggplot(aes(x = b_time_c, y = b_Intercept, color = sigma)) +
  geom_point(alpha = 0.5) +
  geom_density2d(color = "gray30") +
  scale_color_viridis_c() +
  labs(
    title = "Parameter Space - Model 2 (Centered)",
    subtitle = "Intercept Represents the Expected Temp at Mean Time"
  )

p_spaceM2_fig 
#ggsave(filename = "p_spaceM2_fig.png")

p_spaceC_tbl <- 
  post_samplesM2_tbl[1:1000, ] %>%
  ggplot(aes(x = b_time_c, y = b_Intercept, color = sigma)) +
  geom_point(alpha = 0.5) +
  geom_point(data = post_samplesM1_tbl, aes(x = b_time, y = b_Intercept, color = sigma), alpha = 0.5) +
  scale_color_viridis_c() +
  labs(
    title = "Credible Parameter Values for Models 1 and 2",
    subtitle = "Model 1 is Un-Centered, Model 2 is Centered",
    x = expression(beta["time"]),
    y = expression(alpha["Intercept"])) +
  ylim(c(54, 80))

p_spaceC_tbl


par(mfrow = c(1, 1))
plot_gg(p_spaceC_tbl, width = 5, height = 4, scale = 300, multicore = TRUE, windowsize = c(1200, 960),
        fov = 70, zoom = 0.45, theta = 330, phi = 40)

Sys.sleep(0.2)
render_depth(focus = 0.7, focallength = 200)


# If you want more, this code below renders a video guaranteed to impress small children and executives. I borrowed this code from Joey Stanley who borrowed it from Morgan Wall.4
# Se você quiser mais, este código abaixo renderiza um vídeo que impressiona crianças pequenas e executivos. Peguei emprestado esse código de Riley King, que o pegou de  Joey Stanley, que o pegou de Morgan Wall.
# 3D Vowel Plots with Rayshader, http://joeystanley.com/blog/3d-vowel-plots-with-rayshader↩

library(av)

# Set up the camera position and angle
phivechalf = 30 + 60 * 1/(1 + exp(seq(-7, 20, length.out = 180)/2))
phivecfull = c(phivechalf, rev(phivechalf))
thetavec = 0 + 60 * sin(seq(0,359,length.out = 360) * pi/180)
zoomvec = 0.45 + 0.2 * 1/(1 + exp(seq(-5, 20, length.out = 180)))
zoomvecfull = c(zoomvec, rev(zoomvec))

# Actually render the video.
render_movie(filename = "hex_plot_fancy_2", type = "custom", 
            frames = 360,  phi = phivecfull, zoom = zoomvecfull, theta = thetavec)

render_movie(filename = "hex_plot_orbit", type = "orbit",
             phi = 45, theta = 60)
rgl::rgl.close()

#----------------------------------------------------------------------------------------
# Step 6: Push the parameters back through the model
#----------------------------------------------------------------------------------------

#After a lot of work we have finally identified the credible values for our model parameters. We now want to see what sort of predictions our posterior makes. Again, I’ll work with both the centered and un-centered data to try to understand the difference between the approaches. The first step in both cases is to create a sequence of time data to predict off of. For some reason I couldn’t get the predict() function in brms to cooperate so I wrote my own function to predict values. You enter a time value and the function makes a temperature prediction for every combination of mean and standard deviation derived from the parameters in the posterior distribution. Our goal will be to map this function over the sequence of predictor values we just set up.

#sequence of time data to predict off of.  Could use the same for both models but I created 2 for clarity
time_seq_tbl   <- tibble(pred_time   = seq(from = -15, to = 60, by = 1))
time_seq_tbl_2 <- tibble(pred_time_2 = seq(from = -15, to = 60, by = 1))

#function that takes a time value and makes a prediction using model_1 (un-centered) 
rk_predict <- 
  function(time_to_sim){
    rnorm(n = nrow(post_samplesM1_tbl),
          mean = post_samplesM1_tbl$b_Intercept + post_samplesM1_tbl$b_time*time_to_sim,
          sd = post_samplesM1_tbl$sigma
    )
  }

#function that takes a time value and makes a prediction using model_2 (centered)
rk_predict2 <- 
  function(time_to_sim){
    rnorm(n = nrow(post_samplesM2_tbl),
          mean = post_samplesM2_tbl$b_Intercept + post_samplesM2_tbl$b_time_c*time_to_sim,
          sd = post_samplesM2_tbl$sigma
    )
  }

#map the first prediction function over all values in the time sequence
#then calculate the .025 and .975 quantiles in anticipation of 95% prediction intervals
predicts_m1_tbl <- time_seq_tbl %>%
  mutate(preds_for_this_time = map(pred_time, rk_predict)) %>%
  mutate(percentile_2.5  = map_dbl(preds_for_this_time, ~quantile(., .025))) %>%
  mutate(percentile_97.5 = map_dbl(preds_for_this_time, ~quantile(., .975)))

#same for the 2nd prediction function
predicts_m2_tbl <- time_seq_tbl_2 %>%
  mutate(preds_for_this_time = map(pred_time_2, rk_predict2)) %>%
  mutate(percentile_2.5  = map_dbl(preds_for_this_time, ~quantile(., .025))) %>%
  mutate(percentile_97.5 = map_dbl(preds_for_this_time, ~quantile(., .975)))   

#visualize what is stored in the nested prediction cells (sanity check)
test_array <- predicts_m2_tbl[1, 2] %>% unnest(cols = c(preds_for_this_time))
test_array %>% 
  round(digits = 2) %>%
  head(5) %>%
  kable(align = rep("c", 1))


#vAnd now the grand finale - overlay the 95% prediction intervals on the original data along with the credible values of mean. We see there is no difference between the predictions made from centered data vs. un-centered.

big_enchilada <- 
  tibble(h=0) %>%
  ggplot() +
  geom_point(
    data = ablation_dta_tbl, aes(x = time, y = temp),
    colour = "#481567FF",
    size = 2,
    alpha = 0.6
  ) +
  geom_point(
    data = ablation_dta_tbl, aes(x = time_c, y = temp),
    colour = "#55C667FF",
    size = 2,
    alpha = 0.6
  ) +
  geom_abline(aes(intercept = b_Intercept, slope = b_time),
              data = post_samplesM1_tbl,
              alpha = 0.1, color = "gray50"
  ) +
  geom_abline(
    slope = mean(post_samplesM1_tbl$b_time),
    intercept = mean(post_samplesM1_tbl$b_Intercept),
    color = "blue", size = 1
  ) +
  geom_abline(aes(intercept = b_Intercept, slope = b_time_c),
              data = post_samplesM2_tbl,
              alpha = 0.1, color = "gray50"
  ) +
  geom_abline(
    slope = mean(post_samplesM2_tbl$b_time_c),
    intercept = mean(post_samplesM2_tbl$b_Intercept),
    color = "blue", size = 1
  ) +
  geom_ribbon(
    data = predicts_m1_tbl, aes(x = predicts_m1_tbl$pred_time, ymin = predicts_m1_tbl$percentile_2.5, ymax = predicts_m1_tbl$percentile_97.5), alpha = 0.25, fill = "pink", color = "black", size = .3
  ) +
  geom_ribbon(
    data = predicts_m2_tbl, aes(x = predicts_m2_tbl$pred_time_2, ymin = predicts_m2_tbl$percentile_2.5, ymax = predicts_m2_tbl$percentile_97.5), alpha = 0.4, fill = "pink", color = "black", size = .3
  ) +
  labs(
    title = "Regression Line Representing Mean of Slope",
    subtitle = "Centered and Un-Centered Predictor Data",
    x = "Time (s)",
    y = "Temperature (C)"
  ) +
  scale_x_continuous(limits = c(-10, 37), expand = c(0, 0)) +
  scale_y_continuous(limits = c(40, 120), expand = c(0, 0))

big_enchilada
