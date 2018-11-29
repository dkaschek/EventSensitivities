library(dMod)
library(dplyr)
library(ggplot2)
library(cowplot)

## Model text input ---------------------------------------------------------

# Reactions
reactions <- eqnlist() %>% 
  addReaction("0", "Input", "0") %>% 
  addReaction("0", "Ad", "Favail*Input") %>% 
  addReaction("Ad", "Ac", "KA*Ad") %>%
  addReaction("Ac", "0", "CL/V*Ac")

# Events
events <- eventlist() %>% 
  addEvent(var = "Input", time = "tlag", value = "rate", method = "replace") %>% 
  addEvent(var = "Input", time = "tend", value = 0, method = "replace")

# Variables/Observations
variables <- eqnvec(
  Effect = "E0*(1 - (Ac/V)*IMAX/(IC50 + (Ac/V)))"
)

# Parameters
parameters <- c(
  KA = 1,
  CL = 6,
  V = 60,
  IMAX = 1,
  IC50 = 1,
  E0 = 15,
  Favail = 1
)

## Model functional output --------------------------------------------------------

# ODE model
model <- reactions %>% 
  odemodel(events = events, modelname = "PK", compile = FALSE)

# Prediction function
x <- model %>% Xs()

# Observation function
g <- variables %>% 
  Y(f = x, modelname = "variables", compile = FALSE)

# Parameter transformation function
p <- eqnvec() %>% 
  define("x~x", x = c(names(parameters), "tlag", "rate", "tend")) %>% 
  define("x~0", x = c("Ad", "Ac", "Input")) %>% 
  insert("x~exp(x)", x = names(parameters)) %>% 
  define("tend~tlag+tinf") %>% 
  define("rate~dose/tinf") %>% 
  P(condition = "Standard", modelname = "parameters", compile = FALSE)

# Compile into one .so file and remove temporary files
compile(g, x, p, output = "PKPD")
tmpfiles <- list.files(pattern = c("\\.[oc]"))
unlink(tmpfiles)

## Model simulation ---------------------------------------------------------------

# Set parameter values for simulation
pars <- c(log(parameters), tlag = 10, dose = 200, tinf = 10)

# Simulation times
times <- seq(0, 60, .1)

# Simulate and restrict output to Input and Effect
out1 <- (g*x*p)(times, pars) %>% 
  as.data.frame() %>% 
  filter(name %in% c("Input", "Effect"))

# Simulate, get derivatives (wrt outer parameters) and restrict to subset
out2 <- (g*x*p)(times, pars) %>% 
  getDerivs() %>% 
  as.data.frame() %>% 
  filter(name %in% c("Effect.tlag", "Effect.dose", "Effect.tinf")) %>% 
  mutate(name = c("Effect.tlag" = "dEffect / dtlag", "Effect.dose" = "dEffect / ddose", "Effect.tinf" = "dEffect / dtinf")[as.character(name)])


## Generate plots -----------------------------------------------------------------

# Plot model simulation
plot1 <- out1 %>% 
  ggplot(aes(x = time, y = value, color = name, lty = name)) +
  geom_line() +
  # Set colors and linetype
  scale_color_manual(name = NULL, values = c("black", "firebrick2")) +
  scale_linetype_manual(name = NULL, values = c(1, 2)) +
  # Annotate tinf
  annotate("segment", x = pars["tlag"], xend = pars["tlag"] + pars["tinf"], y = 22, yend = 22, arrow = arrow(ends = "both", length = unit(2, "mm"))) +
  annotate("text", x = pars["tlag"] + 0.5*pars["tinf"], y = 24, label = "tinf", size = 2.5) +
  # Annotate tlag
  annotate("segment", x = 0, xend = pars["tlag"], y = -2, yend = -2, arrow = arrow(ends = "both", length = unit(2, "mm"))) +
  annotate("text", x = 0.5*pars["tlag"], y = -4, label = "tlag", size = 2.5) +
  # Annotate dose
  annotate("segment", x = pars["tlag"] - 3, xend = pars["tlag"] - 3, y = 0, yend = 20, arrow = arrow(ends = "both", length = unit(2, "mm"))) +
  annotate("text", x = pars["tlag"] - 6, y = 10, label = "dose/tinf", angle = 90, size = 2.5) +
  # Set theme
  theme_dMod(base_size = 8) +
  theme(panel.grid = element_blank(), 
        legend.position = c(.95, .15), 
        legend.justification = c(1, 0), 
        legend.background = element_blank(), 
        legend.key = element_blank(), 
        legend.key.height = unit(3, "mm"))



plot1


# Plot Effect sensitivities
plot2 <- out2 %>% 
  ggplot(aes(x = time, y = value)) + 
  facet_wrap(~name, ncol = 1, scales = "free_y", labeller = label_parsed) +
  geom_hline(yintercept = 0, color = "darkgray", lty = 2) +
  geom_line() +
  geom_vline(xintercept = c(10, 20), lty = 2, color = "firebrick2") +
  # Set theme
  theme_dMod(base_size = 8) +
  theme(panel.grid = element_blank(), legend.position = c(.9, .9))
  
plot2

# Create grid plot using cowplot
plot12 <- plot_grid(plot1, plot2, align = "v", ncol = 1, labels = c("A", "B"), rel_heights = c(.35, .65), label_size = 9)

# Export to pdf using Latin Modern Roman font for consistency with figure caption font
ggplot2::ggsave("example.pdf", plot12, width = 7.5, height = 13, units = "cm", device = cairo_pdf, family = "Latin Modern Roman")
