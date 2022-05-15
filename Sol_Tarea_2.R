#Soluciones Tarea 2 | Primavera 2022
# Creada por Miguel Negrete
# Revisada y editada por Arturo Aguilar

# Para ejecutar este script desde su computadora pueden guardar
# las bases de datos que emplea en un folder en su computadora
# y establecer dicho folder como el working directory cambiando
# la siguiente linea [ojo: utilicen las diagonales igual que 
# abajo sin importar si usan PC o Mac]. Alternativamente pueden
# crear un proyecto referenciado a dicho folder como vimos en 
# clase

#Cargamos los paquetes que utilizaremos en la tarea

library(sandwich)
library(cobalt) 
library(vtable)
library(tidyverse)
library(stargazer)
library(ggplot2)
library(car)
library(jtools) 
library(dplyr)
library(multcomp)

#Cargamos las base de datos frmgham

frmgham<-read.csv("frmgham.csv")
frmgham <- data.frame(frmgham)

# ====// Pregunta 1 \\====

# ====// Inciso (a)

#Calculamos los promedios de edad para los diferentes valores de CVD
(mean(frmgham[frmgham$CVD == 0, 'AGE']))
(mean(frmgham[frmgham$CVD == 1, 'AGE']))

#Realizamos los histogramas
label_CVD <- c("0" = "CVD = 0","1" = "CVD = 1")
ggplot(frmgham,aes(x = AGE , y = ..count../sum(..count..))) +
  geom_histogram(color = "black", alpha = 0.5, bins = 30) + 
  facet_grid(CVD ~ ., labeller = labeller(CVD = label_CVD)) +
  labs(x = "Edad", y = "Densidad")  +
  theme_classic()

ggsave("Hist_1a.jpeg",  width = 5.54, height = 4.95)

# ====// Inciso (b)

#Constriumos la variable dummy menor_50
frmgham <- frmgham %>%
  mutate( menor_50 = 
            ifelse (AGE < 50, 1,0 ))

#Realizamos nuestra estimación de mco
mco <- lm(CVD ~ menor_50, 
            data = frmgham)
es_mco<- sqrt(diag(vcovHC(mco, type = "HC1")))
stargazer(mco,type = "text",se=es_mco)

# ====// Inciso (c)

# Apartado (III)

#Calculamos todos los datos que necesitamos para estimar el esstadítico t
(X_menor <- mean(frmgham[frmgham$menor_50 == 1, 'CVD']))
(X_mayor <- mean(frmgham[frmgham$menor_50 == 0, 'CVD']))
(s_menor <- var(frmgham[frmgham$menor_50 == 1, 'CVD']))
(s_mayor <- var(frmgham[frmgham$menor_50 == 0, 'CVD']))
(n_menor <- length(frmgham[frmgham$menor_50 == 1, 'CVD']))
(n_mayor <- length(frmgham[frmgham$menor_50 == 0, 'CVD']))

#Estimamos t
(t_num <- (X_menor-X_mayor))
(t_den <- sqrt((s_menor/n_menor)+(s_mayor/n_mayor)))
(t<-(t_num/t_den))

#Varianza heterocedástica MCO
(var_het_mco<-(s_menor/n_menor)*((n_menor-1)/n_menor)+(s_mayor/n_mayor)*((n_mayor-1)/n_mayor))
(se_het_mco <- sqrt(var_het_mco))

# ====// Pregunta 2 \\====

#Creamos las variables que necesitaremos
frmgham <- frmgham %>% mutate(AGE_2= AGE^2,
                              SMKR_CIG = SMKR*CIGDAY)

#Especificación (1)
mco_1 <- lm(CVD ~ AGE + log(TOTCHOL) + HDLC + SYSBP + BPMEDS + SMKR + SMKR_CIG + DIABETES, 
            data = subset(frmgham, WOMEN == 1))
es_mco_1<- sqrt(diag(vcovHC(mco_1, type = "HC1")))

#Especificación (2)
mco_2 <- lm(CVD ~ AGE + log(TOTCHOL) + HDLC + SYSBP + BPMEDS + SMKR + SMKR_CIG + DIABETES, 
            data = subset(frmgham, WOMEN == 0))
es_mco_2<- sqrt(diag(vcovHC(mco_2, type = "HC1")))

#Especificación (3)
mco_3 <- lm(CVD ~ AGE + WOMEN + log(TOTCHOL) + HDLC + SYSBP + BPMEDS + SMKR + SMKR_CIG + DIABETES, 
            data = frmgham)
es_mco_3<- sqrt(diag(vcovHC(mco_3, type = "HC1")))

#Especificación (4)
mco_4 <- lm(CVD ~ AGE + AGE_2 + WOMEN + HDLC + SYSBP + PREVCHD + BPMEDS + SMKR + DIABETES + log(BMI), 
            data = frmgham)
es_mco_4<- sqrt(diag(vcovHC(mco_4, type = "HC1")))

#Tabla 
stargazer(mco_1,mco_2,mco_3,mco_4, 
          digits = 2,
          title = "Estimaciones de MCO",
          label = "tab:tabla_regresiones",
          omit.stat = c("f","adj.rsq","ser"), 
          font.size = "footnotesize", 
          column.sep.width = "1pt", 
          single.row = F,
          no.space = F,
          se=list(es_mco_1,es_mco_2,es_mco_3,es_mco_4), 
          out = "tabla_regs.tex")

# ====// Pregunta 3 \\====

#Definimos nuestra regresión residual y la incluimos en la Tabla de la pregunta 2
mco_5 <- lm(log(TOTCHOL) ~ AGE + AGE_2 + WOMEN + HDLC + SYSBP + PREVCHD + BPMEDS + SMKR + DIABETES + log(BMI), 
            data = frmgham)
es_mco_5<- sqrt(diag(vcovHC(mco_5, type = "HC1")))
stargazer(mco_1,mco_2,mco_3,mco_4,mco_5, 
          digits = 2,
          title = "Estimaciones de MCO",
          label = "tab:tabla_regresiones",
          omit.stat = c("f","adj.rsq","ser"), 
          font.size = "footnotesize", 
          column.sep.width = "1pt", 
          single.row = F,
          no.space = F,
          se=list(es_mco_1,es_mco_2,es_mco_3,es_mco_4,es_mco_5), 
          out = "tabla_regs.tex")

#Estimamos nuestra regresión larga
reg_larga <- lm(CVD ~ AGE + AGE_2 + WOMEN + log(TOTCHOL) + HDLC + SYSBP + PREVCHD + BPMEDS + SMKR + DIABETES + log(BMI), 
            data = frmgham)
es_reg_larga <- sqrt(diag(vcovHC(reg_larga, type = "HC1")))

#Guardamos los coeficientes que utilizaremos
Alpha <- mco_4$coefficients[4]
Beta <- reg_larga$coefficients[4]
Beta_K <- reg_larga$coefficients[5]
Gamma <- mco_5$coefficients[4]

#Calculamos el sesgo

(Beta_K*Gamma)
(Alpha-Beta)

#Coeficiente WOMEN regresión larga

(Alpha-Beta_K*Gamma)

# ====// Pregunta 5 \\====

# [ a ]
#Desviación estándar de SYSBP
(sd_sysbp <- sd(frmgham$SYSBP, na.rm = T))

#Datos promedio 
(CVD_promedio <- mean(frmgham$CVD, na.rm = T))
(SYSBP_promedio <- mean(frmgham$SYSBP, na.rm = T))

#Cambio de CVD asociado a la columna 3
(Ch_CVD <- sd_sysbp*(mco_3$coef["SYSBP"]))
(Ch_percent_CVD <- 100*Ch_CVD/CVD_promedio)

#Error Estándar
l <- c(0,0,0,0,0,sd_sysbp,0,0,0,0) 
(mvcov_mco_3 <- vcovHC(mco_3, type = "HC1"))
(se_5a <- sqrt(t(l) %*% mvcov_mco_3 %*% l))

#Intervalo de confianza
(IC_5a <- (100/CVD_promedio)*c(t(l) %*% mco_3$coef - 1.65*se_5a,
                               t(l) %*% mco_3$coef + 1.65*se_5a))  


# ====// Pregunta 6 \\====

# ====// Inciso (a)

#Cargamos los datos médicos de Oswald
Os_data_1_2 <- c(1,43,log(240),35,160,0,1,3,0)
Os_data_3 <-c(1,43,0,log(240),35,160,0,1,3,0)
Os_data_4 <-c(1,43,43**2,0,35,160,0,0,1,0,log(36))

#Evaluamos las especificaciones
(eval_1<-t(Os_data_1_2)%*%mco_1$coefficients)
(eval_2<-t(Os_data_1_2)%*%mco_2$coefficients)
(eval_3<-t(Os_data_3)%*%mco_3$coefficients)
(eval_4<-t(Os_data_4)%*%mco_4$coefficients)

# ====// Inciso (b)

#Construimos y evaluamos los datos clínicos de Oswald
Os_today <- Os_data_4
Os_tomorrow <- c(1,48,48**2,0,40,130,0,1,1,0,log(28))
(eval_today<-t(Os_today)%*%mco_4$coefficients)
(eval_tomorrow<-t(Os_tomorrow)%*%mco_4$coefficients)

#Diferencia
diff<-Os_tomorrow-Os_today

#Predicción diferencia
(diff_p<-t(diff)%*%mco_4$coefficients)

#Evaluación de la hipótesis
hyp_1 <-linearHypothesis(mco_4,diff,white.adjust = 'hc1', rhs=-0.02)
(t_1<-sqrt(hyp_1$F[2]))
(pval_1<-2*(pnorm(t_1,0,1,lower.tail = F)))

# ====// Inciso (c)

# Utilizaremos la función glht para calcular el intervalo de confianza

#Para poder utilizar la función, renombramos la variable log(BMI) de la estimación 4 a log_BMI
log_BMI <- log(frmgham$BMI)
mco_4_bis <- lm(CVD ~ AGE + AGE_2 + WOMEN + HDLC + SYSBP + PREVCHD + BPMEDS + SMKR + DIABETES + log_BMI, 
            data = frmgham)
es_mco_4_bis<- sqrt(diag(vcovHC(mco_4_bis, type = "HC1")))

# Guardamos la combinacion lineal como un objeto llamado "lincom".
lincom <- glht(mco_4_bis, linfct = c("5*AGE + 455*AGE_2 + 5*HDLC - 30*SYSBP + BPMEDS - 0.25*log_BMI = 0.02"), vcov=vcovHC(mco_4, type="HC1"))

# Obtenemos el intervalo de confianza.
confint(lincom)

# ====// Pregunta 7 \\====

# ====// Inciso (a)

linearHypothesis(mco_3,c("SMKR=0","SMKR_CIG=0"),white.adjust="hc1")

# ====// Inciso (b)

#Obtenemos la muestra de la especificación (3)
frmgham_bis <- frmgham %>% 
      drop_na(RANDID,CVD, AGE, WOMEN, TOTCHOL, HDLC, SYSBP, BPMEDS, SMKR, SMKR_CIG, DIABETES)

#Creamos el vector donde guardaremos los resultados de Bootstrap
betas_diff <-c()

for (i in 1:1000) {
  
  muestra_sim <- sample_n(frmgham_bis, size = nrow(frmgham_bis), replace = T)
  
  #Especificación mujeres
  mc_1 <- lm(CVD ~ AGE + log(TOTCHOL) + HDLC + SYSBP + BPMEDS + SMKR + SMKR_CIG + DIABETES, 
              data = subset(muestra_sim, WOMEN == 1))
  es_mc_1<- sqrt(diag(vcovHC(mc_1, type = "HC1")))
  
  #Especificación hombres
  mc_2 <- lm(CVD ~ AGE + log(TOTCHOL) + HDLC + SYSBP + BPMEDS + SMKR + SMKR_CIG + DIABETES, 
              data = subset(muestra_sim, WOMEN == 0))
  es_mc_2<- sqrt(diag(vcovHC(mc_2, type = "HC1")))
  
  #Calculamos la diferencia y la guardamos
  dif_sim <- mc_1$coefficients[2]-mc_2$coefficients[2]
  betas_diff <-c(betas_diff,dif_sim)
  
}

difs_df <- as.data.frame(betas_diff)

#Histograma de las diferencias
ggplot() +
  geom_histogram(data = difs_df, aes(betas_diff, ..count../sum(..count..)), color = "#2b8cbe", fill = "#a6bddb") +
  labs(x = "Diferencia en coeficientes para la variable AGE entre hombres y mujeres", y = "Densidad") +
  theme_minimal()

ggsave("hist_1c_7.jpeg",  width = 5.54, height = 4.95)

#Prueba de hipótesis
(m_1c <- mean(difs_df$betas_diff, na.rm = T))
var_1c <- var(difs_df$betas_diff, na.rm = T)

#Obtenemos el estadístico t y el valor-p
(t_1c_iii <- m_1c/sqrt(var_1c))
(p_value_1c_iii <- 2*(1-pnorm(abs(t_1c_iii),0,1)))

