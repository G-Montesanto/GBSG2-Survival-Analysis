# ANALISI DI SOPRAVVIVENZA - DATASET GBSG2
#
# Obiettivo:
# Analizzare il tempo libero da recidiva o morte attraverso:
# 1. Analisi descrittiva ed esplorativa
# 2. Curve di Kaplan-Meier e test log-rank
# 3. Valutazione della forma funzionale delle covariate
# 4. Modelli di Cox multivariati
# 5. Verifica dell'ipotesi di rischi proporzionali
# 6. Confronto tra modelli lineari, logaritmici e spline
# 7. Analisi delle osservazioni influenti
# 8. Analisi di sensibilità
#
# Dataset: TH.data::GBSG2
#
# Nota metodologica:
# L'evento rappresenta l'endpoint composito registrato nel dataset.
# I risultati descrivono associazioni statistiche e non dimostrano
# automaticamente relazioni causali.
# Le analisi di sensibilità servono a valutare la stabilità delle
# stime, non a giustificare automaticamente l'esclusione di pazienti.


# 1. CARICAMENTO DEI PACCHETTI

library(TH.data)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(splines)
library(car)


# 2. CARICAMENTO E CONTROLLO INIZIALE DEI DATI

data("GBSG2")

# Esaminiamo struttura, tipi delle variabili e prime osservazioni.
str(GBSG2)
head(GBSG2)

# Verifichiamo la presenza di dati mancanti.
colSums(is.na(GBSG2))

# Distribuzione delle variabili categoriche e dell'indicatore di evento.
table(GBSG2$horTh)
table(GBSG2$menostat)
table(GBSG2$tgrade)
table(GBSG2$cens, useNA = "ifany")

# Distribuzione del tempo osservato, espresso originariamente in giorni.
summary(GBSG2$time)

# Controlliamo alcune righe per verificare la codifica di tempo ed evento.
head(GBSG2[, c("time", "cens")], 20)

# Il dataset contiene 686 pazienti e non presenta valori mancanti.
# Prima di interpretare i risultati, è comunque importante consultare
# la documentazione del dataset per verificare il significato esatto
# dell'indicatore cens e dell'endpoint analizzato.


# 3. PREPARAZIONE DEL DATASET

# Convertiamo il tempo da giorni a mesi utilizzando 30,44 giorni
# come durata media convenzionale di un mese.
#
# Manteniamo la codifica originale di cens nella variabile event.
# In survival::Surv(), 1 indica l'evento e 0 indica la censura.

d <- GBSG2 |>
  mutate(
    months = time / 30.44,
    event = cens
  )

head(d)

# Numero e percentuale di eventi e osservazioni censurate.
table(d$event)
round(prop.table(table(d$event)) * 100, 2)

# Distribuzione del tempo osservato in mesi.
summary(d$months)

# Il numero di eventi è fondamentale per comprendere quanta
# informazione è disponibile per stimare i modelli di sopravvivenza.


# 4. ANALISI DESCRITTIVA DEL CAMPIONE

# Statistiche descrittive delle variabili numeriche.
summary(
  d[, c(
    "age", "tsize", "pnodes",
    "progrec", "estrec", "months"
  )]
)

# Calcoliamo media e deviazione standard.
d |>
  summarise(
    age_mean = mean(age),
    age_sd = sd(age),
    tsize_mean = mean(tsize),
    tsize_sd = sd(tsize),
    pnodes_mean = mean(pnodes),
    pnodes_sd = sd(pnodes),
    progrec_mean = mean(progrec),
    progrec_sd = sd(progrec),
    estrec_mean = mean(estrec),
    estrec_sd = sd(estrec),
    months_mean = mean(months),
    months_sd = sd(months)
  )

# Frequenze assolute delle variabili categoriche.
table(d$horTh)
table(d$menostat)
table(d$tgrade)

# Percentuali delle variabili categoriche.
round(prop.table(table(d$horTh)) * 100, 2)
round(prop.table(table(d$menostat)) * 100, 2)
round(prop.table(table(d$tgrade)) * 100, 2)

# Le statistiche descrittive permettono di conoscere il campione
# prima di procedere alla modellazione.
# Media e deviazione standard vanno interpretate insieme a mediana,
# quantili e grafici, soprattutto in presenza di distribuzioni asimmetriche.


# 5. DISTRIBUZIONE DELLE VARIABILI NUMERICHE

# Gli istogrammi permettono di individuare:
# - asimmetria delle distribuzioni;
# - concentrazione delle osservazioni;
# - possibili valori estremi;
# - variabili per cui può essere utile esplorare trasformazioni.

ggplot(d, aes(x = age)) +
  geom_histogram(bins = 20) +
  labs(
    title = "Distribuzione dell'età",
    x = "Età (anni)",
    y = "Frequenza"
  ) +
  theme_minimal()

ggplot(d, aes(x = tsize)) +
  geom_histogram(bins = 20) +
  labs(
    title = "Distribuzione della dimensione del tumore",
    x = "Dimensione del tumore",
    y = "Frequenza"
  ) +
  theme_minimal()

ggplot(d, aes(x = pnodes)) +
  geom_histogram(bins = 25) +
  labs(
    title = "Distribuzione dei linfonodi positivi",
    x = "Numero di linfonodi positivi",
    y = "Frequenza"
  ) +
  theme_minimal()

ggplot(d, aes(x = progrec)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribuzione del recettore del progesterone",
    x = "Progesterone receptor",
    y = "Frequenza"
  ) +
  theme_minimal()

ggplot(d, aes(x = estrec)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribuzione del recettore degli estrogeni",
    x = "Estrogen receptor",
    y = "Frequenza"
  ) +
  theme_minimal()

# In particolare, pnodes, progrec ed estrec presentano distribuzioni
# asimmetriche a destra.
# Per queste variabili valuteremo trasformazioni logaritmiche e forme
# funzionali alternative nei modelli di Cox.
#
# Attenzione: la forma della distribuzione di una covariata non determina
# da sola la forma della sua relazione con il log-hazard.


# 6. QUANTILI DELLE VARIABILI ASIMMETRICHE

# I quantili descrivono la distribuzione e rendono visibile
# quanto siano concentrate le osservazioni nella parte bassa
# e quanto sia lunga la coda destra.

quantile(
  d$pnodes,
  probs = c(0, .25, .5, .75, .9, .95, .99, 1),
  na.rm = TRUE
)

quantile(
  d$progrec,
  probs = c(0, .25, .5, .75, .9, .95, .99, 1),
  na.rm = TRUE
)

quantile(
  d$estrec,
  probs = c(0, .25, .5, .75, .9, .95, .99, 1),
  na.rm = TRUE
)

# Per pnodes la mediana è inferiore alla media e sono presenti
# pochi pazienti con un numero molto elevato di linfonodi positivi.
# Questa caratteristica diventerà importante nella valutazione
# della forma funzionale e nell'analisi delle osservazioni influenti.


# 7. CORRELAZIONE TRA PROGREC ED ESTREC

# La correlazione di Spearman valuta l'associazione monotona tra
# due variabili senza richiedere una relazione lineare.
#
# È utile per descrivere l'associazione tra i due indicatori,
# ma non misura un rapporto causale.

cor.test(
  d$progrec,
  d$estrec,
  method = "spearman"
)

# Rappresentiamo la relazione usando scale logaritmiche.
# La trasformazione log1p equivale a log(1 + x) e permette
# di visualizzare anche valori uguali a zero.

ggplot(d, aes(x = progrec, y = estrec)) +
  geom_point(alpha = 0.4) +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "Relazione tra recettori del progesterone ed estrogeni",
    x = "Progesterone receptor (scala log1p)",
    y = "Estrogen receptor (scala log1p)"
  ) +
  theme_minimal()

# Il risultato precedentemente ottenuto indica una correlazione
# positiva di Spearman di circa 0,60.
# I due indicatori tendono quindi ad aumentare insieme.
#
# Nel modello multivariato li includeremo entrambi, valutando
# l'associazione di ciascuno con l'hazard a parità delle altre covariate.
# La correlazione non implica necessariamente che uno dei due debba
# essere eliminato, ma è utile tenerla presente nell'interpretazione.


# 8. TRASFORMAZIONI LOGARITMICHE

# Creiamo le trasformazioni log1p delle variabili positive.
# La trasformazione può ridurre l'asimmetria e rappresentare
# una relazione non lineare con la variabile originale.
#
# Non garantisce, da sola, che la relazione con il log-hazard
# sia correttamente specificata.

d <- d |>
  mutate(
    log_pnodes = log1p(pnodes),
    log_progrec = log1p(progrec),
    log_estrec = log1p(estrec),
    log_tsize = log1p(tsize)
  )

# Distribuzione delle variabili trasformate.

ggplot(d, aes(x = log_pnodes)) +
  geom_histogram(bins = 20) +
  labs(
    title = "Distribuzione di log1p(pnodes)",
    x = "log1p(pnodes)",
    y = "Frequenza"
  ) +
  theme_minimal()

ggplot(d, aes(x = log_progrec)) +
  geom_histogram(bins = 25) +
  labs(
    title = "Distribuzione di log1p(progrec)",
    x = "log1p(progrec)",
    y = "Frequenza"
  ) +
  theme_minimal()

ggplot(d, aes(x = log_estrec)) +
  geom_histogram(bins = 25) +
  labs(
    title = "Distribuzione di log1p(estrec)",
    x = "log1p(estrec)",
    y = "Frequenza"
  ) +
  theme_minimal()

ggplot(d, aes(x = log_tsize)) +
  geom_histogram(bins = 20) +
  labs(
    title = "Distribuzione di log1p(tsize)",
    x = "log1p(tsize)",
    y = "Frequenza"
  ) +
  theme_minimal()

# Confrontiamo media e mediana prima e dopo la trasformazione.

d |>
  summarise(
    pnodes_mean = mean(pnodes),
    pnodes_median = median(pnodes),
    log_pnodes_mean = mean(log_pnodes),
    log_pnodes_median = median(log_pnodes),
    
    progrec_mean = mean(progrec),
    progrec_median = median(progrec),
    log_progrec_mean = mean(log_progrec),
    log_progrec_median = median(log_progrec),
    
    estrec_mean = mean(estrec),
    estrec_median = median(estrec),
    log_estrec_mean = mean(log_estrec),
    log_estrec_median = median(log_estrec),
    
    tsize_mean = mean(tsize),
    tsize_median = median(tsize),
    log_tsize_mean = mean(log_tsize),
    log_tsize_median = median(log_tsize)
  )

# Le trasformazioni modificano la scala di misura e possono ridurre
# l'influenza numerica della coda destra.
# La scelta finale della forma funzionale verrà però basata anche
# sul confronto dei modelli di sopravvivenza.


# 9. KAPLAN-MEIER COMPLESSIVA

# La funzione Surv combina il tempo osservato e l'indicatore dell'evento.
# I soggetti censurati contribuiscono all'informazione disponibile
# fino al momento in cui vengono censurati.

surv_obj <- Surv(
  time = d$months,
  event = d$event
)

# Stima Kaplan-Meier della probabilità di rimanere liberi dall'evento.
km_fit <- survfit(
  surv_obj ~ 1,
  data = d
)

summary(km_fit)

# Grafico con intervalli di confidenza e tabella dei soggetti a rischio.

ggsurvplot(
  km_fit,
  data = d,
  conf.int = TRUE,
  risk.table = TRUE,
  xlab = "Tempo (mesi)",
  ylab = "Probabilità di rimanere liberi dall'evento",
  title = "Kaplan-Meier complessiva",
  surv.median.line = "hv"
)

# La curva descrive come cambia nel tempo la probabilità stimata
# di non avere ancora sperimentato l'evento.
#
# La mediana di sopravvivenza è il tempo in cui la curva scende a 0,5,
# purché tale livello venga raggiunto.
# La stima è descrittiva e non identifica, da sola, gli effetti
# delle singole covariate.


# 10. KAPLAN-MEIER PER TERAPIA ORMONALE

# Confrontiamo le curve dei pazienti con e senza terapia ormonale.

km_horTh <- survfit(
  Surv(months, event) ~ horTh,
  data = d
)

ggsurvplot(
  km_horTh,
  data = d,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Tempo (mesi)",
  ylab = "Probabilità di rimanere liberi dall'evento",
  title = "Kaplan-Meier per terapia ormonale",
  legend.title = "Terapia ormonale",
  legend.labs = levels(d$horTh)
)

# Test log-rank per verificare l'uguaglianza delle curve.

survdiff(
  Surv(months, event) ~ horTh,
  data = d
)

# Il log-rank verifica se le curve di sopravvivenza differiscono
# nel complesso.
# Un p-value piccolo fornisce evidenza di una differenza tra i gruppi,
# ma non dimostra che la terapia sia la causa della differenza.
# Per un confronto aggiustato per le altre caratteristiche dei pazienti
# useremo successivamente il modello di Cox multivariato.


# 11. KAPLAN-MEIER PER GRADO TUMORALE

# Confrontiamo la sopravvivenza tra i diversi gradi tumorali.

km_grade <- survfit(
  Surv(months, event) ~ tgrade,
  data = d
)

ggsurvplot(
  km_grade,
  data = d,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Tempo (mesi)",
  ylab = "Probabilità di rimanere liberi dall'evento",
  title = "Kaplan-Meier per grado tumorale",
  legend.title = "Grado tumorale",
  legend.labs = levels(d$tgrade)
)

survdiff(
  Surv(months, event) ~ tgrade,
  data = d
)

# Curve differenti indicano che la sopravvivenza osservata varia
# tra i gruppi di grado tumorale.
# Il test log-rank non corregge per le altre variabili e non indica
# automaticamente quali coppie di gruppi differiscano tra loro.


# 12. KAPLAN-MEIER PER NUMERO DI LINFONODI POSITIVI

# Creiamo tre gruppi per rendere più leggibile il confronto descrittivo.
# Il raggruppamento comporta perdita di informazione: nei modelli
# multivariati pnodes verrà analizzata anche come variabile continua.

d <- d |>
  mutate(
    pnodes_group = cut(
      pnodes,
      breaks = c(0, 3, 7, Inf),
      labels = c("1-3", "4-7", ">7")
    )
  )

table(d$pnodes_group)

km_pnodes <- survfit(
  Surv(months, event) ~ pnodes_group,
  data = d
)

ggsurvplot(
  km_pnodes,
  data = d,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Tempo (mesi)",
  ylab = "Probabilità di rimanere liberi dall'evento",
  title = "Kaplan-Meier per numero di linfonodi positivi",
  legend.title = "Linfonodi positivi"
)

# Test log-rank per confrontare le tre curve.

survdiff(
  Surv(months, event) ~ pnodes_group,
  data = d
)

# Nell'analisi precedente il test log-rank ha mostrato differenze
# marcate tra i gruppi.
# Il gruppo con più linfonodi positivi presenta una prognosi
# complessivamente meno favorevole.
#
# Questo risultato è coerente con l'ipotesi che un numero maggiore
# di linfonodi positivi sia associato a un rischio più elevato.
# Per quantificare l'associazione e aggiustare per altre covariate
# passiamo ai modelli di Cox.


# 13. COX UNIVARIATO: PNODES LINEARE

# Il modello di Cox descrive l'hazard istantaneo in funzione
# delle covariate, senza richiedere di specificare una forma
# parametrica per la funzione di rischio di base.
#
# In questa prima analisi consideriamo pnodes come effetto lineare.

cox_pnodes <- coxph(
  Surv(months, event) ~ pnodes,
  data = d
)

summary(cox_pnodes)

# L'hazard ratio associato a pnodes rappresenta il rapporto tra
# gli hazard per un aumento di un'unità del numero di linfonodi,
# assumendo una relazione lineare tra pnodes e log-hazard.

# Dai risultati precedenti:
# HR circa 1,06 per ogni linfonodo positivo aggiuntivo.
# L'associazione è statisticamente significativa.
#
# Il modello lineare è una prima approssimazione: non garantisce
# che l'effetto di ogni linfonodo aggiuntivo sia costante
# lungo tutta la distribuzione di pnodes.


# Residui di Martingala per esplorare la forma funzionale.

d$martingale_pnodes <- residuals(
  cox_pnodes,
  type = "martingale"
)

ggplot(
  d,
  aes(x = pnodes, y = martingale_pnodes)
) +
  geom_point(alpha = 0.35) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "Numero di linfonodi positivi",
    y = "Residui di Martingala",
    title = "Esplorazione della forma funzionale di pnodes"
  ) +
  theme_minimal()

# Un andamento sistematico della curva smussata può suggerire
# che una relazione lineare sia insufficiente.
# I residui sono uno strumento esplorativo: le conclusioni vanno
# confrontate con modelli alternativi e relativi criteri di adattamento.


# 14. CONFRONTO UNIVARIATO DELLA FORMA FUNZIONALE DI PNODES

# Confrontiamo tre specificazioni:
# 1. pnodes lineare;
# 2. log1p(pnodes);
# 3. spline naturale con 4 gradi di libertà.
#
# L'AIC bilancia adattamento e complessità: tra modelli confrontabili,
# un valore più basso indica il modello preferibile secondo questo criterio.

cox_linear <- coxph(
  Surv(months, event) ~ pnodes,
  data = d
)

cox_log <- coxph(
  Surv(months, event) ~ log_pnodes,
  data = d
)

cox_spline <- coxph(
  Surv(months, event) ~ ns(pnodes, df = 4),
  data = d
)

AIC(
  cox_linear,
  cox_log,
  cox_spline
)

# Il test del rapporto di verosimiglianza confronta il modello lineare
# con quello spline.
# Poiché il modello lineare è annidato nella spline, il test verifica
# se i termini aggiuntivi della spline migliorano significativamente
# l'adattamento.

anova(
  cox_linear,
  cox_spline,
  test = "LRT"
)

summary(cox_linear)
summary(cox_log)
summary(cox_spline)

# Il confronto univariato è esplorativo.
# La forma funzionale scelta verrà rivalutata nei modelli multivariati,
# mantenendo costanti le altre covariate e la stratificazione.


# 15. FORMA FUNZIONALE DELLE ALTRE VARIABILI CONTINUE

# Confrontiamo modelli univariati con covariate originali
# e trasformate, dove previsto.

cox_age <- coxph(
  Surv(months, event) ~ age,
  data = d
)

cox_tsize <- coxph(
  Surv(months, event) ~ tsize,
  data = d
)

cox_progrec <- coxph(
  Surv(months, event) ~ progrec,
  data = d
)

cox_estrec <- coxph(
  Surv(months, event) ~ estrec,
  data = d
)

cox_tsize_log <- coxph(
  Surv(months, event) ~ log_tsize,
  data = d
)

cox_progrec_log <- coxph(
  Surv(months, event) ~ log_progrec,
  data = d
)

cox_estrec_log <- coxph(
  Surv(months, event) ~ log_estrec,
  data = d
)

# Confronto degli AIC.
AIC(
  cox_age,
  cox_tsize,
  cox_tsize_log,
  cox_progrec,
  cox_progrec_log,
  cox_estrec,
  cox_estrec_log
)

summary(cox_age)
summary(cox_tsize)
summary(cox_tsize_log)
summary(cox_progrec)
summary(cox_progrec_log)
summary(cox_estrec)
summary(cox_estrec_log)

# Un AIC inferiore suggerisce una specificazione preferibile
# tra i modelli confrontati per la stessa variabile.
# I confronti univariati non sono sufficienti per determinare
# da soli la forma funzionale definitiva del modello multivariato.


# 16. MODELLO DI COX MULTIVARIATO PRELIMINARE

# Includiamo simultaneamente le covariate selezionate.
# La spline su pnodes consente di rappresentare un'associazione
# potenzialmente non lineare.
#
# In questa prima versione menostat e tgrade sono incluse come
# covariate ordinarie, così da poter esaminare la loro compatibilità
# con l'ipotesi di rischi proporzionali.

cox_full <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    ns(pnodes, df = 4) +
    log_progrec +
    log_estrec +
    horTh +
    menostat +
    tgrade,
  data = d
)

summary(cox_full)

# Nel modello multivariato ogni coefficiente descrive l'associazione
# tra la covariata corrispondente e l'hazard, a parità delle altre
# covariate incluse nel modello.
#
# Per una covariata continua lineare, exp(coef) è l'HR associato
# a un incremento unitario.
# Per le variabili trasformate, l'incremento va interpretato sulla
# scala trasformata.
# Per una spline, i singoli coefficienti della base non rappresentano
# direttamente HR clinicamente interpretabili.


# 17. VERIFICA DELL'IPOTESI DI RISCHI PROPORZIONALI

# Il test di Schoenfeld valuta se l'effetto delle covariate
# è compatibile con l'ipotesi di rischi proporzionali.
#
# L'ipotesi PH richiede che gli hazard ratio siano costanti nel tempo.

ph_test <- cox.zph(cox_full)

print(ph_test)

# Il grafico globale e quelli per singola covariata aiutano
# a individuare eventuali andamenti temporali dei residui.

plot(ph_test)

# Grafici selezionati per alcune covariate.
# I nomi devono corrispondere a quelli riportati da cox.zph.

par(mfrow = c(2, 3))

plot(ph_test, var = "age")
plot(ph_test, var = "log_progrec")
plot(ph_test, var = "log_estrec")
plot(ph_test, var = "menostat")
plot(ph_test, var = "tgrade")

par(mfrow = c(1, 1))

# Un p-value piccolo può indicare una violazione dell'ipotesi PH.
# È importante considerare sia i test sia i grafici, perché i test
# possono essere sensibili alla numerosità campionaria e non descrivono
# da soli la rilevanza pratica di un'eventuale violazione.
#
# Se una covariata categorica presenta una violazione rilevante,
# una possibile soluzione è stratificare il modello per quella variabile.


# 18. MODELLO DI COX STRATIFICATO

# Stratificando per menostat e tgrade permettiamo alla funzione
# di rischio di base di essere diversa nei vari strati.
#
# Per le variabili inserite in strata() non viene stimato un singolo
# coefficiente HR, perché il loro effetto viene gestito attraverso
# funzioni di rischio di base distinte.

cox_strat <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    ns(pnodes, df = 4) +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d
)

summary(cox_strat)

# Ricontrolliamo l'ipotesi PH sulle covariate rimaste nel modello.

ph_strat <- cox.zph(cox_strat)

print(ph_strat)

par(mfrow = c(1, 2))

plot(
  ph_strat,
  var = "age",
  main = "Residui di Schoenfeld: age"
)

plot(
  ph_strat,
  var = "log_estrec",
  main = "Residui di Schoenfeld: log_estrec"
)

par(mfrow = c(1, 1))

# Il modello stratificato è una soluzione possibile quando le
# variabili categoriche non rispettano l'ipotesi PH.
# La stratificazione non dimostra che il problema sia risolto per
# ogni altra covariata: è necessario controllare nuovamente i residui.
#
# La scelta di stratificare deve inoltre essere motivata dai risultati
# diagnostici e dalla domanda scientifica, perché impedisce di stimare
# un HR diretto per le variabili utilizzate come strati.


# 19. MODELLO SPLINE CON NODI PREFISSATI

# Definiamo una spline naturale con nodi interni 3, 7 e 15
# e nodi di bordo 1 e 51.
#
# La posizione dei nodi determina la flessibilità della curva.
# È opportuno verificare la stabilità dei risultati e non interpretare
# eccessivamente la forma stimata nelle regioni con pochi pazienti.

cox_spline_fixed <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    ns(
      pnodes,
      knots = c(3, 7, 15),
      Boundary.knots = c(1, 51)
    ) +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d,
  x = TRUE
)

summary(cox_spline_fixed)

# Verifichiamo nuovamente l'ipotesi di rischi proporzionali.

ph_fixed <- cox.zph(cox_spline_fixed)

print(ph_fixed)

par(mfrow = c(1, 2))

plot(
  ph_fixed,
  var = "age",
  main = "Residui di Schoenfeld: age"
)

plot(
  ph_fixed,
  var = "log_estrec",
  main = "Residui di Schoenfeld: log_estrec"
)

par(mfrow = c(1, 1))

# Questo modello rappresenta una specificazione flessibile
# dell'associazione tra pnodes e hazard.
# La sua adeguatezza verrà confrontata con modelli alternativi
# usando AIC, test del rapporto di verosimiglianza e C-index.


# 20. CURVA HAZARD RATIO PER PNODES NEL MODELLO SPLINE

# Calcoliamo l'HR associato a diversi valori di pnodes rispetto
# al riferimento di 3 linfonodi positivi.
#
# Manteniamo le altre covariate costanti, usando valori rappresentativi.
# Per le variabili categoriche scegliamo livelli espliciti.
#
# Le bande di confidenza vengono calcolate dalla differenza tra
# i predittori lineari e tengono conto della covarianza tra stime.

new_pnodes <- data.frame(
  pnodes = seq(1, 30, by = 0.1)
)

new_pnodes$age <- median(d$age)
new_pnodes$log_tsize <- median(d$log_tsize)
new_pnodes$log_progrec <- median(d$log_progrec)
new_pnodes$log_estrec <- median(d$log_estrec)

new_pnodes$horTh <- factor(
  "no",
  levels = levels(d$horTh)
)

new_pnodes$menostat <- factor(
  "Pre",
  levels = levels(d$menostat)
)

new_pnodes$tgrade <- factor(
  "II",
  levels = levels(d$tgrade),
  ordered = is.ordered(d$tgrade)
)

# Funzione che calcola HR e intervalli di confidenza rispetto
# a un valore di riferimento, utilizzando la matrice del modello.
#
# La matrice viene costruita con la stessa formula del modello,
# mantenendo esclusivamente le colonne associate ai coefficienti
# effettivamente stimati.

calcola_hr_relativo <- function(modello, nuovi_dati, riferimento = 3) {
  
  # Estraiamo la formula e rimuoviamo la risposta Surv().
  formula_termini <- delete.response(terms(modello))
  
  # Costruiamo la matrice del modello.
  X <- model.matrix(
    formula_termini,
    data = nuovi_dati
  )
  
  # Rimuoviamo l'intercetta, se presente.
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  
  # Manteniamo le colonne nello stesso ordine dei coefficienti.
  nomi_beta <- names(coef(modello))
  X <- X[, nomi_beta, drop = FALSE]
  
  beta <- coef(modello)
  V <- vcov(modello)
  
  # Individuiamo il valore più vicino al riferimento richiesto.
  indice_ref <- which.min(
    abs(nuovi_dati$pnodes - riferimento)
  )
  
  # Differenza tra il vettore dei regressori corrente
  # e quello corrispondente al valore di riferimento.
  X_diff <- sweep(
    X,
    2,
    X[indice_ref, ],
    FUN = "-"
  )
  
  # Log-HR rispetto al riferimento.
  log_hr <- as.vector(X_diff %*% beta)
  
  # Varianza della differenza tra predittori lineari:
  # Var(X_diff * beta) = X_diff * V * t(X_diff).
  se_log_hr <- sqrt(
    pmax(
      0,
      rowSums((X_diff %*% V) * X_diff)
    )
  )
  
  data.frame(
    HR = exp(log_hr),
    low = exp(log_hr - 1.96 * se_log_hr),
    high = exp(log_hr + 1.96 * se_log_hr)
  )
}

# Calcoliamo HR e intervalli di confidenza.
hr_pnodes <- calcola_hr_relativo(
  cox_spline_fixed,
  new_pnodes,
  riferimento = 3
)

new_pnodes <- bind_cols(
  new_pnodes,
  hr_pnodes
)

# Grafico della curva HR.
ggplot(new_pnodes, aes(x = pnodes, y = HR)) +
  geom_ribbon(
    aes(ymin = low, ymax = high),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 3,
    linetype = "dotted"
  ) +
  labs(
    title = "Effetto stimato dei linfonodi positivi",
    subtitle = "Hazard Ratio relativo a 3 linfonodi positivi",
    x = "Numero di linfonodi positivi",
    y = "Hazard Ratio"
  ) +
  theme_minimal()

# Interpretazione:
# HR > 1 indica un hazard stimato maggiore rispetto al riferimento;
# HR < 1 indica un hazard stimato minore.
# Se l'intervallo di confidenza include 1, la differenza rispetto
# al riferimento non è statisticamente significativa al livello
# approssimativo del 5% per quel particolare confronto.
#
# La curva descrive un'associazione aggiustata per le altre covariate.
# Non rappresenta una probabilità di sopravvivenza e non dimostra
# un effetto causale dei linfonodi.


# 21. CONFRONTO MULTIVARIATO: LINEARE, LOGARITMICO E SPLINE

# Confrontiamo tre forme funzionali di pnodes mantenendo uguali
# le altre covariate e la stratificazione.
#
# In questo modo il confronto è più interpretabile rispetto
# al confronto tra modelli con insiemi diversi di covariate.

cox_pnodes_linear <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    pnodes +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d,
  x = TRUE
)

cox_pnodes_log <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    log_pnodes +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d,
  x = TRUE
)

cox_pnodes_spline <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    ns(
      pnodes,
      knots = c(3, 7, 15),
      Boundary.knots = c(1, 51)
    ) +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d,
  x = TRUE
)

# Confronto AIC.
AIC(
  cox_pnodes_linear,
  cox_pnodes_log,
  cox_pnodes_spline
)

# Confronto del rapporto di verosimiglianza tra lineare e spline.
# Questi modelli sono annidati: il test valuta se i termini aggiuntivi
# della spline migliorano l'adattamento rispetto alla forma lineare.

anova(
  cox_pnodes_linear,
  cox_pnodes_spline,
  test = "LRT"
)

# Confronto della capacità discriminante tramite C-index.
c(
  Linear = summary(cox_pnodes_linear)$concordance[1],
  Log = summary(cox_pnodes_log)$concordance[1],
  Spline = summary(cox_pnodes_spline)$concordance[1]
)

# Interpretazione dei risultati precedenti:
# - Il test LRT lineare vs spline ha dato p circa 0,0001.
#   Questo suggerisce che la spline migliora significativamente
#   l'adattamento rispetto alla forma lineare.
# - Il modello logaritmico ha ottenuto AIC circa 2597,77.
# - Il modello spline ha ottenuto AIC circa 2600,18.
# - Il C-index è circa 0,682 per il modello logaritmico e 0,684
#   per la spline.
#
# Questi risultati non sono contraddittori: il confronto LRT riguarda
# il modello lineare contro la spline, mentre l'AIC confronta tutte
# e tre le specificazioni.
#
# Nel confronto complessivo l'AIC favorisce il modello logaritmico,
# mentre il C-index favorisce leggermente la spline.
# La differenza di C-index è molto piccola e non va interpretata
# automaticamente come un vantaggio predittivo sostanziale.
#
# La scelta finale deve considerare insieme adattamento, complessità,
# diagnostica e stabilità delle stime.


# 22. CONFRONTO DELLE CURVE HR: LOGARITMO VS SPLINE

# Costruiamo una griglia comune di valori di pnodes.
# Le altre covariate rimangono costanti.
#
# Il confronto delle curve consente di valutare dove le due
# specificazioni producono stime simili o divergenti.

new_compare <- data.frame(
  pnodes = seq(1, 30, by = 0.1)
)

new_compare$age <- median(d$age)
new_compare$log_tsize <- median(d$log_tsize)
new_compare$log_progrec <- median(d$log_progrec)
new_compare$log_estrec <- median(d$log_estrec)

new_compare$horTh <- factor(
  "no",
  levels = levels(d$horTh)
)

new_compare$menostat <- factor(
  "Pre",
  levels = levels(d$menostat)
)

new_compare$tgrade <- factor(
  "II",
  levels = levels(d$tgrade),
  ordered = is.ordered(d$tgrade)
)

new_compare$log_pnodes <- log1p(new_compare$pnodes)

# Calcoliamo HR e intervalli di confidenza per entrambi i modelli.

hr_log <- calcola_hr_relativo(
  cox_pnodes_log,
  new_compare,
  riferimento = 3
)

hr_spline <- calcola_hr_relativo(
  cox_pnodes_spline,
  new_compare,
  riferimento = 3
)

new_compare <- new_compare |>
  mutate(
    HR_log = hr_log$HR,
    HR_log_low = hr_log$low,
    HR_log_high = hr_log$high,
    HR_spline = hr_spline$HR,
    HR_spline_low = hr_spline$low,
    HR_spline_high = hr_spline$high
  )

# Trasformiamo i dati in formato lungo per ggplot.

plot_data <- bind_rows(
  new_compare |>
    transmute(
      pnodes,
      model = "Logaritmico",
      HR = HR_log,
      low = HR_log_low,
      high = HR_log_high
    ),
  
  new_compare |>
    transmute(
      pnodes,
      model = "Spline",
      HR = HR_spline,
      low = HR_spline_low,
      high = HR_spline_high
    )
)

# Grafico comparativo delle curve.

ggplot(
  plot_data,
  aes(
    x = pnodes,
    y = HR,
    linetype = model
  )
) +
  geom_ribbon(
    aes(
      ymin = low,
      ymax = high,
      group = model,
      fill = model
    ),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 3,
    linetype = "dotted"
  ) +
  labs(
    title = "Confronto tra forma logaritmica e spline",
    subtitle = "Hazard Ratio relativo a 3 linfonodi positivi",
    x = "Numero di linfonodi positivi",
    y = "Hazard Ratio",
    linetype = "Forma funzionale",
    fill = "Forma funzionale"
  ) +
  theme_minimal()

# Interpretazione:
# Curve vicine indicano che i due modelli producono stime puntuali
# simili per quei valori di pnodes.
# Una divergenza nelle code va interpretata con cautela, perché
# le osservazioni sono meno numerose e le stime possono essere instabili.
#
# Gli intervalli di confidenza descrivono l'incertezza delle stime
# all'interno di ciascun modello; non costituiscono un test formale
# della differenza tra le due specificazioni.


# 23. CONFRONTO QUANTITATIVO DELLE CURVE

# Quantifichiamo la differenza tra le stime puntuali nella regione
# compresa tra 1 e 20 linfonodi.
#
# Questo confronto è descrittivo e non sostituisce un test inferenziale.

compare_1_20 <- new_compare |>
  filter(pnodes <= 20) |>
  mutate(
    differenza_HR = HR_spline - HR_log,
    rapporto_HR = HR_spline / HR_log
  )

max_diff <- max(
  abs(compare_1_20$differenza_HR),
  na.rm = TRUE
)

max_diff_rel <- max(
  abs(compare_1_20$rapporto_HR - 1),
  na.rm = TRUE
)

cat(
  "Massima differenza assoluta tra HR:",
  round(max_diff, 3),
  "\n"
)

cat(
  "Massima differenza relativa tra HR:",
  round(max_diff_rel * 100, 1),
  "%\n"
)

# Selezioniamo alcuni valori per facilitare l'interpretazione.

punti_confronto <- c(1, 3, 5, 7, 10, 15, 20)

indici <- sapply(
  punti_confronto,
  function(x) which.min(abs(new_compare$pnodes - x))
)

compare_punti <- new_compare[indici, ] |>
  select(
    pnodes,
    HR_log,
    HR_spline,
    HR_log_low,
    HR_log_high,
    HR_spline_low,
    HR_spline_high
  )

print(compare_punti)

# La tabella consente di verificare in quali punti le stime sono
# maggiormente differenti.
# La conclusione deve considerare sia l'entità della differenza
# sia l'incertezza, evitando di attribuire significato clinico
# a scostamenti piccoli senza ulteriori elementi.


# 24. DIAGNOSTICA DELLE OSSERVAZIONI INFLUENTI

# I residui DFBETA descrivono quanto le stime dei coefficienti
# cambiano in relazione alle singole osservazioni.
#
# Un caso influente non è necessariamente un dato errato:
# può essere un'osservazione rara ma valida.

dfbeta_final <- residuals(
  cox_pnodes_spline,
  type = "dfbeta"
)

# Calcoliamo il massimo valore assoluto dei DFBETA per ogni paziente.
# È un indicatore sintetico dell'influenza, non una regola universale
# per decidere se eliminare un'osservazione.

dfbeta_max <- apply(
  abs(dfbeta_final),
  1,
  max
)

# Individuiamo le dieci osservazioni con il valore più elevato.

influential <- order(
  dfbeta_max,
  decreasing = TRUE
)[1:10]

influential_table <- data.frame(
  row = influential,
  max_dfbeta = dfbeta_max[influential],
  pnodes = d$pnodes[influential],
  time = d$months[influential],
  event = d$event[influential]
)

print(influential_table)

# Grafico dell'indicatore sintetico di influenza.

plot(
  dfbeta_max,
  type = "h",
  main = "Influenza delle osservazioni sul modello di Cox",
  xlab = "Osservazione",
  ylab = "Massimo valore assoluto DFBETA"
)

abline(
  h = 2 / sqrt(nrow(d)),
  lty = 2
)

# La soglia grafica è un riferimento orientativo e non una regola
# automatica di esclusione.
# Nelle analisi precedenti le osservazioni più influenti comprendevano
# pazienti con valori elevati di pnodes, fino a 35-51.
#
# Questi casi si trovano nella coda della distribuzione, dove i dati
# sono meno numerosi e la forma della relazione è più difficile
# da stimare con precisione.
# Prima di ipotizzare un errore, è necessario verificare la qualità
# e la validità delle osservazioni.


# 25. ANALISI DI SENSIBILITÀ

# Escludiamo temporaneamente le cinque osservazioni con il massimo
# indicatore di influenza.
#
# Lo scopo è valutare quanto cambiano le stime del modello.
# Non stiamo affermando che queste osservazioni siano errate
# né che debbano essere eliminate dall'analisi principale.

top5_influential <- influential[1:5]

d_sensitivity <- d[-top5_influential, ]

cox_sensitivity <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    ns(
      pnodes,
      knots = c(3, 7, 15),
      Boundary.knots = c(1, 51)
    ) +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d_sensitivity,
  x = TRUE
)

cat("\n========== MODELLO ORIGINALE ==========\n")
print(summary(cox_pnodes_spline))

cat("\n========== MODELLO DI SENSIBILITÀ ==========\n")
print(summary(cox_sensitivity))

# Confronto del C-index.

cat(
  "\nC-index originale:",
  summary(cox_pnodes_spline)$concordance[1],
  "\n"
)

cat(
  "C-index modello di sensibilità:",
  summary(cox_sensitivity)$concordance[1],
  "\n"
)

# Confrontiamo il numero di pazienti e gli eventi effettivamente
# utilizzati dai modelli.
#
# In coxph, $n indica il numero di osservazioni usate dal modello,
# mentre $nevent indica il numero di eventi.
# nobs() non va utilizzato qui come sinonimo della numerosità
# campionaria, perché per coxph può restituire il numero di eventi.

cat(
  "Pazienti nel modello originale:",
  cox_pnodes_spline$n,
  "\n"
)

cat(
  "Eventi nel modello originale:",
  cox_pnodes_spline$nevent,
  "\n"
)

cat(
  "Pazienti nel modello di sensibilità:",
  cox_sensitivity$n,
  "\n"
)

cat(
  "Eventi nel modello di sensibilità:",
  cox_sensitivity$nevent,
  "\n"
)

# Dai risultati precedenti, dopo l'esclusione temporanea di cinque casi,
# il C-index passa approssimativamente da 0,684 a 0,687.
# La variazione è contenuta.
#
# La stabilità di alcuni coefficienti e del C-index suggerisce che
# i risultati principali non dipendono interamente da quei cinque casi.
# Tuttavia, la curva nelle code può cambiare maggiormente: per questo
# motivo l'interpretazione dei valori elevati di pnodes richiede cautela.
#
# Il modello originale rimane l'analisi principale, salvo motivazioni
# documentate e indipendenti dai risultati per escludere osservazioni.


# 26. TEST CONGIUNTO DELL'EFFETTO DI PNODES

# Il test congiunto valuta l'ipotesi che tutti i coefficienti
# della base spline associata a pnodes siano contemporaneamente zero.
#
# Un risultato significativo indica un'associazione complessiva
# tra pnodes e hazard.
# Non dimostra, da solo, che la relazione sia non lineare.
#
# Per valutare la non linearità è necessario confrontare la spline
# con una specificazione lineare tramite un test appropriato.

# Individuiamo i coefficienti della spline nel modello originale.

termini_spline <- grep(
  "^ns\\(pnodes",
  names(coef(cox_pnodes_spline)),
  value = TRUE
)

ipotesi_spline <- paste0(
  termini_spline,
  " = 0"
)

cat("\n========== TEST CONGIUNTO: MODELLO ORIGINALE ==========\n")

test_pnodes_originale <- car::linearHypothesis(
  cox_pnodes_spline,
  ipotesi_spline,
  test = "Chisq"
)

print(test_pnodes_originale)

# Ripetiamo il test nel modello di sensibilità.

termini_spline_sens <- grep(
  "^ns\\(pnodes",
  names(coef(cox_sensitivity)),
  value = TRUE
)

ipotesi_spline_sens <- paste0(
  termini_spline_sens,
  " = 0"
)

cat("\n========== TEST CONGIUNTO: MODELLO DI SENSIBILITÀ ==========\n")

test_pnodes_sensibilita <- car::linearHypothesis(
  cox_sensitivity,
  ipotesi_spline_sens,
  test = "Chisq"
)

print(test_pnodes_sensibilita)

# Nei risultati precedenti il test congiunto è significativo
# sia nel modello originale sia in quello di sensibilità.
#
# Ciò fornisce evidenza di un'associazione complessiva tra pnodes
# e hazard anche dopo l'esclusione temporanea dei casi influenti.
# Non costituisce, però, una prova autonoma di non linearità:
# il test congiunto verifica che l'effetto sia diverso da zero,
# non che la forma spline sia superiore alla forma lineare.
#
# Per la non linearità resta centrale il confronto tra modello
# lineare e spline tramite il test LRT.


# 27. CONFRONTO DELLE CURVE: MODELLO ORIGINALE VS SENSIBILITÀ

# Confrontiamo le curve HR stimate dai due modelli usando gli stessi
# valori delle covariate e lo stesso riferimento.
#
# Limitiamo la visualizzazione a 1-30 linfonodi per non enfatizzare
# eccessivamente la coda estrema, dove le stime sono meno precise.

pnodes_grid <- seq(1, 30, by = 1)

reference_data <- data.frame(
  pnodes = pnodes_grid,
  age = median(d$age),
  log_tsize = median(d$log_tsize),
  log_progrec = median(d$log_progrec),
  log_estrec = median(d$log_estrec),
  
  horTh = factor(
    rep("no", length(pnodes_grid)),
    levels = levels(d$horTh)
  ),
  
  menostat = factor(
    rep("Pre", length(pnodes_grid)),
    levels = levels(d$menostat)
  ),
  
  tgrade = factor(
    rep("II", length(pnodes_grid)),
    levels = levels(d$tgrade),
    ordered = is.ordered(d$tgrade)
  )
)

# Calcoliamo le curve HR relative a 3 linfonodi positivi.
# Utilizziamo la stessa funzione impiegata nel confronto log-spline.

hr_original <- calcola_hr_relativo(
  cox_pnodes_spline,
  reference_data,
  riferimento = 3
)

hr_sensitivity <- calcola_hr_relativo(
  cox_sensitivity,
  reference_data,
  riferimento = 3
)

reference_data$HR_original <- hr_original$HR
reference_data$HR_sensitivity <- hr_sensitivity$HR

# Rappresentiamo le curve.

plot(
  reference_data$pnodes,
  reference_data$HR_original,
  type = "l",
  lwd = 2,
  ylim = range(
    c(
      reference_data$HR_original,
      reference_data$HR_sensitivity
    ),
    finite = TRUE
  ),
  xlab = "Numero di linfonodi positivi",
  ylab = "Hazard Ratio relativo a 3 linfonodi",
  main = "Analisi di sensibilità della relazione con pnodes"
)

lines(
  reference_data$pnodes,
  reference_data$HR_sensitivity,
  lwd = 2,
  lty = 2
)

abline(h = 1, lty = 3)

legend(
  "topleft",
  legend = c(
    "Modello originale",
    "Modello senza 5 osservazioni influenti"
  ),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

# Interpretazione:
# Se le curve sono simili, le stime puntuali sono relativamente stabili
# rispetto all'esclusione temporanea delle osservazioni selezionate.
# Se divergono nella coda, la forma stimata in quella regione è più
# sensibile alla presenza dei casi influenti.
#
# Questa analisi non identifica quale modello sia corretto e non
# giustifica automaticamente l'esclusione dei pazienti.
# Serve a mostrare ai lettori dove le conclusioni sono più robuste
# e dove richiedono maggiore prudenza.


# 28. CONFRONTO FINALE DEI MODELLI

# Ricostruiamo i tre modelli finali con le stesse covariate,
# lo stesso campione e la stessa stratificazione.
#
# In questo modo possiamo confrontare in maniera coerente le tre
# forme funzionali considerate per pnodes.

cox_linear_pnodes <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    pnodes +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d,
  x = TRUE
)

cox_log_pnodes <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    log_pnodes +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d,
  x = TRUE
)

cox_spline_fixed <- coxph(
  Surv(months, event) ~
    age +
    log_tsize +
    ns(
      pnodes,
      knots = c(3, 7, 15),
      Boundary.knots = c(1, 51)
    ) +
    log_progrec +
    log_estrec +
    horTh +
    strata(menostat) +
    strata(tgrade),
  data = d,
  x = TRUE
)

# Verifichiamo numerosità ed eventi nei modelli finali.
# Le numerosità devono essere uguali perché i modelli utilizzano
# le stesse osservazioni e non ci sono dati mancanti.

cat(
  "Pazienti nel modello lineare:",
  cox_linear_pnodes$n,
  "| Eventi:",
  cox_linear_pnodes$nevent,
  "\n"
)

cat(
  "Pazienti nel modello logaritmico:",
  cox_log_pnodes$n,
  "| Eventi:",
  cox_log_pnodes$nevent,
  "\n"
)

cat(
  "Pazienti nel modello spline:",
  cox_spline_fixed$n,
  "| Eventi:",
  cox_spline_fixed$nevent,
  "\n"
)

# Confronto finale tramite AIC.
AIC(
  cox_linear_pnodes,
  cox_log_pnodes,
  cox_spline_fixed
)

# Confronto LRT tra modello lineare e spline.
anova(
  cox_linear_pnodes,
  cox_spline_fixed,
  test = "LRT"
)

# Confronto del C-index.
c(
  Linear = summary(cox_linear_pnodes)$concordance[1],
  Log = summary(cox_log_pnodes)$concordance[1],
  Spline = summary(cox_spline_fixed)$concordance[1]
)

# Riepilogo completo del modello logaritmico.
summary(cox_log_pnodes)

# Livelli delle variabili categoriche.
levels(d$horTh)
levels(d$menostat)
levels(d$tgrade)

# Distribuzione finale di pnodes e della sua trasformazione.
summary(d$pnodes)
summary(d$log_pnodes)

# CONSIDERAZIONI FINALI
#
# 1. L'analisi Kaplan-Meier evidenzia differenze nella sopravvivenza
#    tra gruppi definiti dal numero di linfonodi positivi.
#
# 2. Il modello di Cox permette di valutare l'associazione tra
#    le covariate e l'hazard, tenendo conto simultaneamente
#    delle altre caratteristiche incluse nel modello.
#
# 3. Il test LRT tra modello lineare e spline suggerisce che una
#    relazione flessibile per pnodes si adatti meglio della forma
#    lineare. Questo risultato riguarda il confronto specifico
#    tra questi due modelli.
#
# 4. Nel confronto complessivo, il modello logaritmico ha un AIC
#    leggermente più basso rispetto alla spline prefissata,
#    mentre la spline presenta un C-index marginalmente superiore.
#    La differenza di capacità discriminante è piccola.
#
# 5. Le osservazioni influenti si concentrano anche nella coda alta
#    di pnodes. L'analisi di sensibilità mostra perché sia importante
#    non sovrainterpretare la forma della curva in regioni con pochi
#    pazienti.
#
# 6. Il test congiunto dei termini spline conferma un'associazione
#    complessiva tra pnodes e hazard, ma non costituisce da solo
#    una prova di non linearità.
#
# 7. La scelta finale del modello deve tenere conto di adattamento,
#    complessità, ipotesi PH, stabilità delle stime e obiettivo
#    dell'analisi. Il modello logaritmico è favorito dall'AIC,
#    mentre la spline permette di rappresentare una relazione
#    più flessibile.
#
# 8. I risultati descrivono associazioni statistiche nel campione
#    analizzato. Non dimostrano causalità e non sostituiscono
#    una valutazione clinica.