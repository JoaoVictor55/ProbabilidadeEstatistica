#============================================================================
# Nome    : funcão geradora de dados para fina de avaliação
# Autor   : José Cláudio Faria/DCET/UESC
# Data    : 10/12/2021 21:38:44
# Objetivo: Gerar dados (data.frame, com estrutura de covariância) para as
#           avaliações práticas de análise exploratória de dados dos cursos
#           introdutórios de estatística
# Require : Matrix, mvtnorm
# email   : <<<joseclaudio.faria@gmail.com>>>
#============================================================================
#.. Observações:
# 1- Muito cuidado ao informar as matrículas: a geração dos dados para
# análise (e subsequente correção) dependem dessa informação correta.
# As matrículas informadas devem ser obrigatoriamente as do grupo.
#
# 2- Grupos com menos que 3 alunos: repetir a(s) última(s) matrícula(s)
# para a(s) matrícula(s) restante(s).
#
# 3- Se encontrar algum problema com as matrículas informadas, tente alterar
# alguma das matrículas. Nesse caso não esquecer de comunicar ao professor
# (por escrito) no corpo da prova.
#
# 4- A função está intencionalmente mal documentada.
#============================================================================

#. Função: gerar_dados
gerar_dados <- function(m1=NULL, m2=NULL,nm3=NULL, n=2e3)
{
  stopifnot(is.numeric(m1) &
              is.numeric(m2) &
              is.numeric(m3))
  
  set.seed(m1 +
             m2 +
             m3)
  
  m_1 <- runif(1,
               min=20,
               max=40)
  
  m_2 <- runif(1,
               min=20,
               max=40)
  
  ## Categórica
  n_cat_1 <- n/10 * sample(4:8,
                           1)
  
  ## Matriz de variâncias e covariâncias
  sigma_1 <- matrix(c(m_1,
                      m_1 / 1.1,
                      m_1 / 1.1,
                      m_2),
                    ncol=2)
  
  sigma_2 <- matrix(c(m_1,
                      -1 * (m_2 / 1.2),
                      -1 * (m_2 / 1.2),
                      m_2),
                    ncol=2)
  
  require(Matrix) # S4
  near_1 <- nearPD(sigma_1)
  
  near_2 <- nearPD(sigma_2)
  
  sigma_1 <- matrix(near_1[['mat']]@x,
                    nc=ncol(sigma_1))
  
  sigma_2 <- matrix(near_2[['mat']]@x,
                    nc=ncol(sigma_2))
  
  ## Escala proporcional
  require(mvtnorm)
  v_pro_1 <- round(rmvnorm(n=n_cat_1,
                           mean=c(m_1,
                                  m_2),
                           sigma=sigma_1),
                   2)
  
  v_pro_2 <- round(rmvnorm(n=(n - n_cat_1),
                           mean=c(m_1,
                                  m_2),
                           sigma=sigma_2),
                   2)
  
  ## Escala categórica
  cat_1 <- rep('M',
               n_cat_1)
  
  cat_2 <- rep('F',
               n - n_cat_1)
  
  v_pro <- c('v_pro_1',
             'v_pro_2')
  
  v_cat <- c('cat_1',
             'cat_2')
  
  ord <- sample(1:2,
                2)
  
  sexo <- c(eval(parse(text=v_cat[ord[1]])),
            eval(parse(text=v_cat[ord[2]])))
  
  
  ## Frame de dados
  res <- as.data.frame(rbind(eval(parse(text=v_pro[ord[1]])),
                             eval(parse(text=v_pro[ord[2]]))))
  
  res <- cbind(res,
               sexo)
  
  colnames(res) <- c('Y1',
                     'Y2',
                     'Sexo')
  
  ## Outlier v_pro_1
  n_out_v1 <- sample(10:20,
                     1)
  
  out_v1 <- sample(1:length(res[, 1]),
                   n_out_v1)
  
  res[, 1][out_v1] <- sample(730:999,
                             n_out_v1)
  
  ## Outlier v_pro_2
  n_out_v2 <- sample(10:30,
                     1)
  
  out_v2 <- sample(1:length(res[, 2]),
                   n_out_v2)
  
  res[, 2][out_v2] <- sample(200:300,
                             n_out_v2)
  
  ## NAs
  res[, 1][sample(1:n, 
                  sample(10:20, 
                         1))] <- NA
  
  res[, 2][sample(1:n, 
                  sample(10:20, 
                         1))] <- NA
  
  res[, 3][sample(1:n, 
                  sample(10:20, 
                         1))] <- NA
  
  ## Negativos
  res[, 1][sample(1:n, 
                  sample(10:20, 
                         1))] <- -999
  
  res[, 2][sample(1:n, 
                  sample(10:20, 
                         1))] <- -999
  
  invisible(res)
}

remove_outlier <- function(x, type=2)
{
  ## Remover NAs
  x <- na.omit(x)
  
  while(1) 
  {
    ## Quartis
    q.t <- quantile(x$Y1,
                    type=type)[2:4]
    
    q.r <- quantile(x$Y2,
                    type=type)[2:4]
    
    ## Distâncias interquartílicas - iqr
    iqr.t <- q.t[3] - q.t[1]
    iqr.r <- q.r[3] - q.r[1]
    
    ## Identificação de outliers
    out <- subset(x,
                  Y1 >= q.t[3] + 1.5 * iqr.t |
                    Y1 <= q.t[1] - 1.5 * iqr.t |
                    Y1 < 0
                  |
                    Y2 >= q.r[3] + 1.5 * iqr.r |
                    Y2 <= q.r[1] - 1.5 * iqr.r |
                    Y2 < 0)
    
    ## Condição de saída
    if (dim(out)[1] == 0)
      return(x)
    
    ## Remoção de outliers
    x <- subset(x,
                Y1 < q.t[3] + 1.5 * iqr.t &
                  Y1 > q.t[1] - 1.5 * iqr.t &
                  Y1 >= 0
                &
                  Y2 < q.r[3] + 1.5 * iqr.r &
                  Y2 > q.r[1] - 1.5 * iqr.r &
                  Y2 >= 0)
  }
}

m1 = 202120000
m2 = 202120000
m3 = 202120000

m1 = 201920240
m2 = 201820212
m3 = 201920072

#1.1.1

png(filename="q1.1", width=900, height=450)

  dad <- gerar_dados(m1, m2, m3)
  par(mfrow=c(1, 2))
  boxplot(dad$Y1, 
          dad$Y2, 
          main='Antes', 
          names=c('Y1', 'Y2'))
  #Remove outliers
  dad <- remove_outlier(dad)
  boxplot(dad$Y1, 
          dad$Y2, 
          main='Após', 
          names=c('Y1', 'Y2'))
  
dev.off()

#1.1.2

png(filename="q1.2", width=900, height=450)

  dad <- gerar_dados(m1, m2, m3)
  dad <- remove_outlier(dad)
  masculinoDados <- subset(dad, sub = (Sexo == "M"))
  femininoDados <- subset(dad, sub = (Sexo == "F"))
  
  par(mfrow=c(1, 2))
  boxplot(masculinoDados$Y1, 
          masculinoDados$Y2, 
          main='Masculino', 
          names=c('Y1', 'Y2'))
  
  boxplot(femininoDados$Y1, 
          femininoDados$Y2, 
          main='Feminino', 
          names=c('Y1', 'Y2'))

dev.off()

#2.1.1

#Seleciona os valores das variáveis de acordo com  sexo
Y1Masculino <- masculinoDados$Y1
Y1Feminino <- femininoDados$Y1

Y2Masculino <- masculinoDados$Y2
Y2Feminino <- femininoDados$Y2                                          

#Calcula a média dos valores.
mediaMasculinoY1 <- round(mean(Y1Masculino), 2)
mediaMasculinoY2 <- round(mean(Y2Masculino), 2)

mediaFemininoY1 <- round(mean(Y1Feminino), 2)
mediaFemininoY2 <- round(mean(Y2Feminino), 2)

#Calcula a mediana dos valores.
medianMasculinoY1 <- round(median(Y1Masculino), 2)
medianMasculinoY2 <- round(median(Y2Masculino), 2)

medianFemininoY1 <- round(median(Y1Feminino), 2)
medianFemininoY2 <- round(median(Y2Feminino), 2)

#Calcula a moda dos valores.
modeMasculinoY1 <- round(mfv(Y1Masculino), 2)
modeMasculinoY2 <- round(mfv(Y2Masculino), 2)

modeFemininoY1 <- round(mfv(Y1Feminino), 2)
modeFemininoY2 <- round(mfv(Y2Feminino), 2)

#Numero
qtdMasculinoY1 <- length(Y1Masculino)
qtdFemininoY1 <- length(Y1Feminino)

qtdMasculinoY2 <- length(Y2Masculino)
qtdFemininoY2 <- length(Y2Feminino)

#Constroi uma matriz e nomea as colunas com os nomes das mediddas de tendência central,
#e as linhas com os nomes das variáveis.
matrizMasculino <- matrix(c(qtdMasculinoY1, 
                           mediaMasculinoY1,
                           medianMasculinoY1,
                           modeMasculinoY1,
                           qtdMasculinoY2,
                           mediaMasculinoY2, 
                           medianMasculinoY2,
                           modeMasculinoY2),
                         ncol=4,
                         byrow=TRUE)

colnames(matrizMasculino) <- c("n","m","md","mo")
rownames(matrizMasculino) <- c("Y1","Y2")

#matriz para o sexo feminino
matrizFeminino <- matrix(c(qtdFemininoY1, 
                          mediaFemininoY1,
                          medianFemininoY1,
                          modeFemininoY1,
                          qtdFemininoY2,
                          mediaFemininoY2, 
                          medianFemininoY2,
                          modeFemininoY2),
                        ncol=4,
                        byrow=TRUE)

colnames(matrizFeminino) <- c("n", "m","md","mo")
rownames(matrizFeminino) <- c("Y1","Y2")

#Tabelas masculino e feminino:
as.data.frame(matrizMasculino)
as.data.frame(matrizFeminino)

#1.2.1
dad <- gerar_dados(m1, m2, m3)
dad <- remove_outlier(dad)
masculinoDados <- subset(dad, sub = (Sexo == "M"))
femininoDados <- subset(dad, sub = (Sexo == "F"))

(fdt(masculinoDados$Y1))
(fdt(femininoDados$Y1))

#1.2.2
x11()

#..FIM..#

#.. Exemplo de uso: gerar_dados
# dad <- gerar_dados(m1=202210000,
#                    m2=202210000,
#                    m3=202210000)
# str(dad)


#============================================================================
# Nome    : funcão geradora de dados: gerar_dados_rl
# Autor   : José Cláudio Faria/DCET/UESC
# Data    : 10/12/2021 21:38:44
# Objetivo: Gerar dados (data.frame) para as avaliações práticas
#           de análise exploratória de dados referentes a regressão linear
#           introdutórios de estatística
# email   : <<<joseclaudio.faria@gmail.com>>>
#============================================================================

#.. Função: gerar_dados_rl
gerar_dados_rl <- function(m1=NULL,
                           m2=NULL,
                           m3=NULL,
                           n=10)
{
  stopifnot(is.numeric(m1) &
              is.numeric(m2) &
              is.numeric(m3))
  
  set.seed(sum(m1,
               m2,
               m3))
  
  X <- seq(0, 10, length=n)
  Y <- 1 + 2*X + -.08*X^2 + rnorm(n)
  
  res <- data.frame(X,
                    Y)
  
  invisible(res)
}

#.. Exemplo de uso: gerar_dados_rl
# dad_rl <- gerar_dados_rl(m1=202210000,
#                          m2=202210000,
#                          m3=202210000)
# str(dad_rl)
# plot(dad_rl,
#      pch=20)


#============================================================================
# Nome    : funcão geradora de dados: gerar_tdf
# Autor   : José Cláudio Faria/DCET/UESC
# Data    : 10/12/2021 21:38:44
# Objetivo: Gerar dados para uma tabela de distribuição de frequências para
#           avaliações práticas de análise exploratória de dados dos cursos
#           introdutórios de estatística
# email   : <<<joseclaudio.faria@gmail.com>>>
#============================================================================

#.. Função: gerar_tdf
gerar_tdf <- function(m1=NULL,
                      m2=NULL,
                      m3=NULL)
{
  stopifnot(is.numeric(m1) &
              is.numeric(m2) &
              is.numeric(m3))
  
  set.seed(sum(m1,
               m2,
               m3))
  
  classes <- c("[10, 020)",
               "[20, 030)",
               "[30, 040)",
               "[40, 050)",
               "[50, 060)",
               "[60, 070)",
               "[70, 080)",
               "[80, 090)",
               "[90, 100)")
  
  X <- c(seq(f=10, 
             t=50, 
             b=10), 
         seq(f=40, 
             t=10, 
             b=-10))
  
  Y <- sample(1:3,
              length(X),
              rep=T)
  
  f <- (X - Y)
  
  
  rfp <- round(100*f/sum(f),
               2)
  
  cfp <- round(100*cumsum(f/sum(f)),
               2)
  
  res <- data.frame(classes,
                    f,
                    rfp,
                    cfp)
  
  names(res) <- c('Classes',
                  'f',
                  'rf(%)',
                  'cf(%)')
  
  invisible(res)
}

#.. Exemplo de uso:  gerar_tdf
#tb <- gerar_tdf(m1=202210000,
#                m2=202210000,
#                m3=202210000)

#str(tb)
#tb

