FROM rocker/shiny:3.6.3

RUN apt-get update && apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev default-libmysqlclient-dev

RUN R -e "install.packages(c('devtools'))"

COPY install.R install.R

RUN Rscript install.R

COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
