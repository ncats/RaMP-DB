FROM rocker/tidyverse:4.2.1
# tidyverse:4.2.2+ & 4.3+ is giving an error on apt-get update
# E: Problem executing scripts APT::Update::Post-Invoke 'rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true'

RUN lsb_release -a
RUN sudo apt-get update
RUN apt-get install -y libglpk40
RUN apt-get install -y libxml2
RUN sudo apt-get install -y libmysqlclient21

WORKDIR /ramp-db
RUN mkdir /ramp-db/lib

COPY . /ramp-db

#CMD ["bash"]
CMD ["Rscript", "./tests/travis_test_script.R"]
