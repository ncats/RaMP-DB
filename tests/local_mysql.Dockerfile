# Use the official MySQL 8.0 image from Docker Hub
FROM mysql:8.0

COPY tests/local_mysql.cnf /etc/mysql/conf.d/my.cnf

# these should match what's in local_mysql.dbprops.txt
ENV MYSQL_ROOT_PASSWORD=root_password
ENV MYSQL_DATABASE=ramp
ENV MYSQL_USER=ramp_user
ENV MYSQL_PASSWORD=ramp_password

RUN curl -L -k -o /docker-entrypoint-initdb.d/ramp.sql.gz https://figshare.com/ndownloader/files/44320136

EXPOSE 3306
