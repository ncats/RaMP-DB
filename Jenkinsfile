pipeline {
    options {
        timestamps()
        skipDefaultCheckout()
    }
    agent {
        node { label 'build && aws && linux'}
    }
    triggers {
        pollSCM('H/5 * * * *')
    }
    stages {
        stage('Deploy Application') {
            agent {
                node { label 'rampdb.ncats.io'}
            }
            steps {
                cleanWs()
                configFileProvider([
                    configFile(fileId: 'ramp-db-properties', targetLocation: 'db.properties'),
                    configFile(fileId: 'ramp-db-install-script', targetLocation: 'install.R'),
                    configFile(fileId: 'ramp-db-start-script', targetLocation: 'startup.R')
                ]) {	
                    sh("""
                        Rscript install.R
                        cp db.properties /usr/local/lib/R/site-library/RaMP/shinyApp/db.properties
                        Rscript startup.R &
                    """)
                }
            }
        }
    }
}
