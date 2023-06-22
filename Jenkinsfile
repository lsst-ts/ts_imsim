#!/usr/bin/env groovy

pipeline {

    agent {
        // Use the docker to assign the Python version.
        // Use the label to assign the node to run the test.
        // It is recommended by SQUARE team do not add the label to let the
        // system decide.
        docker {
            image 'lsstts/develop-env:develop'
            args "--entrypoint=''"
        }
    }

    options {
      disableConcurrentBuilds()
      skipDefaultCheckout()
    }

    triggers {
        cron(env.BRANCH_NAME == 'develop' ? '0 4 * * *' : '')
    }

    environment {
        // Position of LSST stack directory
        LSST_STACK = "/opt/lsst/software/stack"
        // Pipeline stack version
        STACK_VERSION = "current"
        // XML report path
        XML_REPORT = "jenkinsReport/report.xml"
        // Module name used in the pytest coverage analysis
        MODULE_NAME = "lsst.ts.imsim"
        // PlantUML url
        PLANTUML_URL = "https://github.com/plantuml/plantuml/releases/download/v1.2022.12/plantuml.jar"
        // Authority to publish the document online
        user_ci = credentials('lsst-io')
        LTD_USERNAME = "${user_ci_USR}"
        LTD_PASSWORD = "${user_ci_PSW}"
        DOCUMENT_NAME = "ts-imsim"
    }

    stages {

        stage('Cloning Repos') {
            steps {
                dir(env.WORKSPACE + '/ts_imsim') {
                    checkout scm
                }
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        git clone https://github.com/lsst-ts/ts_wep.git
                    """
                }
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        git clone https://github.com/lsst-ts/ts_ofc.git
                    """
                }
            }
        }

        stage ('Building the Dependencies') {
            steps {
                // When using the docker container, we need to change
                // the HOME path to WORKSPACE to have the authority
                // to install the packages.
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        source ${env.LSST_STACK}/loadLSST.bash
                        cd ../ts_wep/
                        setup -k -r .
                        scons python
                    """
                }
            }
        }

        stage('Unit Tests and Coverage Analysis') {
            steps {
                // Direct the HOME to WORKSPACE for pip to get the
                // installed library.
                // 'PATH' can only be updated in a single shell block.
                // We can not update PATH in 'environment' block.
                // Pytest needs to export the junit report.
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        source ${env.LSST_STACK}/loadLSST.bash
                        cd ../ts_wep/
                        setup -k -r .
                        cd ../ts_ofc/
                        setup -k -r .
                        cd ../ts_imsim/
                        setup -k -r .
                        pytest --cov-report html --cov=${env.MODULE_NAME} --junitxml=${env.WORKSPACE}/${env.XML_REPORT}
                    """
                }
            }
        }
    }

    post {
        always {
            // The path of xml needed by JUnit is relative to
            // the workspace.
            junit "${env.XML_REPORT}"

            // Publish the HTML report
            publishHTML (target: [
                allowMissing: false,
                alwaysLinkToLastBuild: false,
                keepAll: true,
                reportDir: 'ts_imsim/htmlcov',
                reportFiles: 'index.html',
                reportName: "Coverage Report"
            ])

            script{
              withEnv(["HOME=${env.WORKSPACE}"]) {
                def RESULT = sh returnStatus: true, script: """
                  source ${env.LSST_STACK}/loadLSST.bash

                  curl -L ${env.PLANTUML_URL} -o plantuml.jar
                  export PATH_PLANTUML=${env.WORKSPACE}/plantuml.jar

                  pip install sphinxcontrib-plantuml

                  cd ../ts_wep/
                  setup -k -r .

                  cd ../ts_ofc/
                  setup -k -r .

                  cd ../ts_imsim/
                  setup -k -r .

                  package-docs build

                  pip install ltd-conveyor

                  ltd upload --product ${env.DOCUMENT_NAME} --git-ref ${env.BRANCH_NAME} --dir doc/_build/html
                """

              if ( RESULT != 0 ) {
                  unstable("Failed to push documentation.")
              }
            }
          }
        }
        regression {
            script {
                slackSend(color: "danger", message: "${JOB_NAME} has suffered a regression ${BUILD_URL}", channel: "#aos-builds")
            }

        }
        fixed {
            script {
                slackSend(color: "good", message: "${JOB_NAME} has been fixed ${BUILD_URL}", channel: "#aos-builds")
            }
        }
        cleanup {
            // clean up the workspace
            deleteDir()
        }
    }
}
