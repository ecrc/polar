pipeline {
/*
 * Defining where to run
 */
//// Any:
// agent any
//// By agent label:
//      agent { label 'sandybridge' }

    agent { label 'jenkinsfile' }
    triggers {
        pollSCM('H/10 * * * *')
    }

    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }

    stages {
        stage ('build') {
            steps {
                sh '''#!/bin/bash -le
                    # loads modules
                    module purge
                    module load cmake/3.9.6
                    module load intel/2017
                    module load intelmpi/2017/intel-2017

                    set -x
                    module list

                    export CC=icc # just in case
                    export FC=ifort # just in case
                    export F90=ifort # just in case
                    export I_MPI_CC="$CC"
                    export I_MPI_FC="$FC"
                    export I_MPI_F90="$F90"

                    mkdir -p build
                    cd build && rm -rf ./*
                    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DPOLAR_TESTING:BOOL=ON -DEXTRA_LIBS="ifcore"

                    # build
                    make

                    # install
                    make install

                '''
            }
        }
        stage ('test') {
            steps {
                sh '''#!/bin/bash -le
                    # loads modules
                    module purge
                    module load cmake/3.9.6
                    module load intel/2017
                    module load intelmpi/2017/intel-2017

                    set -x

                    module list

                    # Delete previous CTest results and run tests
                    rm -rf $WORKSPACE/build/Testing
                    cd $WORKSPACE/build
                    export PATH=$PATH:. 
                    ctest --no-compress-output -T Test
                '''
            }
        }
        stage ('package') {
            steps {
                sh 'cd build && make package'
                archiveArtifacts allowEmptyArchive: true, artifacts: 'build/POLAR-2.0.0-Linux.tar.gz'
            }
        }
    }
    // Post build actions
    post {
        //always {
        //}
        //success {
        //}
        //unstable {
        //}
        //failure {
        //}
        unstable {
                emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build is UNSTABLE", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
        }
        failure {
                emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build FAILED", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
        }
    }
}
