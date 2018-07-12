pipeline {
        agent none
        environment {
                QUAY_USER = credentials('quay-robot')
                QUAY_PASS = credentials('quay-robot-token')

        }

        stages {
                stage('Tag Non-master'){
                        agent {label 'master-builder'}
                        when {not {branch 'master'}}
                        steps{
                        //the tag built here is used for references in other stages
                                script{
                                        TAG = //? sh([script: "echo quay.io/encode-dcc/atac-seq-pipeline:${env.BRANCH_NAME}_${env.BUILD_NUMBER}", returnStdout: true]).trim()
                                }
                                echo "On non-master branch"
                        }        
                }
                stage('Tag Master') {
                        agent {label 'master-builder'}
                        when {branch 'master'}
                        steps{
                                script{
                                        TAG = //?  sh([script: "echo quay.io/encode-dcc/atac-seq-pipeline:latest", returnStdout: true]).trim()
                                }
                                echo "On master"
                        }
                }
                stage('Build-nonmaster') {
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'}
                        when {not {branch 'master'}}
                        steps {
                                echo "The tag is $TAG"
                                echo "Going to build a docker image now..."
                                // pull the cache template image (the image is going to stay pretty much the same so it is no need to be dynamic)
                                //sh "docker pull quay.io/encode-dcc/atac-seq-pipeline:v1"
                                sh "docker login -u=${QUAY_USER} -p=${QUAY_PASS} quay.io"
                                //sh "docker build --cache-from quay.io/encode-dcc/atac-seq-pipeline:v1 -f docker_image/Dockerfile -t atac-seq-pipeline ."
                                //sh "docker tag atac-seq-pipeline $TAG"
                                sh "docker push $TAG"
                                sh "docker logout"
                        }
                }
                stage('Build-master') {
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'}
                        when {branch 'master'}
                        steps {
                                echo "Going to build a docker image now..."
                                // pull the cache template image (the image is going to stay pretty much the same so it is no need to be dynamic)
                                //sh "docker pull quay.io/encode-dcc/atac-seq-pipeline:v1"
                                sh "docker login -u=${QUAY_USER} -p=${QUAY_PASS} quay.io"
                               // sh "docker build --cache-from quay.io/encode-dcc/atac-seq-pipeline:v1 -f docker_image/Dockerfile -t atac-seq-pipeline ."
                              //sh "docker tag atac-seq-pipeline quay.io/encode-dcc/atac-seq-pipeline:latest"
                                //sh "docker push quay.io/encode-dcc/atac-seq-pipeline:latest"
                                sh "docker logout"
                        }
                }
                
                //WHAT DOES THIS STAGE DO? ASK IDAN
		stage('Run hiC cromwell') {
                        agent {label 'master-builder'}
			steps { 
				sh "ls -l"
			}
		}
                stage('Run-Task-Level-Tests-Non-Master'){
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'} 
                        steps{
                                 sh """
                                         test/test_task/test.sh test/test_task/test_align.wdl test/test_task/test_align.json $TAG
                                         python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_align.result.json
                                        
                                        
                                         test/test_task/test.sh test/test_task/test_align.wdl test/test_task/test_merge.json $TAG
                                         python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_merge.result.json
                                        
                                         test/test_task/test.sh test/test_task/test_align.wdl test/test_task/test_merge_sort.json $TAG
                                         python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_merge_sort.result.json
                                        
                                         test/test_task/test.sh test/test_task/test_align.wdl test/test_task/test_align.json $TAG
                                         python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_dedup.result.json
                                        
                                        
                                 """
                        }
                }
                //DO I PUT IN WORKFLOW LEVEL TASKS TOO?
        }         
                
	post {
                success {
                        echo 'Post build actions that run on success'
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished successfully."
                }
                failure {
                        echo 'Post build actions that run on failure'
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} failed."
                }
                always {
                        echo 'Post build actions that run always'
                }
	}
}
