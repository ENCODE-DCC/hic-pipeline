pipeline {
        agent none

        stages {
		stage('Run hiC cromwell') {
                        agent {label 'master-builder'}
			steps { 
				sh "ls -l"
			}
		}
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
