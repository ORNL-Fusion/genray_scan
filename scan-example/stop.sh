ps -ef | grep xgenray | grep -v grep | awk '{print $2}' | xargs kill -9
