echo "testPerfectGraphRecognition.sh"
if [ ! -r "$1" ]; then
    wget --no-check-certificate "https://users.cecs.anu.edu.au/~bdm/data/$1"
fi
