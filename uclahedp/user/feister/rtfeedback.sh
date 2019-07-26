#!/bin/bash
# Indefinitely run the python script

while :
do
    python rtfeedback.py
    echo "Synchronized."
    sleep 0.25
done
