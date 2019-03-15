#!/bin/bash

echo "Compute Phi"

python IsingEnergy.py 0 1000 10080 Phi_0_1000.csv
python IsingEnergy.py 1000 2000 10080 Phi_1000_2000.csv
python IsingEnergy.py 2000 3000 10080 Phi_2000_3000.csv
python IsingEnergy.py 3000 4000 10080 Phi_3000_4000.csv
python IsingEnergy.py 4000 5000 10080 Phi_4000_5000.csv
python IsingEnergy.py 5000 6000 10080 Phi_5000_6000.csv
python IsingEnergy.py 6000 7000 10080 Phi_6000_7000.csv
python IsingEnergy.py 7000 8000 10080 Phi_7000_8000.csv
python IsingEnergy.py 8000 9000 10080 Phi_8000_9000.csv
python IsingEnergy.py 9000 10080 10080 Phi_9000_10080.csv

echo "All computation finished!"