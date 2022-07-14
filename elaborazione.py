#il file si deve trovare nella cartella T1 
# script in python che consente di creare un singolo file .xlsx contenente i dati presenti nella cartella T1

import pandas as pd
import os
import glob
files = glob.glob("M*.xlsx")
files.sort()
new = pd.read_excel(files[0])
files = files[1:]

for file in files:
    temp = pd.read_excel(file)
    new = pd.concat([new,temp])

new.to_excel("unito.xlsx",index=False)
