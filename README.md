# Analisi_dati_MRI
Questo repository contiene i dati e i codici utilizzati per svolgere un progetto svolto in collaborazione con il Dipartimento di Matematica dell'Universit\`a di Genova (che ci ha fornito i dati)
Lo scopo principale del progetto consiste nella presentazione di un approccio basato sulla classificazione di dati in pazienti con meningioma ed e' stato svolto insieme a Tommaso Tennna per l'esame di Data Mining svolto alla Sapienza

Contenuto della repository
analisi_dati.m-> file di MATLAB che contiene tutte le funzioni utili per le elaborazioni dei dati 
Dati.mat -> Dati leggibili da MATLAB (leggi l'appendice del tex per la spiegazione)
histological_data_feasible.xlsx -> dati istologici dei pazienti
Nelle cartelle T1_* sono contenuti i dati grezzi (con vari algoritmi di quantizzazione) 
I file *_dati contengono i dati presenti nella cartelle T1_* concatenati. Per concatenare abbiamo usato lo script python elaborazione.py
Le cartelle Plot_* contengono i vari plot sia in formato jpg sia utilizzando tikz
