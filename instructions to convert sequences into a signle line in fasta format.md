1. Open the terminal: if it is already open, open a new window using cdm+t
2. Install the program seqtk: type in the terminal brew install seqtk, press enter
3. Once it is install set the working directory by typing the files path. EXAMPLE:
   carolinagomez@Carolinas-MacBook-Pro-7 ~ % cd /Users/carolinagomez/Documents/Bombus\ project/Datos/Plant-volatile-analysis/Plant-volatile-data/
4. USE "LS" to check what files are inside of that folder
5. use the command nano to create a new file in the terminal and paste the sequences you want to convert. example: nano test.fasta
6. use "ls" to check the nano file is created
7. once created type the following: `seqtk seq test.fasta > Plant_seq_may5_singleline.fasta`


# instructions to align the sequences

1. access to clustal omega: https://www.ebi.ac.uk/jdispatcher/msa/clustalo
2. in sequence type select DNA
3. Attached the .fasta file
4. Select Pearson/FASTA in the output format
5. press Submit
6. Download the alingment selecting the option: ## The alignment in FASTA format converted by Seqret

## Instructions to build the tree in IQ-tree

1. Open the terminal
2. paste the following command:
   /Users/carolinagomez/Desktop/IQ-Tree-trial/iqtree -s /Users/carolinagomez/Desktop/IQ-TREE-trial3/all_plants_alingment.fasta -m MFP
3. press enter
4. Once it finishes running look for Best-fit model 
5. Once you find the best fit model
6. Now paste the following command:
```bash
/Users/carolinagomez/Desktop/IQ-Tree-trial/iqtree -s /Users/carolinagomez/Desktop/IQ-TREE-trial3/all_plants_alingment.fasta -m GTR+F+R7 -bb 1000
```
7. In case it need to over write files type the following:
```bash
/Users/carolinagomez/Desktop/IQ-TREE-trial3/iqtree -s /Users/carolinagomez/Desktop/IQ-TREE-trial3/all_plants_alingment.fasta -m GTR+F+R7 -bb 1000 -redo
```


Instructions to build the tree

