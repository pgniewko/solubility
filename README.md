## Solubility Challange

### Background
Intrinsic solubility (water solubility): solubility of non-charged molecules, i.e. free acid and base free form. It is required that the solubility is determined in the presence of solid substance.     

Put here an example (table) where I show how solubility changes with pH.     

The project is motivated by the [challange](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00345). 


### This work
To be written.

### Results

### References

[Datasets](#datasets)    
[Methods](#methods)   

### Datasets

1. **Can You Predict Solubilities of Thirty-Two Molecules Using a Database of One Hundred Reliable Measurements?**     
Antonio Llinàs, Robert C. Glen and Jonathan M. Goodman    
*J. Chem. Inf. Modeling 2008, 48, 1289-1303*      
[[paper]](https://pubs.acs.org/doi/10.1021/ci800058v)     
[[website]](http://www-jmg.ch.cam.ac.uk/data/solubility/)     
**Note: In the test set, SMILES strings for probenecid and pseudoephedrine swapped. Use only soldataswap.xls file.**

2. **ESOL: Estimating Aqueous Solubility Directly from Molecular Structure**    
John S. Delaney           
*J. Chem. Inf. Comput. Sci. 2004, 44, 1000-1005*      
[[paper]](https://pubs.acs.org/doi/10.1021/ci034243x)    
**Note: There are two files D.2008.JCIC.solubility.v[1-2].txt. These files are the same but come from two 
different sources: (i) [Pat Walters Blog](https://github.com/PatWalters/solubility) (ii) [ChemDB](ftp://ftp.ics.uci.edu/pub/baldig/learning/Delaney/)**      

3. **Can You Predict Solubilities of Thirty-Two Molecules Using a Database of One Hundred Reliable Measurements?**     
Jarmo Huuskonen
*J. Chem. Inf. Comput. Sci. 2000, 40, 773-777*      
[[paper]](https://pubs.acs.org/doi/10.1021/ci9901338)     
[[website]](http://cheminformatics.org/datasets/huuskonen/index.html)     
**Note: Quite a few repeats from Delaney Set. Different measurements though.**     

4. **ADME evaluation in drug discovery. 4. Prediction of aqueous solubility based on atom contribution approach**
Tingjun Hou, Ke Xia, Wei Zhang, Xiaojie Xu     
*Journal of Chemical Information and Computer Sciences, 2004, 44, 266-275*       
[[paper]](https://pubs.acs.org/doi/abs/10.1021/ci034184n)     
[[website]](http://modem.ucsd.edu/adme/databases/databases_logS.htm)      

5. **Development of reliable aqueous solubility models and their application in drug-like analysis**     
Junmei Wang, George Krudy, Tingjun Hou, George Holland, Xiaojie Xu      
*Journal of Chemical Information and Modeling, 2007, 47, 1395-1404*      
[[paper]](https://pubs.acs.org/doi/10.1021/ci700096r)     
[[website]](http://modem.ucsd.edu/adme/databases/databases_logS.htm)     
**Note: In logS database, the aqueous solubility was expressed as logS, where S is the solubility at a temperature of 20-25°C in mol/L. These are two databases for our modeling. In reference [4], the data afforded by Tetko was used. This database includes 1290 organic compounds. The data set was converted from the SMILES flat file representation to the MACCS/sdf structured data file. In reference [5], some new molecules collected from literature were added. This database includes 1708 molecules.**

6. **Can human experts predict solubility better than computers?**      
Samuel Boobier, Anne Osbourn and John B. O. Mitchell    
*Journal of Cheminformatics, 2017, 9:63*       
[[paper]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0250-y)      
[[website]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0250-y)     
**Note: Source code accompany the paper.**



### Methods    

1. MHFP fingerprints; to be tested as an alternative for ECFP4    
**A probabilistic molecular fingerprint for big data settings**    
Daniel Probst, Jean-Louis Reymond    
*Journal of Cheminformatics, 2018, 10:66*     
[[paper]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0321-8)   
[[source code]](https://github.com/reymond-group/mhfp)    

2. SMILES standarizer.   
[[source code]]()    

