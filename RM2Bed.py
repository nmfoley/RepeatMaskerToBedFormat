#### python_homework_20180410

###########################################################################################################################################################################################
#
#Take in a RepeatMasker output file and perform the following manipulations:
#
# 1. Convert to bed format with the following columns.
#        A - chrom/contig/scaffold
#        B - chromStart
#        C - chromEnd
#		 C.1 - ADD ORIENTATION
#        D - TE name
#        E - TE Class/Family
#        F - size of TE insertion in bp
#        G - K2P distance (without "kimura =")
#
# 2. Parse the output into six subcategory files: SINEs, Rolling Circle transposons (DNA/RC), DNA transposons (DNA/<anything other than RC>), NonLTR elements, LTR elements, unknowns.
#        i.e. SINEs_rm.out, RC_rm.out, etc.
#
# input file = aVan_rm.out
#
###########################################################################################################################################################################################

##########################################################################
#					Jenny and Nicole's Masterpiece
##########################################################################

import pandas as pd
from pandas import DataFrame, read_csv

#read in the file with only the columns you want and rename them
df = pd.read_csv("aVan_rm.out", sep="\t", skiprows=(0,1,2), usecols=[4,5,6,8,9,10,15], names =["Chrom", "Start", "End", "Orientation","TE name","TE_family_old", "K2P"], header = None)

#fix the kimura thing
df["K2P"] = df["K2P"].str.replace('kimura=', '')

#Get the length of the TE insertion adding 1 because Python can't count
df['Length'] = df['End'] - df['Start'] - 1

# Split the TE_family column to families and subfamilies in a separate dataframe
df1 = pd.DataFrame(df.TE_family.str.split('/', 1).tolist(), columns = ["TE_family", "TE_subfamily"])

#merge the two dataframes together
df2 = pd.concat([df, df1], axis=1)

#correct the strand data so that bedtools wil recognise the negative orientation
df2["Orientation"] = df2["Orientation"].str.replace('C', '-')

#quick re-jig to get the order we want and get rid of old TE family column
df3 = df2[["Chrom", "Start", "End", "Orientation","TE name","TE_family", "TE_subfamily", 'Length', "K2P"]]

#Split output to TE families/subfamilies of interest and write to a .csv bed file
DNA=df3.loc[df3["TE_family"] == "DNA"]
DNA.to_csv("DNA.bed", sep="\t", header = False, index = False)
SINE=df3.loc[df3["TE_family"] == "SINE"]
SINE.to_csv("SINE.bed", sep="\t", header = False, index = False)
NonLTR=df3.loc[df3["TE_family"] == "NonLTR"]
NonLTR.to_csv("NonLTR.bed", sep="\t", header = False, index = False)
Unknown=df3.loc[df3["TE_family"] == "Unknown"]
Unknown.to_csv("Unknown.bed", sep="\t", header = False, index = False)
LTR=df3.loc[df3["TE_family"] == "LTR"]
LTR.to_csv("LTR.bed", sep="\t", header = False, index = False)
DNA_RC=DNA.loc[DNA["TE_subfamily"] == "RC"]
DNA_RC.to_csv("DNA_RC.bed", sep="\t", header = False, index = False)
NON_DNA_RC=DNA.loc[DNA["TE_subfamily"] != "RC"]
NON_DNA_RC.to_csv("NON_DNA_RC.bed", sep="\t", header = False, index = False)
