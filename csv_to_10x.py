# -*- coding: utf-8 -*-
"""
Created on Thu May 12 19:41:50 2022
# Goal is to make 10x format files FSC files from flow cytometry experients.
# In FlowJo, set gates and export files as csv on compensated parameters and channel values.
# For this experiment, we want all CD45+ cells, so I gated on Lymphocytes>Live>CD45.
# Two "samples" of cells. Output is a large table with cells as rows, variables as columns.
@author: hioki
"""
import pandas as pd
import gzip

###First, import the file
inputdir = "inputdirectory/"
outputdir = "outputdirectory/"
fullfile = "export_CD45+.csv"
file = pd.read_csv(inputdir+fullfile)
ncols = file.shape[1]
nrows = file.shape[0]

###Assign extension of cell barcodes
group_val = "-6"

###Assign row (cell) names
for row in file.index :
    file.at[row, "cellID"] = ("cell_" + str(row+1) + group_val)


###Write barcodes.tsv with the row names
f = open((outputdir+"barcodes.tsv"), "w", newline='')
for row in file.index :
    f.write( file.at[row,"cellID"] + "\n" )
f.close()

with open(outputdir+"barcodes.tsv", "rb") as f_in, gzip.open(outputdir+"barcodes.tsv.gz", "wb") as f_out:
    f_out.writelines(f_in)


###Select columns (varialbes)
# "FSC-A","FSC-H","SSC-A","SSC-B-A","SSC-B-H","SSC-H","Comp-APC-A","Comp-APC-Cy7-A","Comp-APC-Fire 810-A","Comp-BUV395-A","Comp-BUV737-A","Comp-BUV805-A","Comp-BV421-A","Comp-BV510-A","Comp-BV650-A","Comp-DAPI-A","Comp-PE-A","Comp-PE-Cy5-A","Comp-PE-Cy7-A","Comp-PE-Dazzle594-A","Comp-Pacific Blue-A","Comp-PerCP-Cy5.5-A","Comp-Super Bright 600-A","Comp-ZsGreen-A","Time"
#Won't use ZsGreen as a variable when making the UMAP plots. But, want to use it later.
#Export ZsGreen values with rownames.
#Get rid of columns that will never be used.

var_dict = { #dictionary of all variables in file. Make sure they're in the same order as columns.
    "FSC-A" : "FSC-A",
    "FSC-H" : "FSC-H",
    "SSC-A" : "SSC-A",
    "SSC-B-A" : "SSC-B-A",
    "SSC-B-H" : "SSC-B-H",
    "SSC-H" : "SSC-H",
    "Comp-APC-A" : "CD172",
    "Comp-APC-Cy7-A" : "CD8",
    "Comp-APC-Fire 810-A" : "B220",
    "Comp-BUV395-A" : "CD4",
    "Comp-BUV737-A" : "CD24",
    "Comp-BUV805-A" : "CD45",
    "Comp-BV421-A" : "NK1.1",
    "Comp-BV510-A" : "gd",
    "Comp-BV650-A" : "Ly6g",
    "Comp-DAPI-A" : "DAPI",
    "Comp-PE-A" : "SiglecH",
    "Comp-PE-Cy5-A" : "F480",
    "Comp-PE-Cy7-A" : "CD11c",
    "Comp-PE-Dazzle594-A" : "CD3",
    "Comp-Pacific Blue-A" : "Ly6c",
    "Comp-PerCP-Cy5.5-A" : "CD11b",
    "Comp-Super Bright 600-A" : "MHCII",
    "Comp-ZsGreen-A" : "ZsGreen",
    "Time" : "Time"
}

#write values of ZsGreen to a new file
f = open((outputdir+"group"+group_val+"_ZsGreen.csv"), "w", newline='')
for row in file.index :
    f.write( file.at[row,"cellID"] + "," + str(file.at[row,"Comp-ZsGreen-A"]) + "\n" )
f.close()

#write values of DAPI to a new file
f = open((outputdir+"group"+group_val+"_DAPI.csv"), "w", newline='')
for row in file.index :
    f.write( file.at[row,"cellID"] + "," + str(file.at[row,"Comp-DAPI-A"]) + "\n" )
f.close()

#since gate was on CD45+ ("Comp-BUV805-A") cells, normalize values to minimum.
minval = min(file["Comp-BUV805-A"])
tempval = file["Comp-BUV805-A"]
for i in range(0,nrows) :
    file.at[i, "Comp-BUV805-A"] -= minval
tempval = file["Comp-BUV805-A"]

#Delete unneeded columns
file.drop("FSC-H", axis=1, inplace = True)
del var_dict["FSC-H"]
file.drop("SSC-H", axis=1, inplace = True)
del var_dict["SSC-H"]
file.drop("SSC-B-A", axis=1, inplace = True)
del var_dict["SSC-B-A"]
file.drop("SSC-B-H", axis=1, inplace = True)
del var_dict["SSC-B-H"]
#file.drop("Comp-BV510-A", axis=1, inplace = True)
#del var_dict["Comp-BV510-A"]
file.drop("Comp-DAPI-A", axis=1, inplace = True)
file.drop("Comp-ZsGreen-A", axis=1, inplace = True)
file.drop("cellID", axis=1, inplace = True)
file.drop("Time", axis=1, inplace = True)
del var_dict["Time"]

###Write features.tsv with the remaining columns names.
#tab delimited: stain \t Gene-name
f = open((outputdir+"features.tsv"), "w", newline='')
for col in file :
    key = str(col)
    var_val = str(var_dict[key])
    f.write( key + "\t" + var_val + "\n" )
f.close()
with open(outputdir+"features.tsv", "rb") as f_in, gzip.open(outputdir+"features.tsv.gz", "wb") as f_out:
    f_out.writelines(f_in)

###Now onto making the matrix file.
#Sparce matrix, meaning entries with "0" won't be included. However most datapoints will have some value, so will become a large file.
#Iterate each row, then column, and write its value. 1-based coordinates.
#Need to find the number of 0's first, so that the header can be written correctly.
ncols = file.shape[1]
nrows = file.shape[0]
total_values = ncols * nrows
""" If there are zeros in the file, need to subtract from total_values and sparce matrix
zero_count = 0
for col in file :
    for val in col :
        if val == 0 :
            zero_count += 1
"""
f = open((outputdir+"matrix.mtx"), "w", newline='')
#write headers
f.write("%%MatrixMarket matrix coordinate integer general\n")
f.write("%\n")
f.write(str(ncols) + " " + str(nrows) + " " + str(total_values) + "\n")

colnum = 1
rownum = 1
for row in range(0,nrows) :
    for col in file :
        f.write(str(colnum) + " " + str(rownum) + " " + str(file.at[row,col]) + "\n")
        colnum += 1
    colnum = 1
    rownum += 1
f.close()
with open(outputdir+"matrix.mtx", "rb") as f_in, gzip.open(outputdir+"matrix.mtx.gz", "wb") as f_out:
    f_out.writelines(f_in)





