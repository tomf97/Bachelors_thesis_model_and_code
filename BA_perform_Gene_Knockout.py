import matplotlib.pyplot as plt
import numpy as np
import function_outsourcing as fo
#dataframe:
import pandas as pd
#read from excel-file:
import xlrd

#write matrix into np.array
book = xlrd.open_workbook('Examples\\Actual_Matrix.xlsx ')
sheet = book.sheet_by_name('Matrix_with_vb_Coe_Cy')
data = [[sheet.cell_value(r, c) for c in range(sheet.ncols)] for r in range(sheet.nrows)]
A = np.zeros((len(data), len(data[0])))
for i in range(len(data)):
    A[i] = data[i]



#write gene names  and place in pathway in two lists - to pair result of FBA with name of genes
book = xlrd.open_workbook('Examples\\Actual_Matrix.xlsx ')
sheet = book.sheet_by_name('Matrix_with_vb_Cy')
data = [[sheet.cell_value(r, c) for c in range(sheet.ncols)] for r in range(sheet.nrows)]
firstRow = data[0]
place_in_pathways = firstRow[2:]
secondRow = data[1]
flux_Names = secondRow[2:]
#print(len(myList2))

#test if objective function is the right flux
# biomass: 86, python begins to count at 0!
i_ob_fct = 86
print(A[70,i_ob_fct])
print(flux_Names[i_ob_fct])


#number of constrains = number of columns of the matrix A
num_constrains = A.shape[1]

#value of constrains
b = list(range(num_constrains))
#we only have forward reactions i.e. the lower limit is set to 0
for i in range(num_constrains):
    if i<=53:
        b[i] = (i, (0, 50)) # max is 99999
    else:
        b[i] = (i, (-50, 50))

#Additional Constrains

#glucose uptake rate
b[56] = (56,(0,4))

#glycerol input:
b[66] = (66,(0,8.78))

#during growth normally no citrate is extreted
b[70] = (70, (0,0))


#iObfunction can only be used if it consits of only one flux. For more, please adjust it in function_outsourcing
R = fo.geneKnockout_FBA(b,A, c_interval=(0,0), iObfunction=i_ob_fct)

fo.save_knockout_results_as_matrix_in_txt('BA_gene_knockout_results.txt', R)

#now: get the important vector and plot it
#depending on my objective function, this can differ!
resultVector = list(range(R.shape[0]))
for i in resultVector:
    resultVector[i] = R[i][i_ob_fct]
#maximum betragsmäßig
#v_max = max(v, key=abs)

#for plotting the vector we want nonzerofluxes to have a specific colour, so we mark them
colorMarker = []
for i in range(len(resultVector)):
    if resultVector[i] == 0:
        colorMarker.append(1)
    else:
        colorMarker.append(0)


df = pd.DataFrame({'fluxVal': resultVector[0:27],
                   'fluxNonZerVal':colorMarker[0:27],
                   'flName': flux_Names[0:27]})
print(df)
df2 = pd.DataFrame({'fluxVal': resultVector[28:55],
                   'fluxNonZerVal':colorMarker[28:55],
                   'flName': flux_Names[28:55]})
df3 = pd.DataFrame({'fluxVal': resultVector[55:len(resultVector)],
                   'fluxNonZerVal':colorMarker[55:len(resultVector)],
                   'flName': flux_Names[55:len(resultVector)]})
print(df3)

#plot the fluxes of the pathway image in one picture, external fluxes in the other picture
fig, axes = plt.subplots(2,1)
axes[0].scatter(df['flName'], df['fluxVal'], c=df['fluxNonZerVal'], cmap = "bwr") #RdYlBu
axes[1].scatter(df2['flName'], df2['fluxVal'], c=df2['fluxNonZerVal'], cmap = "bwr" ) #RdYlBu

#set name of plot
#axes.set_titel('Fluxes through Genes')
axes[0].set_ylabel('Biomass formation', fontsize = 9)
axes[0].set_xticklabels(df['flName'], fontsize=9)
axes[0].grid()
axes[1].set_xlabel('Genes/Fluxes', fontsize = 9)
axes[1].set_ylabel('Biomass formation', fontsize = 9)
axes[1].grid()

for tick in axes[0].get_xticklabels():
    tick.set_rotation(86)

plt.yticks(fontsize = 10)
plt.xticks(fontsize = 9)
plt.xticks(rotation=86)
plt.ylim([-0.05,1])


#only for the final result of the plots
#axes[0].set_ylim([-0.05,1])
#axes[1].set_ylim([-0.05,1])
#ax = plt.gca()
#ax.text(0.01,0.76,"    non zero value\n    zero value",transform=ax.transAxes, bbox=dict(facecolor='white',edgecolor='black',boxstyle='square'))

fig.subplots_adjust(right = 1.5)
fig.tight_layout()
# Saving the figure
plt.savefig("figures/outputGKA.jpg", bbox_inches='tight')

fig2, axes2 = plt.subplots(1,1)
axes2.scatter(df3['flName'], df3['fluxVal'], c=df3['fluxNonZerVal'], cmap = "bwr") #RdYlBu
axes2.set_xlabel('Genes/Fluxes', fontsize = 9)
axes2.set_ylabel('Biomass formation', fontsize = 9)
axes2.grid()
plt.yticks(fontsize = 10)
plt.xticks(fontsize = 10)
plt.xticks(rotation=86)
fig2.subplots_adjust(right = 10)
fig2.set_size_inches(8,4)
fig2.tight_layout()


#only for the final result of the plots
#axes2.set_ylim([-0.05, 0.9])
#ax = plt.gca()
#ax.text(0.81,0.86,"    non zero value\n    zero value",transform=ax.transAxes, bbox=dict(facecolor='white',edgecolor='black',boxstyle='square'))


# Saving the figure
plt.savefig("figures/outputGKA_v_ext.jpg", bbox_inches='tight')
#print(df)