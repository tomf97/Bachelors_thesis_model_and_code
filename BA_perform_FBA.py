import pulp
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

#see why biomass is produced even though EMP Pathway is disturbed. GY7 is at position 6
#b[6] = (6,(0,0))
#see why cho2 can be knocked out and bimass still be formed
b[48] = (48,(0,0))

#glucose uptake rate
b[56] = (56,(0,4))

#glycerol input:
b[66] = (66,(0,8.78))

#during growth normally no citrate is extreted
b[70] = (70, (0,0))


#The LP-Solver needs the variables to be defined in a exact way
#The function define_variables does that
#Vector a
a = fo.define_Variables(b,A)

#define LP problem, here max or min of objective function is defined
lp_prob = pulp.LpProblem("Maximize",pulp.LpMaximize)

#define the objective function
lp_prob+= a[i_ob_fct]

#get coefficients from the matrix and set Ax = 0
lp_prob = fo.pair_coefficients_with_variables_and_constrains(lp_prob,A,a)
lp_prob.solve(pulp.GLPK_CMD(options=['--nopresol --noscale --std'])) #options=['--nopresol --noscale --std'] for inf solutions

# print(pulp.LpStatus[lp_prob.status])
numberVariables = len(a)
resultVector = list(range(numberVariables))
nonzero_fluxes = list()
zero_fluxes = list()
# save the result, pair flux_Names and place_in_pathways with the values of the flux
for var in lp_prob.variables():
    position = str(var.name).split('a')  # splits a10 into a and 10, where 10 is at position 1
    resultVector[int(position[1])] = (place_in_pathways[int(position[1])],flux_Names[int(position[1])], var.varValue)
    #if flux nonzero, append entry of resultVector
    if(var.varValue != 0):
        nonzero_fluxes.append(resultVector[int(position[1])])
    else: zero_fluxes.append(resultVector[int(position[1])])
    if int(position[1])==76:
        print(var.varValue)


df_values = []
df_nonzeroValues = []
df_names = []

for i in resultVector:
    df_values.append(i[2])
    df_names.append(i[1])

#for plotting the vector we want nonzerofluxes to have a specific colour, so we mark them in a new list
for i in range(len(df_values)):
    if df_values[i] == 0:
        df_nonzeroValues.append(1)
    else:
        df_nonzeroValues.append(0)

#Plot the result vector:
df = pd.DataFrame({'fluxVal': df_values[0:27],
                   'fluxNonZerVal':df_nonzeroValues[0:27],
                   'flName': df_names[0:27]})
print(df)
df2 = pd.DataFrame({'fluxVal': df_values[27:55],
                   'fluxNonZerVal':df_nonzeroValues[27:55],
                   'flName': df_names[27:55]})
print(df2)
df3 = pd.DataFrame({'fluxVal': df_values[55:len(resultVector)],
                   'fluxNonZerVal':df_nonzeroValues[55:len(resultVector)],
                   'flName': df_names[55:len(resultVector)]})
print(df3)

#plot the fluxes of the pathway image in one picture, external fluxes in the other picture
fig, axes = plt.subplots(2,1)
axes[0].scatter(df['flName'], df['fluxVal'], c=df['fluxNonZerVal'], cmap = "bwr") #RdYlBu
axes[1].scatter(df2['flName'], df2['fluxVal'], c=df2['fluxNonZerVal'], cmap = "bwr" ) #RdYlBu

#set name of plot
#axes.set_titel('Fluxes through Genes')
axes[0].set_ylabel('Flux values', fontsize = 9)
axes[0].set_xticklabels(df['flName'], fontsize=9)
axes[0].grid()
axes[1].set_xlabel('Genes/Fluxes', fontsize = 9)
axes[1].set_ylabel('Flux values', fontsize = 9)
axes[1].grid()

for tick in axes[0].get_xticklabels():
    tick.set_rotation(86)

plt.yticks(fontsize = 10)
plt.xticks(fontsize = 9)
plt.xticks(rotation=86)

#only for the final result of the plots
ax = plt.gca()
ax.text(0.01,0.76,"    non zero value\n    zero value",transform=ax.transAxes, bbox=dict(facecolor='white',edgecolor='black',boxstyle='square'))

fig.subplots_adjust(right = 1.5)

fig.tight_layout()
# Saving the figure
plt.savefig("figures/outputFBA.jpg", bbox_inches='tight')

fig2, axes2 = plt.subplots(1,1)
fontSize = 10
axes2.scatter(df3['flName'], df3['fluxVal'], c=df3['fluxNonZerVal'], cmap = "bwr") #RdYlBu
axes2.set_xlabel('Genes/Fluxes', fontsize = fontSize)
axes2.set_ylabel('Flux values', fontsize = fontSize)
axes2.grid()
plt.yticks(fontsize = fontSize)
plt.xticks(fontsize = fontSize)
plt.xticks(rotation=86)
fig2.set_size_inches(8,4)
fig2.tight_layout()

#only for the final result of the plots
ax = plt.gca()
ax.text(0.81,0.86,"    non zero value\n    zero value",transform=ax.transAxes, bbox=dict(facecolor='white',edgecolor='black',boxstyle='square'))

# Saving the figure
plt.savefig("figures/outputFBA_v_ext.jpg", bbox_inches='tight')



#Print Fluxes and specific values in console
print('citExt, b_TAG, CHO2:')
print(resultVector[70])
print(resultVector[i_ob_fct-11])
print(resultVector[48])
print('this should be LRO1')
print(resultVector[38])
print(resultVector[37])
#plt.show()

