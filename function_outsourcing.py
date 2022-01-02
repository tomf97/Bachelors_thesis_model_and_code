import pulp
import numpy as np

#form of b: b is a list with elements b[] = (pos,(lower_bound, upper_bound))
#returns a Lp-Variable with the constrains given by b
def define_Variables(b,S):
    a = list(range(S.shape[1]))
    runVar = 0
    for i in range(S.shape[1]):
        max_length = len(b)
        if runVar >= max_length:
            b_index = 1000
        else:
            b_index = b[runVar][0]
        if i == b_index:
            lower_bound = b[runVar][1][0]
            upper_bound = b[runVar][1][1]
            if lower_bound == 'none' and upper_bound == 'none':
                a[i] = pulp.LpVariable('a' + str(i), cat='Continuous')
            elif lower_bound == 'none':
                a[i] = pulp.LpVariable('a' + str(i), upBound=upper_bound,
                                    cat='Continuous')
            elif upper_bound == 'none':
                a[i] = pulp.LpVariable('a' + str(i), lowBound=lower_bound,
                                    cat='Continuous')
            else:
                a[i] = pulp.LpVariable('a' + str(i), lowBound=lower_bound, upBound=upper_bound,
                                    cat='Continuous')
            runVar = runVar + 1
        else:
            a[i] = pulp.LpVariable('a' + str(i), cat='Continuous')
    return a

#S is the stoichiometrix matrix, a is the vector of variables
def pair_coefficients_with_variables_and_constrains(optProblem,Stoi_matrix,var):
    # pair the coefficients with the corresponding variable
    for r in range(Stoi_matrix.shape[0]):
        index = Stoi_matrix.shape[1]
        koeff = list(range(index))
        for i in range(index):
            koeff[i] = (var[i], Stoi_matrix[r, i])
        koeff_affineExpression = pulp.LpAffineExpression(koeff)
        # print(koeff_affineExpression)
        optProblem += koeff_affineExpression == 0
    return optProblem


#b_constraint is a list with elements b[] = (pos,(lower_bound, upper_bound))
#b_update is a element of the form: b_update = (pos,(lower_bound, upper_bound))
#stoi_matrix is the matrix with stoichiometric coefficients. It is needed, because
#   we do not want more constrains than we have variables
def updateConstrains(b_constraint, b_update,stoi_matrix):
    replaced = False
    for i in range(len(b_constraint)):
        if b_constraint[i][0] == b_update[0]:
            b_constraint[i] = b_update
            replaced = True
    if(replaced==False):
        #maybe the first part is not necessary (len(b_constraint) < stoi_matrix.shape[1])
        if (len(b_constraint) < stoi_matrix.shape[1]) and (b_update[0] <= stoi_matrix.shape[1]-1):
            b_constraint.append(b_update)
            #now: sort the new List
            b_constraint.sort(key=lambda sublst: sublst[0])
        else:
            raise ValueError('There are too much constrains!')
    return b_constraint

#minimization function! Here we have only min a9!
#function takes initial constrains and the corresponding matrix, then updates the constrains via gene knockout
#then implements a LpProblem, solves it and saves the result inside a matrix
def geneKnockout_FBA(b,S,c_interval, iObfunction):
    a = define_Variables(b, S)
    # Matrix with the results
    numberVariables = len(a)
    allSolInMatrix = np.zeros((numberVariables, numberVariables), dtype='float32')
    initial_constrains = b.copy()
    for i in range(numberVariables):
        b = initial_constrains.copy()
        b = updateConstrains(b, (i, c_interval), S)
        a = define_Variables(b, S)
        # define LP problem
        lp_prob = pulp.LpProblem("Maximum", pulp.LpMaximize)
        # define the objective function
        lp_prob += a[iObfunction]
        # get coefficients from the matrix
        lp_prob = pair_coefficients_with_variables_and_constrains(lp_prob, S, a)

        lp_prob.solve(pulp.GLPK_CMD(options=['--nopresol --noscale --std']))

        # print(pulp.LpStatus[lp_prob.status])
        numberVariables = len(a)
        resultVector = list(range(numberVariables))
        # save the result
        for var in lp_prob.variables():
            position = str(var.name).split('a')  # splits a10 into a and 10, where 10 is at position[1]
            resultVector[int(position[1])] = var.varValue
        allSolInMatrix[i, :] = resultVector
    return allSolInMatrix

def save_knockout_results_as_matrix_in_txt(fileName, matrix):
    with open(fileName, 'wb') as f:
        np.savetxt(f, matrix, fmt='%.2f') #fmt='%.2f' puts the result into a readable form
    return