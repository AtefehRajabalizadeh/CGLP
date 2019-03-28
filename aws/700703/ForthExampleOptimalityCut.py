#!/usr/bin/env python
# coding: utf-8

# In[96]:


import numpy as np
from random import randint
from math import floor, fabs
from cplex.callbacks import UserCutCallback
import cplex  as CPX
import cplex.callbacks as CPX_CB
import sys

def SecModelParameter ( VarsValue, m, n, my_rhs, my_coef, BranchedVariables):
#     print 'BranchedVariables',BranchedVariables
    unbranched=list()
    c=[]
    d=(np.zeros(((2 * m) + (4 * n), 2*n)))
    halatest=[]
    a=[]
    newrow=[]
    zer=(np.zeros(((2 * m) + (4 * n), 2*n)))
    for i in range (n):
        if i not in (BranchedVariables):
            unbranched.append(i)
    NodeNum=n-len(unbranched)
    
#     print 'unbranched', unbranched
    for i in unbranched:   
        for h in range (m):
            c.append(0)
        for h in range (m):
            lin=0            
            for j in range (n):
#                 print 'i', i, 'h', h,'j',j
#                 print 'j1',j
                if j in (BranchedVariables):
#                     print 'j2',j
                    sv = VarsValue[BranchedVariables.index(j)]
#                     print 'sv', sv
                    lin =lin+ (sv * my_coef[h][j]) 
#                     print 'lin', lin
#                 print 'sv'
            lin=lin- my_rhs[h]
            c.append(lin)                
        for h in range (n):
            c.append(0)
        for h in range (n):
            c.append(0)

        for jj in range (n):
            lin2=-1
            if jj in (BranchedVariables):
                sv = VarsValue[BranchedVariables.index(jj)]
                lin2 =lin2+ sv
            c.append(lin2)
            
        for jj in range (n):
            lin2=0
            if jj in (BranchedVariables):
                lin2 =-1*VarsValue[BranchedVariables.index(jj)]
            c.append(lin2)

              
    k=NodeNum-1 

    for i in unbranched: 
        k=k+1
        for t in range(m):              
            for j in range(n):    
                d[(k-NodeNum)*(2*m+4*n)+t][j+n]=my_coef[t][j] 
                
        for t in range(m):  
            d[(k-NodeNum)*(2*m+4*n)+t][n+i]=d[(k-NodeNum)*(2*m+4*n)+t][n+i]-my_rhs[t]   
                
        for t in range(m):  
            for kk in unbranched:
                d[(k-NodeNum)*(2*m+4*n)+t+m][kk]=my_coef[t][kk] 
                

        for t in range(m):  
            for j in range(n):
                d[(k-NodeNum)*(2*m+4*n)+t+m][j+n]=-1*my_coef[t][j] 
                 
        for t in range(m):  
            d[(k-NodeNum)*(2*m+4*n)+t+m][n+i]=d[(k-NodeNum)*(2*m+4*n)+t+m][n+i]+my_rhs[t] 
                
        for t in range(n):              
            d[(k-NodeNum)*(2*m+4*n)+t+2*m][i]=-1 
            
        for t in range(n):  
            
            d[(k-NodeNum)*(2*m+4*n)+2*m+t][n+t]=1 
        
            
        for t in range(n):  
            
            d[(k-NodeNum)*(2*m+4*n)+2*m+n+t][n+t]=-1 
        
        for mm in unbranched: 
             
            d[(k-NodeNum)*(2*m+4*n)+2*m+2*n+mm][mm]=1  
            
        for t in range(n):  
            
            d[(k-NodeNum)*(2*m+4*n)+t+2*n+2*m][i]=d[(k-NodeNum)*(2*m+4*n)+t+2*n+2*m][i]+1            
            

        for t in range(n):  
            
            d[(k-NodeNum)*(2*m+4*n)+2*m+2*n+t][n+t]=-1 

        for t in unbranched:  
            d[(k-NodeNum)*(2*m+4*n)+2*m+3*n+t][t]=-1  
               
        
        for t in range(n):  
            
            d[(k-NodeNum)*(2*m+4*n)+2*m+3*n+t][n+t]=1 
            
        d=np.vstack([d, zer])                          
            
    d = np.delete(d, slice(((n - NodeNum)*((2 * m) + (4 * n))),((n-NodeNum+1)*((2 * m) + (4 * n)))), axis=0)
#     print 'c',c
#     print 'unbranched',unbranched
    return c, d,unbranched
#     print c, d,unbranched


# In[97]:


def SecondModel (cc,dd,NotBranched):
#     print 'cc,dd,NotBranched',cc,dd,NotBranched
    numrows = len(dd)    # 3 rows in your example
    numcols = len(dd[0])
    VarNum=numcols/2
    ConstNum=numrows/len(NotBranched)
    my_obj=[]
    my_ub=[]
    my_lb=[]
    my_colnames=[]
    my_rownames=["c1"]
    my_coef=[] 
    my_row=[]
    my_rhs=[]
    my_sense=""
    my_Nmah=list()
    
    for i in range(numrows):
        my_obj=cc
        my_ub.append(1)
        my_lb.append(0)
        my_colnames.append("y"+str(i))

    sumarray=[]     
    NotBranched.sort()    
#     sumshode=list()
    for t in range (VarNum):
        sumlist=list()
        for j in range(numrows):
#             print 't,j',t,j
            if t in NotBranched: 
                indj=NotBranched.index(t)
#                 print 'indj',indj
                indjrange=range(indj*ConstNum,(indj+1)*ConstNum)
#                 print 'indjrange',indjrange
                if j in indjrange:
#                     print 'dd[j,t]+dd[j,t+VarNum]',j,'ll',dd[j,t]+dd[j,t+VarNum]
                    sumlist.append(dd[j,t]+dd[j,t+VarNum])
                    dd[j,t+VarNum]=0
                else:
                    sumlist.append(dd[j,t])
            else:
                sumlist.append(dd[j,t])
        sumarray.append(sumlist)
#     print 'sumarray',sumarray
#     print 'newDD',dd
  
#     print 'notbranched',NotBranched
    lenNotbranched=len(NotBranched)
   
    
    for t in range(lenNotbranched):  
        for i in range (VarNum):
            sumlist=list()
            if t>0:
                for p in range(t):
                    for j in range (ConstNum):
                        sumlist.append(0) 
                        
            for j in range(ConstNum):
                sumlist.append(dd[t*ConstNum+j,i+VarNum])
            if i in NotBranched and i > NotBranched[t]:
                indx=NotBranched.index(i)
                indy=NotBranched.index(NotBranched[t])
                
                for x in range(indx-indy-1):
                    for j in range (ConstNum):
                        sumlist.append(0) 
                for j in range(ConstNum):
                    sumlist.append(dd[(indx)*ConstNum+j,NotBranched[t]+VarNum])
                    dd[(indx)*ConstNum+j,NotBranched[t]+VarNum]=0
#                 sumarray.append([t+indx,indy])
                if t<lenNotbranched-2 :
                    for p in range(lenNotbranched-indx-1):
                        for j in range (ConstNum):
                            sumlist.append(0)
            else:
                if t<lenNotbranched-1 :
#                     print 'lenNotbranched',lenNotbranched
                    for p in range(lenNotbranched-t-1):
                            for j in range (ConstNum):
                                sumlist.append(0) 
#             print 't,i',t,',',i
#             print 'sumlist p2',sumlist
#             print 'lensumlist p2',len(sumlist)
            sumarray.append(sumlist)
#     print 'part2' ,sumarray  
    
    rows=[]
    for j in range(numrows):
        my_row.append(1)
    rows.append([range(numrows),my_row])  
    def populatebyrow(SecProb):
        SecProb.objective.set_sense(SecProb.objective.sense.maximize)
        SecProb.variables.add(obj = my_obj,lb=my_lb, names = my_colnames)
        
        for i in range (len(sumarray)):
            if all(k == 0 for k in sumarray[i])!=True:
#                 print 'i',i
                sumline=[]
                sumline.append([range(len(sumarray[i])),sumarray[i]])
#                 print 'sumline',sumline
                SecProb.linear_constraints.add(lin_expr = sumline, senses =  "E",
                                rhs = [0], names = ["u"+str(i)])
        SecProb.linear_constraints.add(lin_expr = rows, senses = "L",
                                             rhs = [1], names = ["t0"])
    
    my_SecProb = CPX.Cplex()
    handle = populatebyrow(my_SecProb)
    my_SecProb.solve()
#     my_SecProb.write("lpex2.lp")
    SecVariables    = my_SecProb.solution.get_values()
    SecObjective=my_SecProb.solution.get_objective_value()
#     print 'SecVariables,SecObjective,my_colnames',SecVariables, SecObjective,my_colnames
    return SecVariables, SecObjective,my_colnames


# In[98]:


def GeneratintCutParameters ( m, n, my_rhs, my_coef, BranchedVariables,SecVariables): 
    BranchedLen=len(BranchedVariables)
    unbranched=list()
    for i in range (n):
        if i not in (BranchedVariables):
            unbranched.append(i)  
    LeUnbranched=len(unbranched)                   
#     NodeNum=n-len(unbranched)
#     k=NodeNum-1 
    lhs=(np.zeros((LeUnbranched*((2 * m) + (4 * n)),BranchedLen)))
    rhs=(np.zeros((LeUnbranched*((2 * m) + (4 * n)),1)))

    for i in range (LeUnbranched):   
        for h in range (m):
            for j in (BranchedVariables):
                lhs[i*(2*m+4*n)+m+h][BranchedVariables.index(j)]=my_coef[h][j]*SecVariables[i*(2*m+4*n)+m+h]  
                rhs[i*(2*m+4*n)+m+h]=my_rhs[h]*SecVariables[i*(2*m+4*n)+m+h]  
            
        for h in range (n):
            for j in (BranchedVariables):
                lhs[i*(2*m+4*n)+2*m+2*n+j][BranchedVariables.index(j)]=SecVariables[i*(2*m+4*n)+2*m+2*n+j]  
                rhs[i*(2*m+4*n)+2*m+2*n+h]=SecVariables[i*(2*m+4*n)+2*m+2*n+h]  
            
#         for h in range (n):
        for j in (BranchedVariables):
            lhs[i*(2*m+4*n)+2*m+3*n+j][BranchedVariables.index(j)]=-1*SecVariables[i*(2*m+4*n)+2*m+3*n+j]  

    #     print 'lefttt',lhs,'rightttt', rhs
    LenLeft=len(lhs)
    LenRight=len(rhs)
    LenSec=len(SecVariables)
#     print 'LenLeft',LenLeft
#     print'LenRight',LenRight
#     print 'LenSec',LenSec
#     indicesLhs=[]
#     for j in range (8):
#         indicesLhs.append([i for i in range(LenLeft) if lhs[i][j]!=0]) 
#     print indicesLhs
#     indicesRhs = [(i,rhs[i]) for i in range(LenRight) if rhs[i]>0]
#     print 'indicesRhs',indicesRhs
#     indiceSecVariables=[(i , SecVariables[i]) for i in range(LenSec) if SecVariables[i]!=0]
#     ValueSecVariables=[SecVariables[i] for i in range(LenSec) if SecVariables[i]!=0]
#     print 'positiveLeftHand',indicesLhs 
#     print 'indiceSecVariables',indiceSecVariables
#     print 'ValueSecVariables', ValueSecVariables
    FinalLhs=np.sum(lhs, axis=0)
    FinalRhs=np.sum(rhs, axis=0)
#     print 'FinalLhs',FinalLhs, 'FinalRhs',FinalRhs
    return FinalLhs,FinalRhs



# In[99]:


class MyCut(CPX_CB.UserCutCallback):

    def __call__(self):

        branched=list()
        values=list()
        branchedNames=list()
        global BranchedMatrix
        senses= "L"
#         print 'self.get_lower_bounds(s)',self.get_lower_bounds()
#         print 'self.get_upper_bounds(s)',self.get_upper_bounds()
#         print 'objective value', self.get_objective_value()
#         print 'values', self.get_values()
#         print 'get_num_nodes', self.get_num_nodes()
        with open('text_OptimalityCut_700703.txt', 'w') as f:
            f.write('{}:{}\n'.format('Number of explored node',self.get_num_nodes()))
        f.closed
        lb=self.get_lower_bounds()
        ub=self.get_upper_bounds()               
        for j in range(len(lb)):
            if lb[j]==ub[j]:
                branched.append(j)
                values.append(lb[j])
                branchedNames.append(VarsName[j])
                
#         print 'values',values
#         print 'branchedNames',branchedNames
#         print 'branched',branched
#         print 'self.get_incumbent_objective_value()',self.get_incumbent_objective_value()
        if len(values)>0:  
            BranchedMatrix.append([branched,values])
#             BranchedValue.append(values)
#             print 'BranchedMatrix',BranchedMatrix

            if [BranchedMatrix[-1]]!=[BranchedMatrix[-2]] :
                self.my_rhs[-1]=(self.get_incumbent_objective_value())
                self.my_rhs[-1]=-1*self.my_rhs[-1]
#                 print 'values, self.m, self.n, self.my_rhs, self.my_coef,branched', self.my_rhs 
                SecModelParameterOutput=SecModelParameter(values, self.m, self.n, self.my_rhs, self.my_coef,branched)
                SecondModelOutput=SecondModel (SecModelParameterOutput[0],SecModelParameterOutput[1],SecModelParameterOutput[2])
#                 print 'SecondModelOutput[1]',SecondModelOutput[1]
#                 print 'SecondModelOutput[0]',SecondModelOutput[0]
                if SecondModelOutput[1]>0.000001:
    #                 print 'bayad cut bezane'
                    output=GeneratintCutParameters (self.m, self.n, self.my_rhs, self.my_coef,branched, SecondModelOutput[0])
#                     print 'output',output[0]
#                     print 'output[1]',output[1]
                    lhs=[CPX.SparsePair(ind =branchedNames, val = output[0])]
#                     print 'lhs',lhs[0]
                    self.add(cut=lhs[0],rhs=output[1][0],sense=senses)
#                     print 'cut added'
                    self.times_called += 1


# In[100]:


def admipex3(c):
    
    c.parameters.preprocessing.presolve = 0

    c.parameters.preprocessing.reduce = 0

    c.parameters.preprocessing.linear.set(0)

    c.parameters.preprocessing.numpass = 0

    c.parameters.mip.strategy.presolvenode = -1

    c.parameters.preprocessing.linear = 0
    
def CutCplex(c):
    c.parameters.mip.cuts.disjunctive.set(-1)
    c.parameters.mip.cuts.bqp.set(-1)
    c.parameters.mip.cuts.cliques.set(-1)
    c.parameters.mip.cuts.covers.set(-1)
    c.parameters.mip.cuts.flowcovers.set(-1)
    c.parameters.mip.cuts.gomory.set(-1)
    c.parameters.mip.cuts.gubcovers.set(-1)
    c.parameters.mip.cuts.implied.set(-1)
    c.parameters.mip.cuts.liftproj.set(-1)
    c.parameters.mip.cuts.localimplied.set(-1)
    c.parameters.mip.cuts.mcfcut.set(-1)
    c.parameters.mip.cuts.mircut.set(-1)
    c.parameters.mip.cuts.zerohalfcut.set(-1)
    c.parameters.mip.cuts.pathcut.set(-1)
    
c = CPX.Cplex('700703.lp')
Orgm=c.linear_constraints.get_num()
n=c.variables.get_num()
m=Orgm+1
my_rhs=[]
my_coef=[]
for i in range(Orgm):
    my_row=[]
    for j in range(n):
        my_row.append(c.linear_constraints.get_coefficients(i,j))
    my_coef.append(my_row)

for i in range (Orgm):
    my_rhs.append(c.linear_constraints.get_rhs(i))
my_rhs.append(0)
# negativeObjective=-1*(c.objective.get_linear())
# print 'negativeObjective',negativeObjective
my_coef.append(c.objective.get_linear())
for i in range (n):
    my_coef[-1][i]=-1* my_coef[-1][i]
# for x in my_coef[-1]:
#     print 'xxxxxxxxxxx',x
# my_coef[-1]=my_coef[-1]*-1               
# print 'my_coef',my_coef
BranchOrdering=[]
for i in range (n):
    BranchOrdering.append((i,i+1,c.order.branch_direction.default))
c.order.set(BranchOrdering)
admipex3(c)
CutCplex(c)
VarsName=c.variables.get_names()
c.set_log_stream(sys.stdout)
c.set_results_stream(sys.stdout)
# c.set_log_stream(None)
# c.set_error_stream(None)
# c.set_warning_stream(None)
# c.set_results_stream(None)
# c.parameters.mip.display.set(0)


cutinst=c.register_callback(MyCut)
cutinst.times_called = 0
cutinst.m=m
cutinst.n=n
cutinst.my_rhs=my_rhs
cutinst.my_coef=my_coef
cutinst.VarsName=VarsName
BranchedMatrix=[[[-1], [-0.0]]]

c.parameters.mip.interval.set(1)
c.parameters.preprocessing.linear=0
c.parameters.mip.strategy.search.set(c.parameters.mip.strategy.search.values.traditional)
# start_time = c.get_time()
c.solve()
# end_time = c.get_time()
# print("Atefeh Total solve time (sec.):", end_time - start_time)

# c.write("lpex1.lp")
solution = c.solution

# the following line prints the corresponding string
print solution.status[solution.get_status()]
# print "Objective value = " , solution.get_objective_value()
print
x = solution.get_values(0, c.variables.get_num()-1)
# for j in range(c.variables.get_num()):
#     if fabs(x[j]) > 1.0e-10:
#         print "Column %d: Value = %17.10g" % (j, x[j])

# print "Branch callback was called ", cutinst.times_called, "times"
with open('text_OptimalityCut_700703.txt', 'a') as f:
    f.write('{}\n'.format('m=700, n=70'))
    f.write('{}:{}\n'.format('objective value',solution.get_objective_value()))
    f.write('{}:{}\n'.format('Number of user cut added',cutinst.times_called))
    f.write('{}:{}\n'.format('Status',solution.status[solution.get_status()]))
    f.write('{}:{}\n'.format('Number of Running CGlp',CGlpCall))
    f.write('{}:{}\n'.format('ElapsedTime',ElapsedTime))

    for j in range(c.variables.get_num()):
        if fabs(x[j]) > 1.0e-10:
            f.write('{}{}:{}\n'.format('x',j, x[j]))
f.closed


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




