import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

df = pd.read_csv('nuclides.csv')


print(" " * 40,"-------------------------------------")
print(" " * 40,"|         Zakaria Daoudi            |")
print(" " * 40,"|   List of Radioactive Nuclides    |")
print(" " * 40,"-------------------------------------\n")

print("key",df)  

while True :
 j = input("Choose a key of an element : ")
 choice = int(j)
 c = df.loc[choice].at["Nuclide Name"]

 print("=" * 30)
 print("You chose : ",c)
 matrix_res = df.to_numpy()

 t_halflife = matrix_res[choice,2]
 print(f"The half-time of {c} is {t_halflife} years.")

 print("--"*8)
 print("With how many nuclides your simulation start ?")
 y = input()
 N0= int(y)

 print("--"*4)
 print("The length of the arrays are 100000 by default. Do you want to change it ? [yes/no]")
 a = str(input())
 proceed='yes' or 'YES' or 'Yes'
 not_proceed='no' or 'NO' or 'No'
 if a == proceed :
    print("What is the new length ?")
    L1 = input()
 elif a == not_proceed :
    L1 = 100000  

 print("--"*4)    
 print("Over how many years the simulation will be done ?")
 dt1 = input()    

 print("--"*4)     
 tau= t_halflife / math.log(2)
 L=int(L1)
 dt=int(dt1)/L
 N = [0]*(L+1)
 t = [0]*(L+1)
 N[0]=N0
 t[0]=0.0
 for i in range(L) :
  N[i+1]=N[i]-(dt*N[i])/tau
  t[i+1]=i*dt
  print(i,N[i],t[i])
 
 plt.xlabel('Time (years)')
 plt.ylabel('Number of Nuclides')
 plt.grid(visible=True, which='major', axis='both')
 plt.plot(t, N, linewidth=1.5)
 plt.show()
