import numpy as np

class LinearAlgebra:
  
  @staticmethod
  def matrix_multiplication(A,B):
      # A@B = M
      B = np.array(B).T
      A = np.array(A)
      M = np.zeros((len(A), len(B)))
      for i, line in enumerate(A):
          for j, column in enumerate(B):
              M[i,j] = np.sum(line*column)
      return M

  @staticmethod
  def gaussian_elimination(A, B):
    #for square matrix with pivots != 0
      M = np.hstack((A,B))  # Extended matrix [A/B]
      C = M.copy()
      d = [i for i in range(len(A))]
      for i in d:
          pivot = C[i,i]
          if pivot == 0:
              continue
          else:
              for j in d[i+1:]:
                  factor = C[j, i]/pivot
                  C[j] = C[j] - C[i]*factor
  
      #sistem resolution
      c = np.ones(len(A))  # variables vector [x1, x2, x3, ..., xn]
      c[-1] = C[-1,-1]/C[-1,-2]   # gets the value from the variable xn, then using interations discover the other values
      for i in range(1,len(A)):
          if i == len(c):
              break
          else:
              pivot = C[-1-i, -2-i]
              c[-1-i] = (-np.sum(c[-i:]*C[-1-i, -1-i:-1])+C[-1-i, -1])/pivot
      
      return [C, c]
