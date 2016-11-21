"""
Oblicza rozwiązanie układu równań AX = B metodą eliminacji Gaussa.

Korzysta z naiwnego algorytmu zaadaptowanego z książki Numerical Mathematics
and Computing Warda Cheney'a i Davida Kincaida
"""
function gauss_naive!(A, B)
  n = length(B)
  X = Array{eltype(A)}(n)
  for k in 1 : n-1
    for i = k + 1 : n
      mult = A[i, k] / A[k, k]
      A[i, k] = 0 
      for j in k + 1 : n
        A[i, j] -= mult * A[k, j]
      end
      B[i] -= mult * B[k]
    end
  end
  X[n] = B[n] / A[n, n]
  for i in n-1 : -1 : 1
    sum = B[i]
    for j in i + 1 : n
      sum -= A[i, j] * X[j]
    end
    X[i] = sum / A[i, i]
  end
  X
end


"""
Oblicza rozwiązanie układu równań AX = B metodą eliminacji Gaussa.

Korzysta z algorytmu skalowanego elementu głównego, przekształca
macierz wejściową A, tak że w górnej trójkątnej jej części znajduje się
jej postać schodkowa, natomiast w dolnej zapisane są użyte mnożniki
które wykorzystywane są w kroku podstawiania

Algorytm został zaadaptowany z książki Numerical Mathematics and Computing
Warda Cheney'a i Davida Kincaida
"""
function gauss!(A, B, log = false)
  n = length(B)
  l = gauss_eliminate!(A, n)
  if log
    @show A
    @show l
  end
  X = gauss_substitute!(A, B, l, n)
  if log
    @show X
  end
  X
end


"""
Implementuje krok eliminacji, zapisując mnożniki w dolnej trójkątnej
części macierzy A i zwracając permutację wierszy macierzy A
"""
function gauss_eliminate!(A, n)
  s = Array{eltype(A)}(n) # Scale vector
  l = Array{Int}(n)       # Index permutation
  for i in 1 : n
    l[i] = i
    smax = 0
    for j in 1 : n
      smax = max(smax, abs(A[i, j]))
    end
    s[i] = smax
  end
  for k in 1 : n-1
    rmax = 0
    maxidx = k
    for i in k : n
      r = abs(A[l[i], k] / s[l[i]])
      if r > rmax
        rmax = r
        maxidx = i
      end
    end
    t = l[maxidx]
    l[maxidx] = l[k]
    l[k] = t
    for i = k + 1 : n
      mult = A[l[i], k] / A[l[k], k]
      A[l[i], k] = mult 
      for j in k + 1 : n
        A[l[i], j] -= mult * A[l[k], j]
      end
    end
  end
  l
end

"""
Implementuje krok podstawiania
"""
function gauss_substitute!(A, B, l, n)
  X = Array{eltype(A)}(n) # Solution vector
  for k in 1 : n - 1
    for i in k + 1 : n
      B[l[i]] -= A[l[i], k] * B[l[k]]
    end
  end
  X[n] = B[l[n]] / A[l[n], n]
  for i in n-1 : -1 : 1
    sum = B[l[i]]
    for j in i + 1 : n
      sum -= A[l[i], j] * X[j]
    end
    X[i] = sum / A[l[i], i]
  end
  X
end
