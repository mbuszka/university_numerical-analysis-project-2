function build_matrix(x, y)
  n = size(x,1) - 1
  T = eltype(x)
  h = Array{T}(n)
  a, b, f, g = similar(h), similar(h), similar(h), similar(h)

  for i in 1 : n
    h[i] = x[i+1] - x[i]
    f[i] = (y[i+1] - y[i]) / h[i]
  end
  for i in 1 : n - 1
    g[i] = 6(f[i+1] - f[i])
  end
  g[n] = 6(f[1] - f[n])

  for i in 1 : n-1
    a[i] = 2(h[i] + h[i+1])
    b[i] = h[i+1]
  end
  a[n] = 2(h[1] + h[n])
  b[n] = h[1]
  a, b, g, h, f
end

function spline_periodic(x, y; debug = false)
  n = size(x,1) - 1
  T = eltype(x)

  # Tworzymy macierz dla naszego układu równań
  a, b, g, h, f = build_matrix(x, y)

  # Rozwiązujemy układ równań
  k = solve_matrix(a, b, g)
  s = build_spline(x, y, f, h, k, debug)
  s
end

function build_spline(x, y, f, h, k, debug = false)
  T = eltype(x)
  n = size(h, 1)
  w = Array{Array{T}}(n)
  w[1] = Array{T}(4)
  w[1][1] = y[1]
  w[1][2] = -h[1] * k[1] / T(6) - h[1] * k[n] / T(3) + f[1]
  w[1][3] = k[n] / T(2)
  w[1][4] = (k[1] - k[n]) / (6h[1])
  for i in 2 : n
    w[i] = Array{T}(4)
    w[i][1] = y[i]
    w[i][2] = -h[i] * k[i] / T(6) - h[i] * k[i-1] / T(3) + f[i]
    w[i][3] = k[i-1] / T(2)
    w[i][4] = (k[i] - k[i-1]) / (6h[i])
  end
  s = t -> begin
    if debug
      v = pder
    else
      v = pval
    end
    for i in 1 : n
      if t < x[i+1]
        return v(w[i], t - x[i])
      end
    end
    v(w[1], t - x[n+1])
  end
  s
end

function pval(a, x)
  n = size(a, 1)
  y = a[n]
  for i in n - 1 : -1 : 1
    y = y * x + a[i]
  end
  y
end

function pder(a, x)
  n = size(a, 1)
  T = typeof(x)
  y, d1, d2 = a[n], zero(T), zero(T)
  for i in n - 1 : -1 : 1
    d2 = d2 * x + d1
    d1 = d1 * x + y
    y  = y  * x + a[i]
  end
  y, d1, 2d2
end

# function build_matrix(x, y)
#   n = size(x,1) - 1
#   T = eltype(x)
#   h = Array{T}(n)
#   a, b, c, f, g = similar(h), similar(h), similar(h), similar(h), similar(h)
#
#   for i in 1 : n
#     h[i] = x[i+1] - x[i]
#     f[i] = (y[i+1] - y[i]) / h[i]
#   end
#
#   for i in 1 : n - 1
#     g[i] = 3(h[i] * f[i+1] + h[i+1] * f[i])
#   end
#   g[n] = 3(h[n] * f[1] + h[1] * f[n])
#
#   for i in 1 : n - 2
#     a[i] = 2(h[i] + h[i+1])
#     b[i] = h[i]
#     c[i] = h[i+2]
#   end
#   a[n-1], a[n] = 2(h[n-1] + h[n]), h[1] + h[n]
#   b[n-1], b[n] = h[n-1], 2h[n]
#   c[n-1], c[n] = 2h[1], h[2]
#   a, b, c, g, h, f
# end


# function reduce_matrix(a, b, c)
#   n = size(a, 1)
#   d = similar(a)
#
#   d[1] = 1
#   for i in 2 : n
#     d[i] = d[i-1] * sqrt(b[i-1] / c[i-1])
#   end
#
#   for i in 1 : n - 1
#     b[i] = d[i+1] * c[i] / d[i]
#   end
#   b[n] = d[1] * c[n] / d[n]
#
#   a, b
# end

function solve_matrix(a, b, g)
  T, n = eltype(a), size(a, 1)
  α, β, γ, δ  = similar(a), similar(a), Array{T}(n - 2), Array{T}(n - 1)

  α[1] = a[1]
  β[1] = b[n]
  δ[1] = β[1] / α[1]
  α[n] = a[n] - δ[1] * β[1]
  for i in 2 : n - 1
    γ[i-1] = b[i-1] / α[i-1]
    α[i] = a[i] - γ[i-1] * b[i-1]
    if i == n - 1
      β[i] = b[i] - γ[i-1] * β[i-1]
    else
      β[i] = -γ[i-1] * β[i-1]
    end
    δ[i] = β[i] / α[i]
    α[n] -= δ[i] * β[i]
  end

  t = similar(a)
  t[1] = g[1]
  t[n] = g[n] - δ[1] * t[1]
  for i in 2 : n - 1
    t[i] = g[i] - γ[i-1] * t[i-1]
    t[n] -= δ[i] * t[i]
  end

  k = similar(a)
  k[n]   = t[n] / α[n]
  k[n-1] = t[n-1] / α[n-1] - δ[n-1] * k[n]
  for i in n - 2 : -1 : 1
    k[i] = t[i] / α[i] - γ[i] * k[i+1] - δ[i] * k[n]
  end

  k
end

# function build_spline(x, y, h, f, k, debug = true)
#   T = eltype(x)
#   n = size(h, 1)
#   w = Array{Array{T}}(n)
#   w[1] = Array{T}(4)
#   w[1][1] = y[1]
#   w[1][2] = k[n]
#   w[1][3] = (3f[1] - 2k[n] - k[1]) / h[1]
#   w[1][4] = (k[n] + k[1] - 2f[1]) / (h[1] * h[1])
#   for i in 2 : n
#     w[i] = Array{T}(4)
#     w[i][1] = y[i]
#     w[i][2] = k[i-1]
#     w[i][3] = (3f[i] - 2k[i-1] - k[i]) / h[i]
#     w[i][4] = (k[i-1] + k[i] - 2f[i]) / (h[i] * h[i])
#   end
#   s = t -> begin
#     if debug
#       v = pder
#     else
#       v = pval
#     end
#     for i in 1 : n
#       if t < x[i+1]
#         return v(w[i], t - x[i])
#       end
#     end
#     v(w[1], t - x[n])
#   end
#   s
# end



# function spline_cyclic(x, y; debug = false)
#   n = size(x,1) - 1
#   T = eltype(x)
#
#   # Tworzymy macierz dla naszego układu równań
#   a, b, c, g, h, f = build_matrix(x, y)
#   @show a
#   @show b
#   @show c
#   @show g
#
#
#   # Przekształcamy ją do postaci symetrycznej
#   a, b = reduce_matrix(a, b, c)
#   @show b
#
#
#   # Rozwiązujemy układ równań
#   k = solve_matrix(a, b, g)
#   k = solve_pivot(a, b, c, g)
#
#
#   s = build_spline(x, y, h, f, k)
#   s, k
# end



# function solve_pivot(a, b, c, g)
#   n = size(a, 1)
#   T = eltype(a)
#   A = zeros(T, n, n)
#   A[1,1] = a[1]
#   A[1,2] = b[2]
#   A[1,n] = c[n]
#   for i in 2 : n-1
#     A[i, i] = a[i]
#     A[i, i-1] = c[i-1]
#     A[i, i+1] = b[i]
#   end
#   A[n,n] = a[n]
#   A[n,1] = b[n]
#   A[n,n-1] = c[n-1]
#
#   k = gauss!(A,g)
# end
