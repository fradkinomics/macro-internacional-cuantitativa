# Problem Set 1 — Macroeconomía Internacional Cuantitativa
using Pkg; Pkg.activate(@__DIR__)   # usa los paquetes de esta carpeta (Project.toml)
cd(@__DIR__)                        # trabaja en la carpeta del archivo

include("mccall.jl")                # carga el código base del profesor
using Statistics                    # mean y quantile
using Random; Random.seed!(1234)    # fija el azar: las simulaciones dan siempre lo mismo
const NW = 2000    # resolución de la grilla para todo el trabajo

mc = McCall()                       # modelo con los parámetros por defecto
vfi!(mc, verbose = false)           # lo resuelve
mc.w_star                           # muestra el salario de reserva



# ===================== 1.1 Estáticas comparadas =====================

# Resuelve el modelo para cada valor del parámetro y guarda w*
function estatica(param, valores; Nw = NW)
    V = similar(valores)
    for (j, x) in enumerate(valores)
        mc = param == :b ? McCall(b = x, Nw = Nw) : McCall(β = x, Nw = Nw)
        vfi!(mc, verbose = false)
        V[j] = mc.w_star
    end
    return V
end

# Gráfico de puntos
  grafico(x, y, xt, tit) = plot(scatter(x = x, y = y, mode = "markers"),
                                Layout(xaxis_title = xt, yaxis_title = "w*", margin = attr(t = 20)))

mkpath("figures")                               # crea la carpeta de gráficos
ℬb = collect(range(0.5, 1.5, length = 41))      # valores de b
ℬβ = collect(range(0.8, 0.99, length = 41))     # valores de β

for Nw in (50, NW)
    Vb = estatica(:b, ℬb; Nw = Nw)
    savefig(grafico(ℬb, Vb, "b", "Salario de reserva vs b (Nw = $Nw)"), "figures/wstar_b_N$Nw.pdf")
    Vβ = estatica(:β, ℬβ; Nw = Nw)
    savefig(grafico(ℬβ, Vβ, "β", "Salario de reserva vs β (Nw = $Nw)"), "figures/wstar_beta_N$Nw.pdf")
    println("Nw = $Nw | w*(b=0.5) = $(round(Vb[1], digits=3)), w*(b=1.5) = $(round(Vb[end], digits=3)) | w*(β=0.8) = $(round(Vβ[1], digits=3)), w*(β=0.99) = $(round(Vβ[end], digits=3))")
end



# ===================== 1.2 Simulaciones =====================

# Repite la simulación K veces y devuelve el vector de T's
function simular_T(mc::McCall, K::Int)
    𝒯 = Vector{Int}(undef, K)
    for k in 1:K
        𝒯[k] = simul(mc, verbose = false)
    end
    return 𝒯
end

mc = McCall(Nw = NW)
vfi!(mc, verbose = false)
𝒯 = simular_T(mc, 10_000)

p = sum(mc.pw[mc.wgrid .>= mc.w_star])      # prob. de aceptar en cada período
println("w* = ", round(mc.w_star, digits = 3), " | p = ", round(p, digits = 4))
println("E[T] simulado = ", mean(𝒯), " | teórico 1/p = ", round(1/p, digits = 3))
println("Cuantiles 25/50/75/90/99: ", [quantile(𝒯, q) for q in (0.25, 0.5, 0.75, 0.9, 0.99)])
println("Máximo: ", maximum(𝒯))

h = plot(histogram(x = 𝒯), Layout(xaxis_title = "T", yaxis_title = "Frecuencia"))
savefig(h, "figures/hist_T.pdf")
h



# ===================== 1.3 E[T] en función de β =====================

ℬ = collect(range(0.9, 0.99, length = 25))
ET_sim = similar(ℬ)      # media simulada
ET_teo = similar(ℬ)      # 1/p teórico
for (j, x) in enumerate(ℬ)
    m = McCall(β = x, Nw = NW)
    vfi!(m, verbose = false)
    ET_sim[j] = mean(simular_T(m, 10_000))
    ET_teo[j] = 1 / sum(m.pw[m.wgrid .>= m.w_star])
end

p3 = plot([scatter(x = ℬ, y = ET_sim, mode = "markers", name = "Simulado"),
           scatter(x = ℬ, y = ET_teo, mode = "lines", name = "Teórico 1/p")],
          Layout(xaxis_title = "β", yaxis_title = "E[T]"))
savefig(p3, "figures/ET_beta.pdf")
println("E[T]: β=0.90 → ", round(ET_sim[1], digits=2), " | β=0.99 → ", round(ET_sim[end], digits=2))
p3