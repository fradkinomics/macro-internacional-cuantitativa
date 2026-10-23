# Problem Set 2 — Problema de la torta
using Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()   # instala los paquetes si faltan
cd(@__DIR__)                        # trabaja en la carpeta del archivo

include("itpcake.jl")               # carga itpcake.jl, que a su vez incluye cakeeating.jl
const NK = 1000                     # grilla fina, única para todo el trabajo
mkpath("figures")                   # carpeta de gráficos




# ===================== 1.1 Funciones de consumo =====================

# Crea y resuelve un modelo; metodo = :vfi (discreto) o :itp (interpolación)
function resolver(; metodo = :vfi, Nk = NK, kw...)
    ce = CakeEating(; Nk = Nk, kw...)
    metodo == :vfi ? vfi!(ce) : vfi_itp!(ce)
    println()
    return ce
end

# Fracción de la torta que se come cada período según la forma cerrada (con r)
theta(; β = 0.96, r = 0.02, γ = 2) = (1 + r) - (β * (1 + r))^(1 / γ)

# Colores y orden en la leyenda, fijos por serie
const COLOR = Dict("Sin interpolación (vfi!)" => "#636efa", "Con interpolación (vfi_itp!)" => "#EF553B",
                   "Forma cerrada θ" => "#00cc96")
const RANGO = Dict("Sin interpolación (vfi!)" => 1, "Con interpolación (vfi_itp!)" => 2, "Forma cerrada θ" => 3)

# Gráfico de varias series contra la misma x (sin título). La última serie queda dibujada arriba.
function grafico(x, ys, nombres; xt = "k", yt = "")
    sc = [scatter(x = x, y = y, mode = "lines+markers", name = n, line_color = COLOR[n],
                  marker_color = COLOR[n], legendrank = RANGO[n]) for (y, n) in zip(ys, nombres)]
    return plot(sc, Layout(xaxis_title = xt, yaxis_title = yt, margin = attr(t = 20)))
end

# Imprime los números clave de un modelo resuelto
function resumen(ce, nombre; kmin = 0.2)
    ck = ce.gc ./ ce.kgrid
    m  = ce.kgrid .>= kmin
    println(rpad(nombre, 14), "| v(1) = ", round(ce.v[end], digits = 2),
            " | c(1) = ", round(ce.gc[end], digits = 4),
            " | c/k con k≥$kmin: mín ", round(minimum(ck[m]), digits = 4),
            ", máx ", round(maximum(ck[m]), digits = 4))
end



#BLOQUE C
ce20  = resolver(Nk = 20)                  # 1.1.1: parámetros por defecto, método discreto
ce20i = resolver(Nk = 20, metodo = :itp)   # mismo modelo, con interpolación
ceN   = resolver()                         # grilla fina, discreto
ceNi  = resolver(metodo = :itp)            # grilla fina, interpolación

println("θ teórico = ", round(theta(), digits = 4))
for (ce, n) in ((ce20, "Discreto 20"), (ce20i, "Interp. 20"), (ceN, "Discreto $NK"), (ceNi, "Interp. $NK"))
    resumen(ce, n)
end



#BLOQUE D
met = ["Sin interpolación (vfi!)", "Con interpolación (vfi_itp!)"]

# 1.1.2: consumo en niveles (objeto por defecto, 20 puntos)
savefig(grafico(ce20.kgrid, [ce20.gc, ce20i.gc], met, yt = "c(k)"), "figures/c_k.pdf")

# 1.1.4: ahorro relativo (objeto por defecto, 20 puntos)
savefig(grafico(ce20.kgrid, [ce20.gk ./ ce20.kgrid, ce20i.gk ./ ce20i.kgrid], met, yt = "k'/k"), "figures/kk_k.pdf")

# 1.1.3 + Modos de resolución: c/k (θ abajo, luego vfi!, luego vfi_itp! arriba) y v(k)
for (a, b, N) in ((ce20, ce20i, 20), (ceN, ceNi, NK))
    θs = fill(theta(), length(a.kgrid))
    savefig(grafico(a.kgrid, [θs, a.gc ./ a.kgrid, b.gc ./ b.kgrid], ["Forma cerrada θ"; met], yt = "c/k"),
            "figures/ck_k_N$N.pdf")
    m = a.kgrid .>= 0.2                                   # corta la zona donde v → -∞
    savefig(grafico(a.kgrid[m], [a.v[m], b.v[m]], met, yt = "v(k)"), "figures/v_k_N$N.pdf")
end




# ===================== Para pensar: opción 1 =====================

u0(c, γ) = c^(1 - γ) / (1 - γ)          # CRRA sin la constante -1 (la de la consigna)

# Encuentra ṽ iterando la ecuación de Bellman evaluada en k = 1
function vtilde(; β = 0.96, r = 0.02, γ = 2, tol = 1e-10, maxiter = 5000)
    vt = 0.0
    for it in 1:maxiter
        res   = Optim.maximize(kp -> u0(1 + r - kp, γ) + β * vt * kp^(1 - γ), 1e-6, 1 + r - 1e-10)
        nuevo = Optim.maximum(res)
        if abs(nuevo - vt) < tol * max(1, abs(vt))
            println("ṽ numérico = ", round(nuevo, digits = 4), " en $it iteraciones")
            return nuevo
        end
        vt = nuevo
    end
    error("ṽ no convergió")
end

# Regla de consumo eligiendo k' en forma continua con la v exacta, en los puntos de la grilla de ce
function politica_vt(ce, vt)
    β, r, γ = ce.pars[:β], ce.pars[:r], ce.pars[:γ]
    gc = similar(ce.kgrid)
    for (j, kv) in enumerate(ce.kgrid)
        res   = Optim.maximize(kp -> u(kv * (1 + r) - kp, γ) + β * vt * kp^(1 - γ),
                               1e-12 * kv, kv * (1 + r) * (1 - 1e-12))
        gc[j] = kv * (1 + r) - Optim.maximizer(res)
    end
    return gc
end

vt = vtilde()
θ  = theta(); g = 1 + 0.02 - θ
println("ṽ forma cerrada = ", round(-1 / (θ * (1 - 0.96 / g)), digits = 4))

gc20_vt = politica_vt(ce20, vt)
gcN_vt  = politica_vt(ceN,  vt)
for (gc, ce, N) in ((gc20_vt, ce20, 20), (gcN_vt, ceN, NK))
    ck = gc ./ ce.kgrid
    println("c/k con v exacta, N = $N: mín ", round(minimum(ck), digits = 6), ", máx ", round(maximum(ck), digits = 6))
end

# Figura: θ abajo, luego vfi!, vfi_itp! y v exacta arriba
COLOR["Con v exacta"] = "#AB63FA"; RANGO["Con v exacta"] = 3; RANGO["Forma cerrada θ"] = 4
savefig(grafico(ce20.kgrid, [fill(θ, 20), ce20.gc ./ ce20.kgrid, ce20i.gc ./ ce20i.kgrid, gc20_vt ./ ce20.kgrid],
                ["Forma cerrada θ"; met; "Con v exacta"], yt = "c/k"), "figures/ck_vt.pdf")

                


# ===================== Para pensar: opción 2 (x = log k) =====================

struct CakeLog
    pars::Dict{Symbol, Float64}
    xgrid::Vector{Float64}       # grilla equiespaciada en x = log k
    kgrid::Vector{Float64}       # los mismos puntos en niveles, k = exp(x)
    v::Vector{Float64}
    gc::Vector{Float64}
    gk::Vector{Float64}          # torta de mañana en niveles, exp(x')
end

function CakeLog(; β = 0.96, r = 0.02, γ = 2, Nx = 20, xmin = -10.0, xmax = 0.0)
    xgrid = collect(range(xmin, xmax, length = Nx))
    pars  = Dict(:β => β, :r => r, :γ => γ)
    return CakeLog(pars, xgrid, exp.(xgrid), zeros(Nx), zeros(Nx), zeros(Nx))
end

# Una iteración de Bellman en x. Con extrap = true, x' puede salir de la grilla y v se prolonga con una recta
function vfi_log_iter!(new_v, cl::CakeLog, extrap::Bool)
    β, r, γ = cl.pars[:β], cl.pars[:r], cl.pars[:γ]
    itp   = interpolate((cl.xgrid,), cl.v, Gridded(Linear()))
    itp_v = extrap ? extrapolate(itp, Line()) : itp
    for (j, xv) in enumerate(cl.xgrid)
        x_hi = xv + log(1 + r) + log(1 - 1e-10)                      # deja consumo positivo
        x_lo = extrap ? first(cl.xgrid) - 5 : first(cl.xgrid)         # con extrapolación, hasta 5 por debajo
        extrap || (x_hi = min(x_hi, last(cl.xgrid)))
        obj = xp -> u(exp(xv) * (1 + r) - exp(xp), γ) + β * itp_v(xp)
        res = Optim.maximize(obj, x_lo, x_hi)
        xp  = Optim.maximizer(res)
        new_v[j] = Optim.maximum(res)
        cl.gk[j] = exp(xp)
        cl.gc[j] = exp(xv) * (1 + r) - exp(xp)
    end
end

function vfi_log!(cl::CakeLog; extrap = false, tol = 1e-8, maxiter = 2000)
    new_v = similar(cl.v); dist, iter = 1 + tol, 0
    while dist > tol && iter < maxiter
        iter += 1
        vfi_log_iter!(new_v, cl, extrap)
        dist = norm(new_v - cl.v) / max(1, norm(cl.v))
        cl.v .= new_v
    end
    println("log k, N = $(length(cl.xgrid)), extrap = $extrap: $iter iteraciones. Distance = $dist")
    return cl
end

logs = Dict((N, e) => vfi_log!(CakeLog(Nx = N); extrap = e) for N in (20, NK) for e in (false, true))

# Resumen: valor en k = 1, c/k en los extremos, desvío lejos del cero y zona donde c/k se aparta más de 1% de θ
θ = theta()
zona(k, ck; tol = 0.01) = (m = abs.(ck ./ θ .- 1) .> tol; any(m) ? maximum(k[m]) : 0.0)
println(rpad("Interp. en k, N=$NK", 30), "| zona >1%: k < ", round(zona(ceNi.kgrid, ceNi.gc ./ ceNi.kgrid), sigdigits = 2))
for N in (20, NK), e in (false, true)
    cl = logs[(N, e)]; ck = cl.gc ./ cl.kgrid; m = cl.kgrid .>= 0.2
    println(rpad("log k, N=$N, " * (e ? "con" : "sin") * " extrap.", 30),
            "| v(1) = ", round(cl.v[end], digits = 2),
            " | c/k(1) = ", round(ck[end], digits = 5),
            " | c/k(k mín) = ", round(ck[1], digits = 4),
            " | desvío k≥0.2: ", round(100 * maximum(abs.(ck[m] ./ θ .- 1)), digits = 2), "%",
            " | zona >1%: k < ", round(zona(cl.kgrid, ck), sigdigits = 2))
end

# Gráfico con una x distinta para cada serie (las grillas en k y en log k no coinciden)
function grafico_xs(xs, ys, nombres; xt = "k", yt = "", xtipo = "linear")
    sc = [scatter(x = x, y = y, mode = "lines+markers", name = n, line_color = COLOR[n],
                  marker_color = COLOR[n], legendrank = RANGO[n]) for (x, y, n) in zip(xs, ys, nombres)]
    return plot(sc, Layout(xaxis_title = xt, yaxis_title = yt, xaxis_type = xtipo, margin = attr(t = 20)))
end

COLOR["log k, sin extrapolación"] = "#FFA15A"; RANGO["log k, sin extrapolación"] = 3
COLOR["log k, con extrapolación"] = "#AB63FA"; RANGO["log k, con extrapolación"] = 4
RANGO["Forma cerrada θ"] = 5
nombres = ["Forma cerrada θ", "Con interpolación (vfi_itp!)", "log k, sin extrapolación", "log k, con extrapolación"]
for (N, ci) in ((20, ce20i), (NK, ceNi))
    a, b = logs[(N, false)], logs[(N, true)]
    xs = [ci.kgrid, ci.kgrid, a.kgrid, b.kgrid]
    ys = [fill(θ, length(ci.kgrid)), ci.gc ./ ci.kgrid, a.gc ./ a.kgrid, b.gc ./ b.kgrid]
    savefig(grafico_xs(xs, ys, nombres, yt = "c/k"), "figures/ck_log_N$N.pdf")
    N == NK && savefig(grafico_xs(xs, ys, nombres, yt = "c/k", xtipo = "log"), "figures/ck_log_N$(N)_logx.pdf")
end




# ===================== 1.2 Simulador =====================
const T  = 250          # horizonte, como en el gráfico de clase
const K0 = 1.0          # torta inicial = máximo de la grilla

# Pasos 2–4: simula {c_t, k_t} para t = 0,…,T con las reglas de un modelo resuelto
function simular(ce::CakeEating, T, k0)
    k0 <= maximum(ce.kgrid) || error("k0 fuera de la grilla")
    𝒞, 𝒦 = zeros(T + 1), zeros(T + 1)                              # posición t+1 = período t
    itp_c = interpolate((ce.kgrid,), ce.gc, Gridded(Linear()))      # regla de consumo c(k)
    itp_k = interpolate((ce.kgrid,), ce.gk, Gridded(Linear()))      # regla de ahorro k'(k)
    𝒦[1] = k0
    for t in 1:T+1
        𝒞[t] = itp_c(𝒦[t])                                         # consumo de hoy
        t <= T && (𝒦[t+1] = itp_k(𝒦[t]))                           # torta de mañana
    end
    return 𝒞, 𝒦
end

# Trayectoria de la forma cerrada con interés: k_t = k0 g^t, c_t = θ k_t
function forma_cerrada(ce::CakeEating, T, k0)
    β, r, γ = ce.pars[:β], ce.pars[:r], ce.pars[:γ]
    θ = theta(β = β, r = r, γ = γ); g = 1 + r - θ
    𝒦 = [k0 * g^t for t in 0:T]
    return θ .* 𝒦, 𝒦
end

𝒞c, 𝒦c = forma_cerrada(ceNi, T, K0)
modelos = ((ce20, "Discreto 20"), (ce20i, "Interp. 20"), (ceN, "Discreto $NK"), (ceNi, "Interp. $NK"))
sims = Dict(n => simular(ce, T, K0) for (ce, n) in modelos)

println(rpad("Forma cerrada", 14), "| c_0 = ", round(𝒞c[1], digits = 4), " | c_T = ", round(𝒞c[end], digits = 4),
        " | k_T = ", round(𝒦c[end], digits = 4), " | Σc = ", round(sum(𝒞c), digits = 3))
for (_, n) in modelos
    𝒞, 𝒦 = sims[n]
    e = abs.(𝒞 ./ 𝒞c .- 1) .* 100                                  # error relativo en %
    println(rpad(n, 14), "| c_0 = ", round(𝒞[1], digits = 4), " | c_T = ", round(𝒞[end], digits = 4),
            " | k_T = ", round(𝒦[end], digits = 4), " | error c_t: máx ", round(maximum(e), digits = 2),
            "%, medio ", round(sum(e) / length(e), digits = 2), "%")
end

# Figuras: forma cerrada abajo, luego vfi!, vfi_itp! arriba
COLOR["Forma cerrada"] = "#00cc96"; RANGO["Forma cerrada"] = 3
tt = collect(0:T)
savefig(grafico(tt, [𝒞c, sims["Discreto $NK"][1], sims["Interp. $NK"][1]], ["Forma cerrada"; met],
                xt = "t", yt = "c<sub>t</sub>"), "figures/sim_c.pdf")
savefig(grafico(tt, [𝒦c, sims["Discreto $NK"][2], sims["Interp. $NK"][2]], ["Forma cerrada"; met],
                xt = "t", yt = "k<sub>t</sub>"), "figures/sim_k.pdf")




                # Extra: error relativo con signo (%) de c_t y k_t respecto de la forma cerrada, N = NK
err(x, xc) = 100 .* (x ./ xc .- 1)
savefig(grafico(tt, [err(sims["Discreto $NK"][1], 𝒞c), err(sims["Interp. $NK"][1], 𝒞c)], met,
                xt = "t", yt = "Error en c<sub>t</sub> (%)"), "figures/sim_err_c.pdf")
savefig(grafico(tt, [err(sims["Discreto $NK"][2], 𝒦c), err(sims["Interp. $NK"][2], 𝒦c)], met,
                xt = "t", yt = "Error en k<sub>t</sub> (%)"), "figures/sim_err_k.pdf")

# Error en la torta final, para el cuadro
for (_, n) in modelos
    println(rpad(n, 14), "| error en k_T = ", round(err(sims[n][2][end], 𝒦c[end]), digits = 2), "%")
end