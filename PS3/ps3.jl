# Problem Set 3 — Deuda con default (Arellano)
using Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()   # instala los paquetes si faltan
cd(@__DIR__)                        # trabaja en la carpeta del archivo

include("simuldef.jl")              # carga simuldef.jl, que a su vez incluye arellano.jl y deuda.jl
mkpath("figures")                   # carpeta de gráficos


# ===================== Constantes de todo el trabajo =====================
const SEMILLA = 1234                                   # semilla de las simulaciones
const TOL     = 1e-6                                   # tolerancia de convergencia (la de mpe!)
const MAXITER = 1000                                   # tope de iteraciones
const BETA0, YH0, CHI0 = 0.953, 0.969, 0.01            # parámetros por defecto del constructor
const BETAS  = [0.90, 0.92, 0.94, 0.953, 0.96, 0.97]   # paciencia (incluye el valor base)
const YHS    = [0.93, 0.95, 0.969, 0.98, 0.99]         # costo de default original (incluye el valor base)
const CHIS   = [0.0001, 0.1]                           # χ que fija la consigna en 1.2
const B_MAX_12 = 0.4                                   # deuda máxima que se muestra en los gráficos de 1.2
const BMAX_ROB = 1.2                                   # grilla de deuda más ancha para el chequeo de 1.3
const T_SIM, T_BURN, N_SIM = 10_000, 1_000, 400        # largo de cada simulación, períodos descartados y repeticiones
const UMBRAL = 0.5                                     # probabilidad que define b*(y)
const ROJO, AZUL, VERDE = "#EF553B", "#636efa", "#00cc96"
const RAMPA  = ["#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"]

Random.seed!(SEMILLA)

# Todo lo que se imprime queda también en resultados.txt
LOG = open("resultados.txt", "w")
informar(xs...) = (println(xs...); println(LOG, xs...); flush(LOG))
r4(x) = round(x, digits = 4)


# ===================== Funciones comunes =====================

# Igual que mpe! del curso, pero con la opción de no actualizar el precio de la deuda (para 1.3).
# Devuelve la cantidad de iteraciones y la distancia final.
function resolver!(dd::Default; actualizar_q = true, tol = TOL, maxiter = MAXITER)
    new_v, new_q = similar(dd.v), copy(dd.q)
    dist, iter = 1 + tol, 0
    while dist > tol && iter < maxiter
        iter += 1
        actualizar_q && q_iter!(new_q, dd)                       # precio dada la probabilidad de default
        dist_q = norm(new_q - dd.q) / max(1, norm(dd.q))
        itp_q  = make_itp(dd, new_q)
        vfi_iter!(new_v, itp_q, dd)                              # valor dado el precio
        dist_v = norm(new_v - dd.v) / max(1, norm(dd.v))
        dist   = max(dist_q, dist_v)
        dd.v .= new_v
        dd.q .= new_q
    end
    return iter, dist
end

# Elección de b' robusta para los modelos con default. El opt_value del curso busca un máximo con Optim sobre
# toda la grilla; cuando el objetivo tiene más de un pico (pasa con χ muy chico y con precios fijos) salta de un
# pico a otro entre iteraciones y la función de valor no converge. Con ROBUSTO[] = true primero se evalúa el
# objetivo en una grilla gruesa (uno de cada PASO_ROB puntos de bgrid) y después se refina con Optim alrededor
# del mejor punto. Con ROBUSTO[] = false se usa el opt_value del curso sin cambios.
const ROBUSTO  = Ref(false)
const PASO_ROB = 4
function opt_value(jb, jy, q, itp_v, dd::Arellano)
    ROBUSTO[] || return invoke(opt_value, Tuple{Any, Any, Any, Any, Deuda}, jb, jy, q, itp_v, dd)
    obj_f(bpv) = eval_value(jb, jy, bpv, q, itp_v, dd)[1]
    Nb = length(dd.bgrid)
    jmax, fmax = 1, -Inf
    for j in [1:PASO_ROB:Nb; Nb]
        f = obj_f(dd.bgrid[j])
        f > fmax && ((jmax, fmax) = (j, f))
    end
    lo, hi = dd.bgrid[max(jmax - PASO_ROB, 1)], dd.bgrid[min(jmax + PASO_ROB, Nb)]
    b_star = Optim.maximizer(Optim.maximize(obj_f, lo, hi))
    fmax > obj_f(b_star) && (b_star = dd.bgrid[jmax])       # Optim no evalúa los extremos del intervalo
    vp, c_star = eval_value(jb, jy, b_star, q, itp_v, dd)
    return vp, c_star, b_star
end

# Crea y resuelve un modelo de Arellano con costos originales
function modelo(nombre; actualizar_q = true, robusto = false, kw...)
    dd = Arellano(; kw...)
    t0 = time()
    ROBUSTO[] = robusto
    iter, dist = resolver!(dd; actualizar_q = actualizar_q)
    ROBUSTO[] = false
    informar(rpad(nombre, 34), "| ", iter, " iteraciones | dist = ", round(dist, sigdigits = 3),
             " | ", round(time() - t0, digits = 1), " s")
    return dd
end

# b*(y): nivel de deuda donde la probabilidad de default cruza el umbral (interpolando entre puntos de la grilla).
# NaN si la probabilidad no llega al umbral dentro de la grilla.
function bstar(dd; umbral = UMBRAL)
    bs = fill(NaN, length(dd.ygrid))
    for jy in eachindex(dd.ygrid)
        p  = dd.prob[:, jy]
        jb = findfirst(>=(umbral), p)
        if jb === nothing
            continue
        elseif jb == 1
            bs[jy] = dd.bgrid[1]
        else
            w = (umbral - p[jb-1]) / (p[jb] - p[jb-1])
            bs[jy] = dd.bgrid[jb-1] + w * (dd.bgrid[jb] - dd.bgrid[jb-1])
        end
    end
    return bs
end

# Gráfico de varias series (sin título). La serie número `arriba` se dibuja por encima de las demás.
# Con chico = true usa letra grande y la leyenda debajo, para las figuras que van de a tres por fila.
function grafico(xs, ys, nombres, colores; xt = "", yt = "", arriba = length(xs), chico = false)
    sc = [scatter(x = x, y = y, mode = "lines+markers", name = n, line_color = c, marker_color = c, legendrank = k)
          for (k, (x, y, n, c)) in enumerate(zip(xs, ys, nombres, colores))]
    orden = [setdiff(1:length(sc), arriba); arriba]
    lay = chico ? Layout(xaxis_title = xt, yaxis_title = yt, margin = attr(t = 20), font_size = 22,
                         legend = attr(orientation = "h", x = 0, y = -0.3))  :
                  Layout(xaxis_title = xt, yaxis_title = yt, margin = attr(t = 20))
    return plot(sc[orden], lay)
end
guardar_chico(p, archivo) = savefig(p, archivo, width = 520, height = 520)

# Contour de la probabilidad de default, con la curva b*(y) encima.
# En contour las filas de z van al eje vertical: prob[jb, jy] => x = ygrid, y = bgrid.
function grafico_contour(dd, archivo)
    bs = bstar(dd); m = isfinite.(bs)
    tr = [contour(x = dd.ygrid, y = dd.bgrid, z = dd.prob, colorbar_title = "P(b, y)")]
    any(m) && push!(tr, scatter(x = dd.ygrid[m], y = bs[m], mode = "lines+markers", name = "b*(y)",
                                line_color = "white", marker_color = "white", showlegend = false))
    savefig(plot(tr, Layout(xaxis_title = "y", yaxis_title = "b", margin = attr(t = 20))), archivo)
end

# Gráfico de b*(y) para una familia de modelos que difieren en un parámetro
# Con dif = true grafica la diferencia con el b*(y) del modelo base (sirve cuando las curvas se tapan).
function grafico_bstar(mods, valores, base, simbolo, archivo; dif = false)
    xs, ys, ns, cs = [], [], [], []
    bs0 = bstar(mods[base])
    k = 0
    for val in valores
        dd = mods[val]; bs = bstar(dd)
        if val == base
            dif && continue
            m = isfinite.(bs)
            push!(xs, dd.ygrid[m]); push!(ys, bs[m]); push!(ns, "$simbolo = $val (base)"); push!(cs, ROJO)
        else
            k += 1
            m = dif ? isfinite.(bs .- bs0) : isfinite.(bs)
            push!(xs, dd.ygrid[m]); push!(ys, dif ? (bs .- bs0)[m] : bs[m]); push!(ns, "$simbolo = $val"); push!(cs, RAMPA[k])
        end
    end
    arriba = dif ? length(xs) : findfirst(==(base), valores)
    yt = dif ? "b*(y) − b*(y) con $simbolo = $base" : "b*(y)"
    savefig(grafico(xs, ys, ns, cs, xt = "y", yt = yt, arriba = arriba), archivo)
end

fila(v) = join([lpad(isnan(x) ? "—" : string(round(x, digits = 3)), 7) for x in v])


# ===================== 1.1 Regiones de repago y default =====================
informar("===== 1.1 Regiones de repago y default =====")
dd0 = modelo("Base")
JY = argmin(abs.(dd0.ygrid .- 1))                      # índice del ingreso más cercano a 1
informar("ygrid: ", length(dd0.ygrid), " puntos de ", r4(dd0.ygrid[1]), " a ", r4(dd0.ygrid[end]),
         " | y[", JY, "] = ", r4(dd0.ygrid[JY]))
informar("bgrid: ", length(dd0.bgrid), " puntos de ", dd0.bgrid[1], " a ", dd0.bgrid[end],
         " | paso = ", r4(dd0.bgrid[2] - dd0.bgrid[1]))

grafico_contour(dd0, "figures/prob_contour.pdf")

jys = [1, 6, 11, 16, 21]
informar("P(b,y) base. Columnas: y = ", fila(dd0.ygrid[jys]))
for bv in 0:0.05:0.6
    jb = argmin(abs.(dd0.bgrid .- bv))
    informar("  b = ", rpad(round(dd0.bgrid[jb], digits = 3), 6), fila(dd0.prob[jb, jys]))
end
informar("P creciente en b para todo y: ", all(diff(dd0.prob, dims = 1) .>= -1e-8),
         " | P decreciente en y para todo b: ", all(diff(dd0.prob, dims = 2) .<= 1e-8))
informar("máx de P que sube con y: ", round(maximum(diff(dd0.prob, dims = 2)), sigdigits = 3))
informar("y      = ", fila(dd0.ygrid))
informar("b*(y)  = ", fila(bstar(dd0)))
informar("q(b',y) en y[", JY, "]: ", fila([dd0.q[argmin(abs.(dd0.bgrid .- bv)), JY] for bv in 0:0.1:0.6]))

# Opcional: b*(y) para distintos β y distintos costos de default
informar("--- Opcional 1.1: variar β ---")
mods_beta = Dict(β => (β == BETA0 ? dd0 : modelo("β = $β"; β = β)) for β in BETAS)
informar("--- Opcional 1.1: variar yh ---")
mods_yh = Dict(yh => (yh == YH0 ? dd0 : modelo("yh = $yh"; yh = yh)) for yh in YHS)

informar("b*(y) por β.   y = ", fila(dd0.ygrid))
for β in BETAS
    informar("  β = ", rpad(β, 13), fila(bstar(mods_beta[β])))
end
informar("b*(y) por yh.  y = ", fila(dd0.ygrid))
for yh in YHS
    informar("  yh = ", rpad(yh, 12), fila(bstar(mods_yh[yh])))
end
grafico_bstar(mods_beta, BETAS, BETA0, "β", "figures/bstar_beta.pdf")                     # niveles (las curvas se tapan)
grafico_bstar(mods_beta, BETAS, BETA0, "β", "figures/bstar_beta_dif.pdf", dif = true)     # diferencia con el base
grafico_bstar(mods_yh, YHS, YH0, "yh", "figures/bstar_yh.pdf")


# ===================== 1.2 Mecánica de los shocks de preferencias =====================
informar("===== 1.2 Shocks de preferencias =====")
mods_chi = Dict(χ => modelo("χ = $χ"; χ = χ, robusto = (χ < CHI0)) for χ in CHIS)
mods_chi[CHI0] = dd0

for χ in sort(collect(keys(mods_chi)))
    dd = mods_chi[χ]
    vR, vD, V = dd.vR[:, JY], dd.vD[JY], dd.v[:, JY]
    jb = findfirst(vR .< vD)                               # primer punto donde conviene defaultear
    informar("χ = ", χ, " | vD(y=1) = ", r4(vD), " | vR(0) = ", r4(vR[1]), " | vR(bmax) = ", r4(vR[end]),
             " | V(0) = ", r4(V[1]), " | V(bmax) = ", r4(V[end]))
    if jb === nothing || jb == 1
        informar("   vR y vD no se cruzan dentro de la grilla")
    else
        w  = (vR[jb-1] - vD) / (vR[jb-1] - vR[jb])
        bc = dd.bgrid[jb-1] + w * (dd.bgrid[jb] - dd.bgrid[jb-1])
        informar("   cruce vR = vD en b = ", r4(bc), " | V - max(vR, vD) en el punto de grilla anterior al cruce = ",
                 r4(V[jb-1] - max(vR[jb-1], vD)), " | χ·log 2 = ", r4(χ * log(2)))
    end
    informar("   máx de V - max(vR, vD) en toda la grilla de b = ", round(maximum(V .- max.(vR, vD)), sigdigits = 3))
    for bv in 0:0.1:0.6
        k = argmin(abs.(dd.bgrid .- bv))
        informar("   b = ", rpad(round(dd.bgrid[k], digits = 3), 6), "| vR = ", r4(vR[k]), " | V = ", r4(V[k]),
                 " | P = ", r4(dd.prob[k, JY]))
    end
end

# Un gráfico por χ: vR y vD abajo, 𝒱 arriba
for (χ, nombre) in zip(CHIS, ("valores_chi_bajo", "valores_chi_alto"))
    dd = mods_chi[χ]; m = dd.bgrid .<= B_MAX_12; b = dd.bgrid[m]
    savefig(grafico([b, b, b], [dd.vR[m, JY], fill(dd.vD[JY], length(b)), dd.v[m, JY]],
                    ["v<sup>R</sup>(b, y)", "v<sup>D</sup>(y)", "V(b, y)"], [AZUL, VERDE, ROJO],
                    xt = "b", yt = "valor"), "figures/$nombre.pdf")
end


# ===================== 1.3 Cuántos defaults son por los precios? =====================
informar("===== 1.3 Defaults por precios =====")
# El constructor arranca con q = 1/(1+r) en toda la grilla: el precio de una deuda que se repaga siempre.
# Sin actualizar q, el gobierno resuelve su problema (con opción de default) frente a ese precio.
ddq = modelo("Precios de repago seguro"; actualizar_q = false, robusto = true)
informar("q fijo = ", r4(ddq.q[1, 1]), " en toda la grilla: ", all(ddq.q .≈ 1 / (1 + ddq.pars.r)))
informar("P con precios fijos: mín = ", round(minimum(ddq.prob), sigdigits = 3), " | máx = ", round(maximum(ddq.prob), sigdigits = 3))
informar("y                        = ", fila(dd0.ygrid))
informar("b*(y) equilibrio         = ", fila(bstar(dd0)))
informar("b*(y) precios fijos      = ", fila(bstar(ddq)))
D0, Dq = dd0.prob .> UMBRAL, ddq.prob .> UMBRAL
informar("Puntos (b,y) con P > 0.5: equilibrio = ", sum(D0), " | precios fijos = ", sum(Dq), " | de ", length(D0))
informar("   sólo en equilibrio = ", sum(D0 .& .!Dq), " | sólo con precios fijos = ", sum(Dq .& .!D0), " | en ambos = ", sum(D0 .& Dq))
informar("P_fijos - P_equilibrio: mín = ", round(minimum(ddq.prob - dd0.prob), sigdigits = 3),
         " | máx = ", round(maximum(ddq.prob - dd0.prob), sigdigits = 3))
informar("vD(y) equilibrio    = ", fila(dd0.vD))
informar("vD(y) precios fijos = ", fila(ddq.vD))
informar("b'(b, y=1) equilibrio:    ", fila([dd0.gb[argmin(abs.(dd0.bgrid .- bv)), JY] for bv in 0:0.1:0.6]))
informar("b'(b, y=1) precios fijos: ", fila([ddq.gb[argmin(abs.(ddq.bgrid .- bv)), JY] for bv in 0:0.1:0.6]))
informar("P(b, y) precios fijos. Columnas: y = ", fila(ddq.ygrid[jys]))
for bv in 0:0.1:0.6
    jb = argmin(abs.(ddq.bgrid .- bv))
    informar("  b = ", rpad(round(ddq.bgrid[jb], digits = 3), 6), fila(ddq.prob[jb, jys]))
end

grafico_contour(ddq, "figures/prob_contour_precios.pdf")
m0, mq = isfinite.(bstar(dd0)), isfinite.(bstar(ddq))
savefig(grafico([ddq.ygrid[mq], dd0.ygrid[m0]], [bstar(ddq)[mq], bstar(dd0)[m0]],
                ["Precios de repago seguro", "Precios de equilibrio"], [AZUL, ROJO], xt = "y", yt = "b*(y)"),
        "figures/bstar_precios.pdf")

# Chequeo: con precios fijos el gobierno se endeuda hasta el tope de la grilla, así que b*(y) puede depender de bmax
ddq2 = modelo("Precios fijos, bmax = $BMAX_ROB"; actualizar_q = false, robusto = true, bmax = BMAX_ROB)
informar("b*(y) precios fijos, bmax = ", BMAX_ROB, ": ", fila(bstar(ddq2)))
informar("b'(b, y=1) precios fijos, bmax = ", BMAX_ROB, ": ", fila([ddq2.gb[argmin(abs.(ddq2.bgrid .- bv)), JY] for bv in 0:0.2:1.2]))

informar("===== FIN EJERCICIO 1 =====")


# ===================== 2. Calibración del modelo de default =====================
informar("===== 2 Calibración =====")

# Momentos del curso (deuda/PBI, spread promedio y desvío) más la mediana del spread, la frecuencia anual de
# default (defaults nuevos por cada 100 años con acceso al mercado) y el % del tiempo en default
function momentos(pp::SimulPath)
    M = moments(pp)
    d = pp[:d]
    rep_indices = [t for t in eachindex(pp) if pp[:d, t] == 0]
    nuevos = count(t -> d[t] == 1 && d[t-1] == 0, 2:length(d))
    M[:medspr]  = median(pp[:spr][rep_indices]) * 10_000
    M[:deffreq] = 100 * 4 * nuevos / max(1, count(==(0), d[1:end-1]))
    M[:tdef]    = 100 * mean(d)
    return M
end

# Para un vector de simulaciones: calcula los momentos de cada una y los promedia
function momentos(pps::Vector{SimulPath})
    Ms = [momentos(pp) for pp in pps]
    return Dict(k => mean(M[k] for M in Ms) for k in keys(Ms[1]))
end

# Para un modelo resuelto: N_SIM simulaciones de T_SIM períodos, descartando los primeros T_BURN de cada una
function momentos(dd::Default; semilla = SEMILLA)
    Random.seed!(semilla)
    return momentos([subpath(simul(dd, T = T_SIM + T_BURN), T_BURN + 1, T_SIM + T_BURN) for _ in 1:N_SIM])
end

const CLAVES = [:debtGDP, :avgspr, :stdspr, :medspr, :deffreq, :tdef]

# Momentos de una familia de modelos, con dos semillas distintas para chequear que no dependen del azar
function tabla_momentos(mods, valores, simbolo)
    informar(rpad(simbolo, 8), join(lpad.(string.(CLAVES), 11)), "   | dif. máx. entre semillas (%)")
    Ms = []
    for val in valores
        t0 = time()
        M  = momentos(mods[val])
        M2 = momentos(mods[val], semilla = SEMILLA + 1)
        dif = maximum(abs(M2[k] / M[k] - 1) for k in CLAVES[1:3]) * 100
        informar(rpad(val, 8), join([lpad(round(M[k], digits = 2), 11) for k in CLAVES]), "   | ", round(dif, digits = 2),
                 "  (", round(time() - t0, digits = 1), " s)")
        push!(Ms, M)
    end
    return Ms
end

function grafico_momentos(valores, Ms, simbolo, sufijo)
    x = valores
    col(k) = [M[k] for M in Ms]
    guardar_chico(grafico([x], [col(:debtGDP)], ["Deuda / PBI"], [AZUL], xt = simbolo, yt = "Deuda / PBI anual (%)",
                          chico = true), "figures/mom_deuda_$sufijo.pdf")
    guardar_chico(grafico([x, x], [col(:stdspr), col(:avgspr)], ["Desvío", "Promedio"], [AZUL, ROJO],
                          xt = simbolo, yt = "Spread (pb)", chico = true), "figures/mom_spread_$sufijo.pdf")
    guardar_chico(grafico([x], [col(:deffreq)], ["Frecuencia de default"], [VERDE], xt = simbolo,
                          yt = "Defaults por año (%)", chico = true), "figures/mom_def_$sufijo.pdf")
end

informar("--- Momentos en función de β (yh = ", YH0, ") ---")
Ms_beta = tabla_momentos(mods_beta, BETAS, "β")
grafico_momentos(BETAS, Ms_beta, "β", "beta")

informar("--- Momentos en función de yh (β = ", BETA0, ") ---")
Ms_yh = tabla_momentos(mods_yh, YHS, "yh")
grafico_momentos(YHS, Ms_yh, "yh", "yh")

# Chequeo: cuánto de los spreads y de la frecuencia de default viene del shock de preferencias
informar("--- Momentos en función de χ (β = ", BETA0, ", yh = ", YH0, ") ---")
CHIS_TODOS = sort(collect(keys(mods_chi)))
Ms_chi = tabla_momentos(mods_chi, CHIS_TODOS, "χ")
for χ in CHIS_TODOS
    dd = mods_chi[χ]
    informar("χ = ", rpad(χ, 7), "| P(b = 0, y) con y = ", fila(dd.ygrid[[1, 6, 11]]), " : ", fila(dd.prob[1, [1, 6, 11]]),
             " | b*(y) en los mismos y: ", fila(bstar(dd)[[1, 6, 11]]))
end

# Diagnóstico de una simulación del modelo base
pp = simul(dd0, T = T_SIM)
rep = [t for t in eachindex(pp) if pp[:d, t] == 0]
informar("Simulación base (T = ", T_SIM, "): spread máx = ", round(maximum(pp[:spr][rep]) * 1e4, digits = 1),
         " pb | spread no finito: ", count(!isfinite, pp[:spr][rep]),
         " | b/y medio en repago = ", r4(mean(pp[:b][rep] ./ pp[:y][rep])), " | b máx = ", r4(maximum(pp[:b])),
         " | y medio en default = ", r4(mean(pp[:y][[t for t in eachindex(pp) if pp[:d, t] == 1]])),
         " | y medio en repago = ", r4(mean(pp[:y][rep])))

informar("===== FIN =====")
close(LOG)

# Copia las figuras (y el logo, si falta) a informe/ para que esa carpeta compile sola
mkpath("informe")
for f in readdir("figures")
    cp(joinpath("figures", f), joinpath("informe", f), force = true)
end
isfile("informe/logo_udesa.pdf") || cp("../PS2/informe/logo_udesa.pdf", "informe/logo_udesa.pdf")
