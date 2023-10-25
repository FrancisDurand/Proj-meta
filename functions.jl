import LinearAlgebra
using Logging
using Random



# Définition des structures de données

struct Instance
    name::String
    graphe::BitMatrix # Matrice d'adjacence du graphe
    k::Int
end

mutable struct Solution
    nodecolors::Vector{Int} # coloration
    conflictmatrix::Matrix{Int} # Matrice de conflits si l'on assigne la couleur i au nœud j
    tabuexpiry::Matrix{Int} # Matrice des itérations minimum où une certaine modification est autorisée
    obj::Int # Nombre de conflits
end

Solution(instance::Instance,nodecolors::Vector{Int}) = Solution(
    nodecolors,
    zeros(Int,(instance.k,length(instance))),
    zeros(Int,(instance.k,length(instance))),
    length(instance)^2
)

"""Renvoie le nombre de sommets du graphe de l'instance."""
function Base.length(instance::Instance)
    return size(instance.graphe)[1]
end

function Base.length(solution::Solution)
    return length(solution.nodecolors)
end

# Chargement des données

function recup_numbre_dans_string(str::String,i::Int) #cette fonction recupere le nombre qui commence à l'indice i dans le string str
    j = i
    while length(str) >= j+1 && str[j+1] >= '0' && str[j+1] <= '9' #verifie si on manipule bien un chiffre
        j+=1
    end
    return (parse(Int64,str[i:j]),j)
end

"""Construit une instance de k-coloration à partir du fichier."""
function read_instance(path::String,k::Int)
    io = open(path)
    b = true
    graphe = BitMatrix
    while (b)
        ligne = readline(io)
        if length(ligne) >= 1 && ligne[1] == 'p'
            n = recup_numbre_dans_string(ligne,8)[1] #nombre de sommets dans le graphe
            graphe = falses(n,n)
        end
        if length(ligne) >= 1 && ligne[1] == 'e'
            u,j = recup_numbre_dans_string(ligne,3)
            v = recup_numbre_dans_string(ligne,j+1)[1]
            graphe[u,v] = true
            graphe[v,u] = true
        end
        if length(ligne) < 1
            b = false
        end
    end
    name = split(path,"/")[2]
    return Instance(name,graphe,k)
end


# Définition des instances de base

@static if !@isdefined(instance_list)
    const instance_list = (
        read_instance("graphs/flat300_26_0.col", 26),
        read_instance("graphs/le450_15c.col", 15),
        read_instance("graphs/dsjc125.1.col", 5),
        read_instance("graphs/dsjc125.9.col", 44),
        read_instance("graphs/dsjc250.1.col", 8),
        read_instance("graphs/dsjc250.9.col", 72),
        read_instance("graphs/dsjc250.5.col", 28),
        read_instance("graphs/dsjc1000.5.col", 86),
        read_instance("graphs/dsjc1000.5.col", 85),
        read_instance("graphs/dsjc1000.5.col", 84))
end


# Fonctions de base

"""Calcule le nombre de collisions de la coloration."""
function nbr_collision(instance::Instance,solution::Solution)
    compteur = 0 
    for j = 1:length(instance)
        for i = j+1:length(instance)
            compteur += instance.graphe[i,j] && solution.nodecolors[i] == solution.nodecolors[j]
        end
    end
    solution.obj = compteur
    return compteur
end

"""Calcule le nombre de collisions au sommet position en lui attribuant couleur."""
function nb_collision_sommet_couleur(instance::Instance,solution::Solution,position::Int,couleur::Int)
    n_collision = 0
    for  i = 1:length(solution)
        n_collision += instance.graphe[i,position] && solution.nodecolors[i] == couleur
    end
    return n_collision
end

"""Calcule l'ordre dans lequel traiter les sommets dans l'heuristique gloutonne (degré décroissant)."""
function calcul_de_ordre_des_sommets(instance::Instance) 
    deg = vec(sum(instance.graphe,dims=1))
    return sortperm(deg,rev=true)
end

"""Détermine la meilleure couleur à attribuer au nœud."""
function meilleure_couleur_locale(instance::Instance,solution::Solution,position::Int)
    best_c = 1
    best_nb_collision = length(solution)
    for c = 2:instance.k # Les couleurs sont entre 1 et k
        nb_collision = nb_collision_sommet_couleur(instance,solution,position,c)
        if nb_collision < best_nb_collision
            best_c = c
            best_nb_collision = nb_collision
        end
    end
    return best_c
end



# Algorithme glouton

"""Heuristique gloutonne générant une solution par coloration successive des sommets."""
function glouton(instance::Instance)
    solution = Solution(instance,zeros(Int,length(instance))) # la couleur 0 veut dire non colorée
    ordre_des_sommets = calcul_de_ordre_des_sommets(instance)
    for i ∈ ordre_des_sommets
        c = meilleure_couleur_locale(instance,solution,i)
        solution.nodecolors[i] = c
    end
    return solution
end

"""Génère une solution aléatoire."""
function sol_alea(instance::Instance)
    return Solution(instance,rand(1:instance.k,length(instance)))
end



# Opérateur de voisinage

function fill_conflicts(instance::Instance,solution::Solution)
    nbr_collision(instance,solution)
    for i = 1:length(instance)
        for color = 1:instance.k
            solution.conflictmatrix[color,i] = nb_collision_sommet_couleur(instance,solution,i,color) 
        end
    end
end

function update_conflicts(instance::Instance,solution::Solution,node::Int,color::Int)
    for othernode = 1:length(instance)
        if instance.graphe[othernode,node]
            for othercolor in 1:instance.k
                solution.conflictmatrix[othercolor,othernode] += (othercolor==color) - (othercolor==solution.nodecolors[node])
            end
        end
    end
    solution.obj += solution.conflictmatrix[color,node] - solution.conflictmatrix[solution.nodecolors[node],node]
    solution.nodecolors[node] = color
end

"""Opérateur appliquant une couleur à un nœud de manière à minimiser les conflits."""
function simple_neighbor(instance::Instance,solution::Solution)
    if iszero(solution.conflictmatrix)
        fill_conflicts(instance,solution)
    end
    improvementmatrix = solution.conflictmatrix .-
        sum(solution.conflictmatrix .* (1:instance.k .== permutedims(solution.nodecolors)),dims=1)
    color,node = Tuple(argmin(improvementmatrix))
    update_conflicts(instance,solution,node,color)
end


# Tabou

"""Version de `simple_neighbor` qui implémente le tabou."""
function tabu_neighbor(instance::Instance,solution::Solution,tabulength::Int,it::Int)
    if iszero(solution.conflictmatrix)
        fill_conflicts(instance,solution)
    end
    improvementmatrix = solution.conflictmatrix .-
        sum(solution.conflictmatrix .* (1:instance.k .== permutedims(solution.nodecolors)),dims=1)
    tabupenalizedimprovementmatrix = improvementmatrix + 2*length(instance)*(
        solution.tabuexpiry .> it) # TABOU (application)
    color,node = Tuple(argmin(tabupenalizedimprovementmatrix))
    solution.tabuexpiry[color,node] = it+tabulength # TABOU (ajout)
    update_conflicts(instance,solution,node,color)
end



# Algorithme mémétique

"""Lance une recherche locale simple sur chaque individu de la population. S'arrête dès que l'on atteint un minimum local ou un plateau."""
function simple_local_search(instance::Instance,population::Vector{Solution},MAXIT::Int)
    for solution ∈ population
        old_obj = length(instance)^2
        for t = 1:MAXIT
            simple_neighbor(instance,solution)
            # @debug string("it: ",t,"\t","conflits: ",solution.obj)
            if old_obj==solution.obj
                break
            end
            old_obj=solution.obj
        end
    end
end

"""Renvoie un indice aléatoire selon la règle de la roue de la fortune."""
function roue_fortune(weights::Vector{Float64})
    cumw = cumsum(weights)
    r = rand()*sum(weights)
    return findfirst(r .<= cumw)
end

"""Calcule la permutation minimisant la distance entre deux solutions.\\
On cherche donc la permutation maximisant ∑S(i,perm(i)) avec S la matrice de similarité."""
function permutation_distance(similarity::Matrix{Int})
    # Le problème de la recherche de la meilleure permutation est NP-complet
    # Donc on va pas s'amuser à écrire un algo exact
    # Ceci est donc une heuristique
    perm = zeros(Int,size(similarity)[1])
    m,AB = findmax(similarity) # On choisit la plus grande valeur de S -> perm[b] ≟ a
    M = m+1
    a,b = Tuple(AB)
    while 0∈perm
        similarity[a,b] -= M
        m2,u = findmax(@view similarity[:,b]) # On cherche la deuxième meilleure valeur de S sur la colonne b -> perm[b] ≟ u
        m3,v = findmax(@view similarity[a,:])                                                       # ligne a -> perm[v] ≟ a
        if m2 + m3 ≤ m + similarity[u,v]
            # On se dit qu'il y a des chances que perm[b] = a si S(a,b)+S(u,v)≥S(u,b)+S(a,v), sinon on pourrait échanger et obtenir perm[b] = u et perm[v] = a qui aurait un meilleur objectif
            @view(similarity[:,b]) .-= M
            @view(similarity[a,:]) .-= M
            perm[b] = a
            m,AB = findmax(similarity) # Il y a des chances que les plus grands éléments de S soient atteints par la permutation
            a,b = Tuple(AB)
        else
            # Il y a des chances que perm[b] ≠ a donc il y a des chances qu'il y ait au moins perm[b] = u ou perm[v] = a (il existe un unique élément de S atteint sur chaque ligne et chaque colonne). On commence par tester celui qui atteint la plus grande valeur de S
            if m2≥m3
                m = m2
                a = u
            else
                m = m3
                b = v
            end
        end
    end
    return perm
end

"""Calcule la distance entre deux solutions."""
function distance(instance::Instance,x::Solution,y::Solution)
    similarity = zeros(Int,(instance.k,instance.k))
    for i = 1:length(instance)
        similarity[x.nodecolors[i],y.nodecolors[i]] += 1
    end
    perm = permutation_distance(similarity)
    dist = 0
    for i = 1:length(instance)
        if x.nodecolors[i] == perm[y.nodecolors[i]]
            dist += 1
        end
    end
    return dist
end

"""Copie des classes de couleurs entières des parents dans l'enfant\\
(ce qu'il est possible de copier, on commence par les plus grandes classes)."""
function croisement(instance::Instance,population::Vector{Solution},mom::Int,dad::Int)
    enfant = Solution(instance,zeros(Int,length(instance)))
    parentnodecolors = [population[mom].nodecolors population[dad].nodecolors]
    colorcount = permutedims(stack([vec(sum(parentnodecolors.==i,dims=1)) for i = 1:instance.k]))
    while 0 ∈ enfant.nodecolors
        color,parent = Tuple(argmax(colorcount))
        colorcount[color,parent] = 0
        enfant.nodecolors += color*(parentnodecolors[:,parent] .== color .&& enfant.nodecolors .== 0)
    end
    return enfant
end

"""Génère λ enfants à partir des parents et rejette ceux trop proches d'une solution existante."""
function faire_enfants(instance::Instance,population::Vector{Solution},λ::Int,DISTTHR::Int,BESTOBJ::Int)
    enfants = empty(population)
    for i = 1:λ
        mom=dad=0
        while mom==dad
            weights = [1/sample.obj for sample ∈ population]
            mom,dad = roue_fortune(weights),roue_fortune(weights)
        end
        e = croisement(instance,population,mom,dad)
        # Rejeter si trop proche de solutions existantes
        dmin = minimum(distance(instance,sample,e) for sample ∈ population)
        if dmin ≥ DISTTHR || nbr_collision(instance,e) < BESTOBJ
            push!(enfants,e)
        end
    end
    return enfants
end

"""Opérateur d'élimination, évite la convergence prématurée en maintenant une diversité dans la population."""
function discard_excess(instance::Instance,population::Vector{Solution},DISTTHR::Int,popsize::Int)
    while length(population)>popsize
        distancetriangle = Matrix{Int}(undef,length(population),length(population))
        for j = 1:length(population)
            for i = j:length(population)
                distancetriangle[i,j] = distance(instance,population[i],population[j])
            end
        end
        distancetriangle += LinearAlgebra.UniformScaling(length(instance))
        distancematrix = LinearAlgebra.Symmetric(distancetriangle)
        dmin,IJ = findmin(distancematrix)
        I,J = 0,0
        if dmin < DISTTHR
            # Élimine les solutions trop proches
            I,J = Tuple(IJ)
        else
            # Maintient une distance moyenne entre les solutions
            weights = convert(Vector{Float64},[sample.obj for sample ∈ population])
            I = roue_fortune(weights)
            J = argmin(distancematrix[:,I])
        end
        if population[I].obj ≥ population[J].obj
            popat!(population,I)
        else
            popat!(population,J)
        end
    end
end

"""Algorithme mémétique"""
function genetique(instance::Instance,popsize::Int,nbchildren::Int,DISTTHR::Int,MAXIT::Int,localMAXIT::Int)
    start_time = time()
    population = [sol_alea(instance) for i = 1:popsize]
    # Remplacer par le tabou ?
    simple_local_search(instance,population,localMAXIT)
    BESTOBJ = minimum(sample.obj for sample ∈ population)
    @info string("it: ",0,"\ttemps:",trunc(100*(time()-start_time))/100,"s\tconflits: ",[solution.obj for solution ∈ population])
    for t = 1:MAXIT
        enfants = faire_enfants(instance,population,nbchildren,DISTTHR,BESTOBJ)
        # Remplacer par le tabou ?
        simple_local_search(instance,enfants,localMAXIT)
        BESTOBJ = isempty(enfants) ? BESTOBJ : min(minimum(sample.obj for sample ∈ enfants),BESTOBJ)
        append!(population,enfants)
        discard_excess(instance,population,DISTTHR,popsize)
        @info string("it: ",t,"\ttemps:",trunc(100*(time()-start_time))/100,"s\tconflits: ",[solution.obj for solution ∈ population])
    end
    return population[argmin(solution.obj for solution ∈ population)]
end
