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

# Définit les instances à évaluer
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
    read_instance("graphs/dsjc1000.5.col", 84)
)

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



# Algorithme génétique

# On commence par coder le tabou et si on a le temps on fait le génétique
# ça a l'air de prendre beaucoup plus de temps que je croyais

"""Renvoie une solution aléatoire selon la règle de la roue de la fortune."""
function roue_fortune(instance::Instance,population::Vector{Solution})
    # Solutions supposées évaluées
    weights = [1/sample.obj for sample ∈ population]
    cumw = cumsum(weights)
    r = rand(1:sum(weights))
    return findfirst(r .<= cumw)
end

function genetique(μ::Int)
    
end
