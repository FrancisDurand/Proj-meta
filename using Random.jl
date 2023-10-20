using Random
cd("c:\\Users\\franc\\Proj-meta")
pwd()
io = open("graphs/dsjc125.1.col") 

function recup_numbre_dans_string(str::String,i::Int) #cette fonction recupere le nombre qui commence à l'indice i dans le string str
    j = i
    while length(str) >= j+1 && str[j+1] >= '0' && str[j+1] <= '9' #verifie si on manipule bien un chiffre
        j+=1
    end
    return (parse(Int64,str[i:j]),j)
end

b = true
graphe = BitMatrix
while (b)
    ligne = readline(io)
    if length(ligne) >= 1 && ligne[1] == 'p'
        n = recup_numbre_dans_string(ligne,8)[1] #nombre de sommets dans le graphe
        graphe = falses(n,n)
        display(n)
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


display(graphe)

function is_collision(graphe::BitMatrix,solution::Vector{Int},i::Int,j::Int)
    return graphe[i,j] && solution[i] == solution[j]
end

function nbr_collision(graphe::BitMatrix,solution::Vector{Int})
    conteur = 0 
    for i = 1:size(graphe)[1]
        for j = i+1:size(graphe)[2]
            conteur+= is_collision(graphe,solution,i,j)
        end
    end
end

"""
Calcule le nombre de collisions au sommet position en attribuant la couleur 
"""
function nb_collision_sommet_couleur(graphe::BitMatrix,solution::Vector{Int},position::Int,couleur::Int) 
    n_collision = 0
    for  i = 1:length(solution)
        n_collision += graphe[position,i] && solution[i] == couleur
    end
    return n_collision
end

solution = convert(Vector{Int},bitrand(125))

nb_collision_sommet_couleur(graphe,solution,9,1)


A =  [true false; true true]

sum(A,dims=1)

function calcul_de_ordre_des_sommets(graphe::BitMatrix) 
    deg = vec(sum(graphe,dims=1))
    return sortperm(deg,rev=true)
end

function meilleure_couleur_locale(graphe::BitMatrix,k,solution::Vector{Int},position)
    best_c = 1
    best_nb_collision = length(solution)
    for c = 1:k
        nb_collision = nb_collision_sommet_couleur(graphe,solution,position,c)
        if nb_collision < best_nb_collision
            best_c = c
            best_nb_collision = nb_collision
        end
    end
    return best_c
end



meilleure_couleur_locale(graphe,2,solution,9)

function glouton(graphe::BitMatrix,k::Int)
    solution = zeros(Int,size(graphe)[1]) # la couleur 0 veut dire non colorée
    ordre_des_sommets = calcul_de_ordre_des_sommets(graphe)
    for i in ordre_des_sommets
        c = meilleure_couleur_locale(graphe,k,solution,i)
        solution[i] = c
    end
    return solution
end

glouton(graphe,14)