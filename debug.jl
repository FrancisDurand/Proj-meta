include("functions.jl")


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

begin
    Random.seed!(0)
    ENV["JULIA_DEBUG"] = Main
    instance = instance_list[3]

    # DO
end
